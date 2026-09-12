suppressPackageStartupMessages({library(babelgene); library(jsonlite); library(readxl)})

args <- commandArgs(trailingOnly=FALSE)
script_arg <- args[grep("^--file=", args)]
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."), winslash="/", mustWork=TRUE)
raw_dir <- file.path(root, "01_data_raw", "lcv_proteomics")
human_dir <- file.path(raw_dir, "PMC4286187_supplementary")
out_dir <- file.path(root, "08_evidence_atlas")
cp_dir <- file.path(root, "checkpoints")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(cp_dir, recursive=TRUE, showWarnings=FALSE)

human_wt_file <- file.path(human_dir, "pr500765x_si_001.xlsx")
human_ankb_file <- file.path(human_dir, "pr500765x_si_002.xlsx")
mouse_file <- file.path(raw_dir, "PLOS_S2_LCV_untreated.xlsx")
required <- c(human_wt_file, human_ankb_file, mouse_file)
if (any(!file.exists(required))) stop("Missing LCV proteomics input: ", paste(required[!file.exists(required)], collapse=", "))

read_human <- function(path, condition) {
  z <- as.data.frame(read_excel(path, sheet=1, skip=2, col_names=FALSE))
  z <- z[rowSums(!is.na(z)) > 0, , drop=FALSE]
  gene <- trimws(as.character(z[[1]])); protein <- trimws(as.character(z[[2]]))
  accession <- apply(z, 1, function(r) {
    v <- trimws(as.character(r)); v <- v[!is.na(v) & nzchar(v)]
    if (length(v) >= 3) tail(v, 1) else NA_character_
  })
  keep <- !is.na(gene) & nzchar(gene) & toupper(gene) != "GENE"
  data.frame(
    source_record_id=paste0("HUMAN_", condition, "_", seq_len(sum(keep))),
    human_symbol=toupper(gene[keep]), original_symbol=gene[keep],
    protein_name=protein[keep], accession=accession[keep], host_species="HUMAN",
    condition=condition, evidence_type="LCV_SPATIAL_PROTEOMICS",
    quantitative_value=NA_real_, quantitative_unit=NA_character_,
    source_pmid="25369898", source_pmcid="PMC4286187", stringsAsFactors=FALSE)
}

human_wt <- read_human(human_wt_file, "WT_LCV_U937")
human_ankb <- read_human(human_ankb_file, "ANKB_MUTANT_LCV_U937")

mouse_raw <- as.data.frame(read_excel(mouse_file, sheet="identified proteins (host)"))
gene_col <- "Gene names"; lfq_col <- "Average LFQ intensity in untreated samples"
mouse_expanded <- do.call(rbind, lapply(seq_len(nrow(mouse_raw)), function(i) {
  symbols <- trimws(unlist(strsplit(as.character(mouse_raw[[gene_col]][i]), ";", fixed=TRUE)))
  symbols <- symbols[nzchar(symbols) & !is.na(symbols)]
  if (!length(symbols)) return(NULL)
  data.frame(original_symbol=symbols, protein_name=as.character(mouse_raw[["Protein name"]][i]),
    accession=as.character(mouse_raw[["Uniprot identifier"]][i]),
    quantitative_value=as.numeric(mouse_raw[[lfq_col]][i]), stringsAsFactors=FALSE)
}))
mouse_symbols <- unique(mouse_expanded$original_symbol)
ortho <- orthologs(mouse_symbols, species="mouse", human=FALSE, min_support=3, top=TRUE)
omap <- setNames(toupper(ortho$human_symbol), ortho$symbol)
mouse_expanded$human_symbol <- unname(omap[mouse_expanded$original_symbol])
mouse_expanded <- mouse_expanded[!is.na(mouse_expanded$human_symbol) & nzchar(mouse_expanded$human_symbol), , drop=FALSE]
mouse <- data.frame(
  source_record_id=paste0("MOUSE_UNTREATED_", seq_len(nrow(mouse_expanded))),
  human_symbol=mouse_expanded$human_symbol, original_symbol=mouse_expanded$original_symbol,
  protein_name=mouse_expanded$protein_name, accession=mouse_expanded$accession,
  host_species="MOUSE", condition="UNTREATED_LCV_BMM", evidence_type="LCV_SPATIAL_PROTEOMICS",
  quantitative_value=mouse_expanded$quantitative_value, quantitative_unit="AVERAGE_LFQ_INTENSITY",
  source_pmid="26829557", source_pmcid="PMC4734697", stringsAsFactors=FALSE)

records <- rbind(human_wt, human_ankb, mouse)
records <- records[order(records$human_symbol, records$host_species, records$condition), ]
write.table(records, file.path(out_dir, "lcv_spatial_source_records.tsv"), sep="\t", row.names=FALSE, quote=FALSE, na="")

genes <- sort(unique(records$human_symbol))
evidence <- do.call(rbind, lapply(genes, function(g) {
  d <- records[records$human_symbol == g, , drop=FALSE]
  data.frame(human_symbol=g, lcv_human_wt=any(d$condition == "WT_LCV_U937"),
    lcv_human_ankb=any(d$condition == "ANKB_MUTANT_LCV_U937"),
    lcv_mouse_untreated=any(d$condition == "UNTREATED_LCV_BMM"), lcv_spatial_supported=TRUE,
    lcv_species_n=length(unique(d$host_species)), lcv_source_records=nrow(d),
    lcv_sources=paste(sort(unique(paste(d$source_pmid, d$condition, sep=":"))), collapse=";"), stringsAsFactors=FALSE)
}))
write.table(evidence, file.path(out_dir, "lcv_spatial_evidence.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

truth <- read.delim(file.path(root, "00_protocol", "KNOWN_POSITIVE_CONTROL_HOST_FACTORS.tsv"), stringsAsFactors=FALSE, check.names=FALSE)
truth$human_symbol <- toupper(truth$gene_symbol)
truth_lcv <- merge(truth, evidence, by="human_symbol", all.x=TRUE)
for (cc in c("lcv_human_wt", "lcv_human_ankb", "lcv_mouse_untreated", "lcv_spatial_supported")) truth_lcv[[cc]][is.na(truth_lcv[[cc]])] <- FALSE
write.table(truth_lcv, file.path(out_dir, "known_positive_lcv_recovery.tsv"), sep="\t", row.names=FALSE, quote=FALSE, na="")

stats <- list(human_wt_records=nrow(human_wt), human_ankb_records=nrow(human_ankb),
  mouse_records_raw=nrow(mouse_raw), mouse_gene_records_mapped=nrow(mouse),
  mouse_symbols_submitted=length(mouse_symbols), mouse_orthologs_mapped=nrow(ortho),
  atlas_spatial_genes=nrow(evidence), spatial_genes_human=sum(evidence$lcv_human_wt | evidence$lcv_human_ankb),
  spatial_genes_mouse=sum(evidence$lcv_mouse_untreated), truth_set_spatial_recovered=sum(truth_lcv$lcv_spatial_supported),
  gate6b_lcv_spatial=if (nrow(human_wt) > 1000 && nrow(mouse) > 500 && nrow(evidence) > 1000) "PASS" else "FAIL",
  evidence_boundary="LCV detection is spatial co-localization evidence, not proof of direct effector binding or functional dependency",
  holdout_accessed=FALSE)
writeLines(toJSON(stats, auto_unbox=TRUE, pretty=TRUE), file.path(out_dir, "stage6b_statistics.json"))
report <- c("# Stage 6B LCV spatial proteomics evidence", "",
  paste0("- Human WT U937 LCV records: ", stats$human_wt_records),
  paste0("- Human ankB-mutant U937 LCV records: ", stats$human_ankb_records),
  paste0("- Mouse LCV human-ortholog records: ", stats$mouse_gene_records_mapped),
  paste0("- Unique human-symbol spatial evidence genes: ", stats$atlas_spatial_genes),
  paste0("- Known-positive controls recovered by LCV evidence: ", stats$truth_set_spatial_recovered, "/", nrow(truth_lcv)),
  paste0("- Gate: ", stats$gate6b_lcv_spatial), "",
  "LCV detection is treated as orthogonal spatial co-localization evidence. It is not interpreted as direct physical binding or infection-specific functional dependency. Mouse symbols were mapped to supported human orthologs with babelgene.")
writeLines(report, file.path(out_dir, "STAGE6B_LCV_SPATIAL_PROTEOMICS.md"))
writeLines(toJSON(stats, auto_unbox=TRUE, pretty=TRUE), file.path(cp_dir, "stage6b_lcv_proteomics.done"))
cat(toJSON(stats, auto_unbox=TRUE, pretty=TRUE), "\n")
if (stats$gate6b_lcv_spatial != "PASS") quit(status=2)
