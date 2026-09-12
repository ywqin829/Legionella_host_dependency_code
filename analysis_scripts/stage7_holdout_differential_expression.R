suppressPackageStartupMessages({library(edgeR); library(jsonlite)})

args <- commandArgs(trailingOnly=FALSE)
script_arg <- args[grep("^--file=", args)]
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."), winslash="/", mustWork=TRUE)
raw_file <- file.path(root, "01_data_raw", "holdout", "E-MTAB-16655", "RawCounts.csv")
norm_file <- file.path(root, "01_data_raw", "holdout", "E-MTAB-16655", "FilteredNormLogCPM.csv")
sdrf_file <- file.path(root, "02_metadata", "E-MTAB-16655.sdrf.txt")
freeze_file <- file.path(root, "08_evidence_atlas", "CANDIDATE_FREEZE_MANIFEST.json")
out_dir <- file.path(root, "09_holdout_validation")
cp_dir <- file.path(root, "checkpoints")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(cp_dir, recursive=TRUE, showWarnings=FALSE)
for (p in c(raw_file, norm_file, sdrf_file, freeze_file)) if (!file.exists(p)) stop("Missing input: ", p)
freeze <- fromJSON(freeze_file)
if (!isTRUE(freeze$holdout_expression_unlock_authorized)) stop("Holdout expression is not authorized by the pre-holdout freeze")

raw <- read.csv(raw_file, check.names=FALSE, stringsAsFactors=FALSE)
ann <- read.csv(norm_file, check.names=FALSE, stringsAsFactors=FALSE)
sdrf <- read.delim(sdrf_file, check.names=FALSE, stringsAsFactors=FALSE)
sample_design <- unique(sdrf[, c("Source Name", "Factor Value[stimulus]")])
rownames(sample_design) <- sample_design[["Source Name"]]

sample_cols <- intersect(colnames(raw), rownames(sample_design))
if (length(sample_cols) != 45) stop("Expected 45 holdout samples; observed ", length(sample_cols))
counts <- as.matrix(raw[, sample_cols, drop=FALSE])
storage.mode(counts) <- "integer"
rownames(counts) <- sub("\\.[0-9]+$", "", raw$Geneid)

ann$ensembl_clean <- sub("\\.[0-9]+$", "", ann$ensembl_gene_id)
gene_map <- unique(ann[, c("ensembl_clean", "external_gene_name")])
gene_map <- gene_map[!is.na(gene_map$external_gene_name) & nzchar(gene_map$external_gene_name), ]
gene_map <- gene_map[!duplicated(gene_map$ensembl_clean), ]
map <- setNames(toupper(gene_map$external_gene_name), gene_map$ensembl_clean)
symbols <- unname(map[rownames(counts)])
keep_mapped <- !is.na(symbols) & nzchar(symbols)
counts <- counts[keep_mapped, , drop=FALSE]
symbols <- symbols[keep_mapped]
counts <- rowsum(counts, group=symbols, reorder=FALSE)

group <- sub("[234]$", "", sample_cols)
group <- factor(group, levels=c("UI", "WT", "dotA", "A", "B", "C", "D", "E", "F", "FGHI", "I", "O", "P", "R", "S"))
batch <- factor(sub("^.*([234])$", "\\1", sample_cols))
design <- model.matrix(~0 + group + batch)
colnames(design) <- sub("^group", "", colnames(design))

y <- DGEList(counts=counts)
keep <- filterByExpr(y, design=design)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y, method="TMM")
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y, design, robust=TRUE)

contrast_defs <- c(WT_vs_UI="WT-UI", WT_vs_dotA="WT-dotA",
  A_vs_WT="A-WT", B_vs_WT="B-WT", C_vs_WT="C-WT", D_vs_WT="D-WT",
  E_vs_WT="E-WT", F_vs_WT="F-WT", FGHI_vs_WT="FGHI-WT", I_vs_WT="I-WT",
  O_vs_WT="O-WT", P_vs_WT="P-WT", R_vs_WT="R-WT", S_vs_WT="S-WT")

all_results <- list(); summaries <- list()
for (nm in names(contrast_defs)) {
  cm <- makeContrasts(contrasts=contrast_defs[[nm]], levels=design)
  qlf <- glmQLFTest(fit, contrast=cm)
  tab <- topTags(qlf, n=Inf, sort.by="none")$table
  tab$human_symbol <- rownames(tab); tab$contrast_id <- nm
  tab <- tab[, c("contrast_id", "human_symbol", "logFC", "logCPM", "F", "PValue", "FDR")]
  all_results[[nm]] <- tab
  summaries[[nm]] <- data.frame(contrast_id=nm, numerator=sub("_vs_.*$", "", nm),
    denominator=sub("^.*_vs_", "", nm), genes_tested=nrow(tab),
    fdr05=sum(tab$FDR < .05), fdr05_abs_logfc1=sum(tab$FDR < .05 & abs(tab$logFC) >= 1),
    stringsAsFactors=FALSE)
}
results <- do.call(rbind, all_results)
summary <- do.call(rbind, summaries)
write.table(results, gzfile(file.path(out_dir, "holdout_all_contrasts.tsv.gz")), sep="\t", row.names=FALSE, quote=FALSE)
write.table(summary, file.path(out_dir, "holdout_contrast_summary.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(data.frame(sample=sample_cols, group=group, batch=batch, library_size=colSums(counts)),
  file.path(out_dir, "holdout_sample_qc.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

stats <- list(dataset="E-MTAB-16655/PRJEB108343", samples=length(sample_cols), groups=length(levels(group)),
  mapped_symbols=nrow(counts), retained_symbols=sum(keep), contrasts=length(contrast_defs),
  design="~ 0 + group + replicate_batch", method="edgeR quasi-likelihood, TMM, robust dispersion",
  raw_counts_used=TRUE, author_logCPM_used_for_annotation_only=TRUE,
  candidate_freeze_sha256=freeze$candidate_panel_sha256, holdout_unlocked_after_freeze=TRUE,
  gate7_differential_expression=if (length(sample_cols)==45 && length(contrast_defs)==14 && sum(keep)>10000) "PASS" else "FAIL")
writeLines(toJSON(stats, auto_unbox=TRUE, pretty=TRUE), file.path(out_dir, "stage7_statistics.json"))
writeLines(c("# Stage 7 external holdout differential expression", "",
  paste0("- Samples: ", stats$samples, " across ", stats$groups, " groups"),
  paste0("- Retained mapped genes: ", stats$retained_symbols),
  paste0("- Pre-specified contrasts: ", stats$contrasts),
  paste0("- Gate: ", stats$gate7_differential_expression), "",
  "The holdout was opened only after the candidate panel and validation rules were frozen. Raw counts were modelled with edgeR quasi-likelihood; the author's normalized matrix supplied gene annotation only. Distinct multi-effector deletion groups were not pooled."),
  file.path(out_dir, "STAGE7_HOLDOUT_DIFFERENTIAL_EXPRESSION.md"))
writeLines(toJSON(stats, auto_unbox=TRUE, pretty=TRUE), file.path(cp_dir, "stage7_holdout_de.done"))
cat(toJSON(stats, auto_unbox=TRUE, pretty=TRUE), "\n")
if (stats$gate7_differential_expression != "PASS") quit(status=2)

