# Post-freeze correlation-aware sensitivity for the locked holdout panel.
# This is supplementary only: it does not alter the frozen permutation test,
# candidates, tiers, or holdout rules.
suppressPackageStartupMessages({library(edgeR); library(limma); library(jsonlite)})

args <- commandArgs(trailingOnly=FALSE)
script_arg <- args[grep("^--file=", args)]
if (length(script_arg) == 0) stop("Run with Rscript postfreeze_correlation_aware_gene_set.R")
script_path <- sub("^--file=", "", script_arg[[1]])
ROOT <- normalizePath(file.path(dirname(script_path), ".."), winslash="/", mustWork=TRUE)
OUT <- file.path(ROOT, "13_postfreeze_interpretation")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

raw_file <- file.path(ROOT, "01_data_raw/holdout/E-MTAB-16655/RawCounts.csv")
norm_file <- file.path(ROOT, "01_data_raw/holdout/E-MTAB-16655/FilteredNormLogCPM.csv")
sdrf_file <- file.path(ROOT, "02_metadata/E-MTAB-16655.sdrf.txt")
panel_file <- file.path(ROOT, "08_evidence_atlas/FROZEN_CANDIDATE_PANEL_PRE_HOLDOUT.tsv")
for (p in c(raw_file, norm_file, sdrf_file, panel_file)) if (!file.exists(p)) stop("Missing input: ", p)

raw <- read.csv(raw_file, check.names=FALSE, stringsAsFactors=FALSE)
ann <- read.csv(norm_file, check.names=FALSE, stringsAsFactors=FALSE)
sdrf <- read.delim(sdrf_file, check.names=FALSE, stringsAsFactors=FALSE)
sample_design <- unique(sdrf[, c("Source Name", "Factor Value[stimulus]")])
rownames(sample_design) <- sample_design[["Source Name"]]
sample_cols <- intersect(colnames(raw), rownames(sample_design))
if (length(sample_cols) != 45) stop("Expected 45 samples; observed ", length(sample_cols))

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
contrast <- makeContrasts(WT-UI, levels=design)

y <- DGEList(counts=counts)
keep <- filterByExpr(y, design=design)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y, method="TMM")
v <- voom(y, design, plot=FALSE)

panel <- read.delim(panel_file, check.names=FALSE, stringsAsFactors=FALSE)
panel <- panel[panel$human_symbol %in% rownames(v$E) & panel$expected_wt_vs_uninfected %in% c("UP", "DOWN"), ]
panel$orientation <- ifelse(panel$expected_wt_vs_uninfected == "UP", 1, -1)

scenarios <- list(
  full_54_gene_panel=character(0),
  minus_RAB1A="RAB1A",
  minus_RAB10="RAB10",
  minus_RAB1A_RAB10=c("RAB1A", "RAB10")
)
out <- list()
for (nm in names(scenarios)) {
  sub <- panel[!(panel$human_symbol %in% scenarios[[nm]]), ]
  idx <- match(sub$human_symbol, rownames(v$E))
  E <- v$E
  E[idx, ] <- E[idx, ] * sub$orientation
  cam <- camera(E, index=list(panel=idx), design=design, contrast=contrast, sort=FALSE)
  removed_label <- if (length(scenarios[[nm]]) == 0) "NONE" else paste(sort(scenarios[[nm]]), collapse=";")
  cam_direction <- if ("Direction" %in% names(cam)) as.character(cam$Direction[1]) else NA_character_
  cam_cor <- if ("Correlation" %in% names(cam)) as.numeric(cam$Correlation[1]) else NA_real_
  cam_p <- if ("PValue" %in% names(cam)) as.numeric(cam$PValue[1]) else NA_real_
  out[[nm]] <- data.frame(
    scenario=nm,
    removed_genes=removed_label,
    panel_n=nrow(sub),
    camera_direction=cam_direction,
    inter_gene_cor=cam_cor,
    camera_p=cam_p,
    method="limma::camera on edgeR-voom expression; expected directions encoded by row orientation",
    stringsAsFactors=FALSE
  )
}
out <- do.call(rbind, out)
write.table(out, file.path(OUT, "g_camera_gene_set_sensitivity.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
writeLines(toJSON(list(
  method="correlation-aware competitive gene-set sensitivity",
  model="edgeR TMM + voom; design ~ 0 + group + replicate_batch",
  contrast="WT_vs_UI",
  expected_direction="panel rows multiplied by +1 for UP and -1 for DOWN before camera",
  primary_test_unchanged=TRUE,
  results=out
), pretty=TRUE, auto_unbox=TRUE), file.path(OUT, "g_camera_gene_set_sensitivity.json"))
print(out)
