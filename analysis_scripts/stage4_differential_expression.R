suppressPackageStartupMessages({library(limma); library(edgeR); library(jsonlite)})

args <- commandArgs(trailingOnly=FALSE)
script_arg <- args[grep("^--file=", args)]
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."), winslash="/", mustWork=TRUE)
matrix_dir <- file.path(root, "05_expression", "matrices")
map_dir <- file.path(root, "05_expression", "gene_maps")
out_dir <- file.path(root, "06_differential_expression")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(root, "checkpoints"), recursive=TRUE, showWarnings=FALSE)
design_all <- read.delim(file.path(root, "02_metadata", "stage3_frozen_sample_design.tsv"), check.names=FALSE, stringsAsFactors=FALSE)

read_mat <- function(path) {
  x <- read.delim(gzfile(path), check.names=FALSE, stringsAsFactors=FALSE)
  rn <- as.character(x[[1]]); x <- as.matrix(x[,-1,drop=FALSE]); mode(x) <- "numeric"; rownames(x) <- rn; x
}

collapse_gene <- function(dataset, mat, counts=FALSE) {
  mp <- read.delim(file.path(map_dir, paste0(dataset, "_gene_map.tsv")), stringsAsFactors=FALSE, check.names=FALSE)
  idx <- match(rownames(mat), mp$feature_id)
  symbol <- mp$gene_symbol[idx]; ambiguous <- mp$ambiguous_mapping[idx]
  keep <- !is.na(symbol) & nzchar(symbol) & !is.na(ambiguous) & !ambiguous
  mat <- mat[keep,,drop=FALSE]; symbol <- symbol[keep]
  summed <- rowsum(mat, group=symbol, reorder=FALSE, na.rm=TRUE)
  if (counts) return(round(summed))
  denom <- as.numeric(table(factor(symbol, levels=rownames(summed))))
  summed / denom
}

write_result <- function(tab, dataset, contrast, modality, note="") {
  if (!"AveExpr" %in% names(tab) && "logCPM" %in% names(tab)) tab$AveExpr <- tab$logCPM
  tab$gene_symbol <- rownames(tab)
  tab$dataset <- dataset; tab$contrast_id <- contrast; tab$modality <- modality; tab$boundary_note <- note
  tab <- tab[,c("gene_symbol","dataset","contrast_id","modality","logFC","AveExpr","P.Value","adj.P.Val","boundary_note")]
  path <- file.path(out_dir, paste0(dataset, "_", contrast, ".tsv.gz"))
  write.table(tab, gzfile(path), sep="\t", quote=FALSE, row.names=FALSE)
  data.frame(dataset=dataset, contrast_id=contrast, genes=nrow(tab),
             fdr05=sum(tab$adj.P.Val < .05, na.rm=TRUE),
             fdr05_abslogfc1=sum(tab$adj.P.Val < .05 & abs(tab$logFC)>=1, na.rm=TRUE),
             note=note, stringsAsFactors=FALSE)
}

limma_group <- function(dataset, mat, groups, contrasts, modality, note="") {
  group <- factor(groups); design <- model.matrix(~0+group); colnames(design) <- levels(group)
  fit <- eBayes(lmFit(mat, design), robust=TRUE)
  out <- list()
  for (nm in names(contrasts)) {
    fit2 <- eBayes(contrasts.fit(fit, makeContrasts(contrasts=contrasts[[nm]], levels=design)), robust=TRUE)
    tab <- topTable(fit2, number=Inf, sort.by="none")
    out[[length(out)+1]] <- write_result(tab, dataset, nm, modality, note)
  }
  do.call(rbind,out)
}

edger_group <- function(dataset, mat, groups, contrasts, modality, note="") {
  group <- factor(groups); design <- model.matrix(~0+group); colnames(design) <- levels(group)
  y <- DGEList(counts=mat); keep <- filterByExpr(y, design=design); y <- y[keep,,keep.lib.sizes=FALSE]
  y <- calcNormFactors(y); y <- estimateDisp(y, design, robust=TRUE); fit <- glmQLFit(y, design, robust=TRUE)
  out <- list()
  for (nm in names(contrasts)) {
    q <- glmQLFTest(fit, contrast=makeContrasts(contrasts=contrasts[[nm]], levels=design))
    tab <- topTags(q, n=Inf, sort.by="none")$table
    names(tab)[names(tab)=="FDR"] <- "adj.P.Val"; names(tab)[names(tab)=="PValue"] <- "P.Value"
    out[[length(out)+1]] <- write_result(tab, dataset, nm, modality, note)
  }
  do.call(rbind,out)
}

summaries <- list()

# Human macrophage Affymetrix expression
ds <- "GSE61535"; m <- collapse_gene(ds, read_mat(file.path(matrix_dir,"GSE61535_log2_expression.tsv.gz")))
d <- design_all[design_all$dataset==ds,]; d <- d[match(colnames(m),d$gsm),]
summaries[[length(summaries)+1]] <- limma_group(ds,m,d$infection_group,
  list(A1_WT_vs_UNINFECTED="WT-UNINFECTED",D1_WT_vs_ANKB="WT-ANKB_MUTANT",D2_WT_vs_T2SS="WT-T2SS_MUTANT"),"RNA_ARRAY")

# Human macrophage RNA-seq FPKM; never treated as counts
ds <- "GSE135803"; raw <- collapse_gene(ds, read_mat(file.path(matrix_dir,"GSE135803_fpkm.tsv.gz")))
keep <- rowSums(raw>=1,na.rm=TRUE)>=2; m <- log2(raw[keep,,drop=FALSE]+0.1)
d <- design_all[design_all$dataset==ds,]; d <- d[match(colnames(m),d$gsm),]
summaries[[length(summaries)+1]] <- limma_group(ds,m,d$infection_group,
  list(A2_WT_vs_UNINFECTED="WT-UNINFECTED",D3_WT_vs_ANKH="WT-ANKH_MUTANT"),"RNA_FPKM","FPKM log2(x+0.1) linear model")

# Human THP-1 integer counts
ds <- "GSE287755"; m <- collapse_gene(ds, read_mat(file.path(matrix_dir,"GSE287755_raw_counts.tsv.gz")), TRUE)
d <- design_all[design_all$dataset==ds,]; d <- d[match(colnames(m),d$gsm),]
summaries[[length(summaries)+1]] <- edger_group(ds,m,d$infection_group,
  list(A3_WT_vs_UNINFECTED="WT-UNINFECTED",B1_WT_vs_DOTA="WT-DOTA_MUTANT",C1_WT_vs_DELTA4="WT-DELTA4_TRANSLATION_EFFECTOR"),"RNA_COUNTS")

# Mouse macrophage dual-color ratios; technical replicates collapsed first
ds <- "GSE26490"; m <- collapse_gene(ds, read_mat(file.path(matrix_dir,"GSE26490_normalized_log_ratio.tsv.gz")))
d <- design_all[design_all$dataset==ds,]; d <- d[match(colnames(m),d$gsm),]
keys <- paste(d$infection_group,d$biological_replicate,sep="__")
m2 <- sapply(unique(keys),function(k) rowMeans(m[,keys==k,drop=FALSE],na.rm=TRUE)); if(is.null(dim(m2))) m2 <- matrix(m2,ncol=1)
colnames(m2) <- unique(keys); groups <- sub("__.*$","",colnames(m2))
summaries[[length(summaries)+1]] <- limma_group(ds,m2,groups,
  list(A4_WT_vs_UNINFECTED="WT-UNINFECTED",C2_WT_vs_DELTA5="WT-DELTA5_TRANSLATION_EFFECTOR"),"RNA_ARRAY_RATIO","Technical replicates averaged before inference")

# Mouse RNA counts only, WT genotype and untreated
ds <- "GSE89184"; allm <- collapse_gene(ds, read_mat(file.path(matrix_dir,"GSE89184_raw_counts.tsv.gz")), TRUE)
d <- design_all[design_all$dataset==ds,]; d <- d[match(colnames(allm),d$gsm),]
sel <- d$assay_modality=="RNA" & d$host_genotype=="WT" & d$drug_treatment=="NONE"
m <- allm[,sel,drop=FALSE]; dr <- d[sel,]
summaries[[length(summaries)+1]] <- edger_group(ds,m,dr$infection_group,
  list(A5_FLAA_vs_UNINFECTED="FLAA_MUTANT-UNINFECTED",B2_FLAA_vs_DOTAFLAA="FLAA_MUTANT-DOTA_FLAA_MUTANT",C3_FLAA_vs_SEVEN="FLAA_MUTANT-SEVEN_EFFECTOR_FLAA_MUTANT"),"RNA_COUNTS","Mouse flaA background; NAIP5/NLRC4 boundary")

# Relative RNA-by-RPF interaction; does not estimate absolute global translation
sel <- d$host_genotype=="WT" & d$drug_treatment=="NONE"
m <- allm[,sel,drop=FALSE]; dt <- d[sel,]
g <- factor(paste(dt$assay_modality,dt$infection_group,sep="_"))
design <- model.matrix(~0+g); colnames(design) <- levels(g)
y <- DGEList(counts=m); keep <- filterByExpr(y,design=design); y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y); y <- estimateDisp(y,design,robust=TRUE); fit <- glmQLFit(y,design,robust=TRUE)
te_contrasts <- list(
  E1_TE_FLAA="(RPF_FLAA_MUTANT-RPF_UNINFECTED)-(RNA_FLAA_MUTANT-RNA_UNINFECTED)",
  E2_TE_DOTA="(RPF_DOTA_FLAA_MUTANT-RPF_UNINFECTED)-(RNA_DOTA_FLAA_MUTANT-RNA_UNINFECTED)",
  E3_TE_SEVEN="(RPF_SEVEN_EFFECTOR_FLAA_MUTANT-RPF_UNINFECTED)-(RNA_SEVEN_EFFECTOR_FLAA_MUTANT-RNA_UNINFECTED)")
for(nm in names(te_contrasts)) {
  q <- glmQLFTest(fit,contrast=makeContrasts(contrasts=te_contrasts[[nm]],levels=design)); tab <- topTags(q,n=Inf,sort.by="none")$table
  names(tab)[names(tab)=="FDR"] <- "adj.P.Val"; names(tab)[names(tab)=="PValue"] <- "P.Value"
  summaries[[length(summaries)+1]] <- write_result(tab,ds,nm,"RNA_RPF_INTERACTION","Relative TE interaction only; no absolute global-shutdown claim")
}

# Whole-lung dye-swap: orient both arrays as infected minus PBS
ds <- "GSE60085"; m <- collapse_gene(ds, read_mat(file.path(matrix_dir,"GSE60085_normalized_log_ratio.tsv.gz")))
oriented <- cbind(-m[,"GSM1464787"],m[,"GSM1464788"]); colnames(oriented) <- c("dyeswap1","dyeswap2")
fit <- eBayes(lmFit(oriented,matrix(1,nrow=2,ncol=1)),robust=TRUE); tab <- topTable(fit,coef=1,number=Inf,sort.by="none")
summaries[[length(summaries)+1]] <- write_result(tab,ds,"A6_WT_LP_vs_WT_PBS","RNA_ARRAY_RATIO","Whole lung; dye-swap signs explicitly oriented")

summary <- do.call(rbind,summaries)
write.table(summary,file.path(out_dir,"stage4_contrast_summary.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
stats <- list(
  contrasts_completed=nrow(summary), datasets_completed=length(unique(summary$dataset)),
  total_fdr05=sum(summary$fdr05), total_fdr05_abslogfc1=sum(summary$fdr05_abslogfc1),
  global_translation_claim_authorized=FALSE, holdout_accessed=FALSE,
  gate4_differential_expression=if(nrow(summary)==17 && length(unique(summary$dataset))==6) "PASS" else "FAIL"
)
json <- toJSON(stats,auto_unbox=TRUE,pretty=TRUE)
writeLines(json,file.path(out_dir,"stage4_statistics.json"),useBytes=TRUE)
writeLines(json,file.path(root,"checkpoints","stage4_differential_expression.done"),useBytes=TRUE)
cat(json,"\n")
if(stats$gate4_differential_expression!="PASS") quit(status=2)
