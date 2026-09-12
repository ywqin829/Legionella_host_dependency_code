suppressPackageStartupMessages({library(babelgene); library(jsonlite)})

args <- commandArgs(trailingOnly=FALSE)
script_arg <- args[grep("^--file=", args)]
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."), winslash="/", mustWork=TRUE)
de_dir <- file.path(root,"06_differential_expression")
out_dir <- file.path(root,"07_convergence")
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(root,"checkpoints"),recursive=TRUE,showWarnings=FALSE)

files <- list.files(de_dir,pattern="_A[1-6]_.*\\.tsv\\.gz$",full.names=TRUE)
if(length(files)!=6) stop("Expected six primary infection contrasts; found ",length(files))
rows <- lapply(files,function(f) read.delim(gzfile(f),stringsAsFactors=FALSE,check.names=FALSE))
x <- do.call(rbind,rows)
x$species <- ifelse(x$dataset %in% c("GSE26490","GSE89184","GSE60085"),"MOUSE","HUMAN")
x$original_symbol <- x$gene_symbol

mouse_symbols <- unique(x$gene_symbol[x$species=="MOUSE"])
ortho <- orthologs(mouse_symbols,species="mouse",human=FALSE,min_support=3,top=TRUE)
omap <- setNames(ortho$human_symbol,ortho$symbol)
x$human_symbol <- ifelse(x$species=="HUMAN",toupper(x$gene_symbol),unname(omap[x$gene_symbol]))
x <- x[!is.na(x$human_symbol) & nzchar(x$human_symbol) & is.finite(x$logFC) & is.finite(x$P.Value),]
x$signed_z <- sign(x$logFC)*qnorm(pmax(pmin(x$P.Value,1-1e-15),1e-300)/2,lower.tail=FALSE)
x$study_fdr_hit <- x$adj.P.Val < .05

# Preserve mapped row-level evidence, then collapse multiple mouse orthologues
# mapping to the same human symbol within one dataset. A dataset is one
# independent study unit and must never contribute more than once per gene.
write.table(x,file.path(out_dir,"primary_infection_gene_study_evidence_raw_mapped.tsv.gz"),sep="\t",quote=FALSE,row.names=FALSE)
collapse_dataset_unit <- function(d) {
  z <- mean(d$signed_z)
  data.frame(
    dataset=d$dataset[1], contrast_id=d$contrast_id[1], species=d$species[1],
    original_symbol=paste(sort(unique(d$original_symbol)),collapse=";"),
    gene_symbol=paste(sort(unique(d$gene_symbol)),collapse=";"),
    human_symbol=d$human_symbol[1], logFC=median(d$logFC),
    P.Value=2*pnorm(-abs(z)), adj.P.Val=min(d$adj.P.Val,na.rm=TRUE),
    signed_z=z, study_fdr_hit=any(d$study_fdr_hit),
    n_original_rows_collapsed=nrow(d), stringsAsFactors=FALSE)
}
x <- do.call(rbind,lapply(split(x,interaction(x$human_symbol,x$dataset,drop=TRUE)),collapse_dataset_unit))
rownames(x) <- NULL
write.table(x,file.path(out_dir,"primary_infection_gene_study_evidence.tsv.gz"),sep="\t",quote=FALSE,row.names=FALSE)

summarize_gene <- function(d) {
  direction <- if(sum(d$logFC>0)>=sum(d$logFC<0)) "UP" else "DOWN"
  concordant <- if(direction=="UP") sum(d$logFC>0) else sum(d$logFC<0)
  z <- sum(d$signed_z)/sqrt(nrow(d))
  data.frame(
    human_symbol=d$human_symbol[1], n_studies=nrow(d), n_human=sum(d$species=="HUMAN"), n_mouse=sum(d$species=="MOUSE"),
    consensus_direction=direction, direction_concordance=concordant/nrow(d), median_logFC=median(d$logFC),
    combined_z=z, combined_p=2*pnorm(-abs(z)), min_study_fdr=min(d$adj.P.Val,na.rm=TRUE),
    n_study_fdr05=sum(d$study_fdr_hit), datasets=paste(sort(unique(d$dataset)),collapse=";"), stringsAsFactors=FALSE)
}
split_all <- split(x,x$human_symbol)
summary <- do.call(rbind,lapply(split_all,summarize_gene)); rownames(summary)<-NULL
summary$combined_fdr <- p.adjust(summary$combined_p,"BH")
summary$naip_nlrc4_boundary <- summary$human_symbol %in% c("NAIP","NLRC4","CASP1","IL1B")
summary$stable_cross_cohort <- with(summary,n_studies>=3 & n_human>=2 & direction_concordance>=.75 & combined_fdr<.05 & n_study_fdr05>=1)
summary <- summary[order(summary$combined_fdr,-summary$direction_concordance),]
write.table(summary,file.path(out_dir,"primary_infection_convergence.tsv"),sep="\t",quote=FALSE,row.names=FALSE)

human_x <- x[x$species=="HUMAN",]
human_summary <- do.call(rbind,lapply(split(human_x,human_x$human_symbol),summarize_gene)); rownames(human_summary)<-NULL
human_summary$combined_fdr <- p.adjust(human_summary$combined_p,"BH")
human_summary$stable_human <- with(human_summary,n_studies>=2 & direction_concordance>=2/3 & combined_fdr<.05 & n_study_fdr05>=1)
human_summary <- human_summary[order(human_summary$combined_fdr),]
write.table(human_summary,file.path(out_dir,"primary_infection_human_only_convergence.tsv"),sep="\t",quote=FALSE,row.names=FALSE)

stats <- list(
  primary_contrasts=length(unique(x$contrast_id)), genes_tested=nrow(summary),
  stable_cross_cohort=sum(summary$stable_cross_cohort), stable_human=sum(human_summary$stable_human),
  orthologs_mapped=nrow(ortho),
  combination_method="equal-weight signed Stouffer; independence is an approximation",
  tier_assignment_authorized=FALSE, holdout_accessed=FALSE,
  gate5_expression_convergence=if(length(unique(x$dataset))==6 && nrow(summary)>5000) "PASS" else "FAIL"
)
json <- toJSON(stats,auto_unbox=TRUE,pretty=TRUE)
writeLines(json,file.path(out_dir,"stage5_expression_statistics.json"),useBytes=TRUE)
writeLines(json,file.path(root,"checkpoints","stage5_expression_convergence.done"),useBytes=TRUE)
cat(json,"\n")
if(stats$gate5_expression_convergence!="PASS") quit(status=2)
