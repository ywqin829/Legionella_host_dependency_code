suppressPackageStartupMessages({library(babelgene); library(jsonlite)})
args <- commandArgs(trailingOnly=FALSE); script_arg <- args[grep("^--file=",args)]
root <- normalizePath(file.path(dirname(sub("^--file=","",script_arg)),".."),winslash="/",mustWork=TRUE)
de_dir <- file.path(root,"06_differential_expression"); out_dir <- file.path(root,"07_convergence")
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)

files <- list.files(de_dir,pattern="_(B[12]|C[123]|D[123]|E[123])_.*\\.tsv\\.gz$",full.names=TRUE)
x <- do.call(rbind,lapply(files,function(f) read.delim(gzfile(f),stringsAsFactors=FALSE,check.names=FALSE)))
x$route <- ifelse(grepl("^B",x$contrast_id),"DOT_ICM",
           ifelse(grepl("^C",x$contrast_id),"TRANSLATION_EFFECTOR_MODULE",
           ifelse(grepl("^D",x$contrast_id),"SINGLE_EFFECTOR_OR_T2SS","RNA_RPF_INTERACTION")))
x$species <- ifelse(x$dataset %in% c("GSE26490","GSE89184"),"MOUSE","HUMAN")
ortho <- orthologs(unique(x$gene_symbol[x$species=="MOUSE"]),species="mouse",human=FALSE,min_support=3,top=TRUE)
omap <- setNames(ortho$human_symbol,ortho$symbol)
x$human_symbol <- ifelse(x$species=="HUMAN",toupper(x$gene_symbol),unname(omap[x$gene_symbol]))
x <- x[!is.na(x$human_symbol)&nzchar(x$human_symbol)&is.finite(x$logFC)&is.finite(x$P.Value),]
x$signed_z <- sign(x$logFC)*qnorm(pmax(pmin(x$P.Value,1-1e-15),1e-300)/2,lower.tail=FALSE)
write.table(x,file.path(out_dir,"mechanism_gene_study_evidence.tsv.gz"),sep="\t",quote=FALSE,row.names=FALSE)

summ <- function(d){
 dir <- if(sum(d$logFC>0)>=sum(d$logFC<0)) "UP" else "DOWN"; conc <- if(dir=="UP") sum(d$logFC>0) else sum(d$logFC<0)
 z <- sum(d$signed_z)/sqrt(nrow(d))
 data.frame(human_symbol=d$human_symbol[1],route=d$route[1],n_studies=length(unique(d$dataset)),n_contrasts=nrow(d),
   consensus_direction=dir,direction_concordance=conc/nrow(d),combined_z=z,combined_p=2*pnorm(-abs(z)),
   min_study_fdr=min(d$adj.P.Val,na.rm=TRUE),n_study_fdr05=sum(d$adj.P.Val<.05,na.rm=TRUE),
   datasets=paste(sort(unique(d$dataset)),collapse=";"),stringsAsFactors=FALSE)
}
key <- paste(x$route,x$human_symbol,sep="::"); s <- do.call(rbind,lapply(split(x,key),summ)); rownames(s)<-NULL
s$combined_fdr <- ave(s$combined_p,s$route,FUN=function(p)p.adjust(p,"BH"))
s$mechanism_convergent <- with(s,
  (route=="DOT_ICM" & n_studies>=2 & direction_concordance>=.75 & combined_fdr<.05 & n_study_fdr05>=1) |
  (route=="TRANSLATION_EFFECTOR_MODULE" & n_studies>=2 & direction_concordance>=2/3 & combined_fdr<.05 & n_study_fdr05>=1) |
  (route=="SINGLE_EFFECTOR_OR_T2SS" & n_studies>=1 & combined_fdr<.05 & n_study_fdr05>=1) |
  (route=="RNA_RPF_INTERACTION" & n_studies>=1 & combined_fdr<.05 & n_study_fdr05>=1))
s$negative_interpretation <- ifelse(s$route=="SINGLE_EFFECTOR_OR_T2SS" & !s$mechanism_convergent,
  "NO_SIGNAL_IS_NOT_PROOF_OF_NO_EFFECT; REDUNDANCY_OR_POWER_MAY_MASK","")
s <- s[order(s$route,s$combined_fdr),]
write.table(s,file.path(out_dir,"mechanism_convergence.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
stats <- list(contrasts=length(unique(x$contrast_id)),genes=nrow(s),convergent=sum(s$mechanism_convergent),
  by_route=as.list(table(s$route)),holdout_accessed=FALSE,
  gate5_mechanism=if(length(unique(x$contrast_id))==11) "PASS" else "FAIL")
json <- toJSON(stats,auto_unbox=TRUE,pretty=TRUE); writeLines(json,file.path(out_dir,"stage5_mechanism_statistics.json"),useBytes=TRUE)
writeLines(json,file.path(root,"checkpoints","stage5_mechanism_convergence.done"),useBytes=TRUE); cat(json,"\n")
if(stats$gate5_mechanism!="PASS") quit(status=2)
