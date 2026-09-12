suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
})

args <- commandArgs(trailingOnly=FALSE)
script_arg <- args[grep("^--file=", args)]
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."), winslash="/", mustWork=TRUE)
matrix_dir <- file.path(root, "05_expression", "matrices")
raw_meta <- file.path(root, "01_data_raw", "metadata")
out_dir <- file.path(root, "05_expression", "gene_maps")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(root, "checkpoints"), recursive=TRUE, showWarnings=FALSE)

read_features <- function(path) {
  x <- read.delim(gzfile(path), check.names=FALSE, stringsAsFactors=FALSE)
  as.character(x[[1]])
}

read_platform <- function(path) {
  x <- readLines(gzfile(path), warn=FALSE)
  begin <- grep("^!platform_table_begin", x)
  end <- grep("^!platform_table_end", x)
  if (!length(begin) || !length(end)) stop("Platform table missing: ", path)
  text <- paste(x[(begin[1]+1):(end[1]-1)], collapse="\n")
  read.delim(text=text, quote="", comment.char="", check.names=FALSE, stringsAsFactors=FALSE, fill=TRUE)
}

clean_symbol <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "---", "NA", "null")] <- NA_character_
  x
}

write_map <- function(dataset, feature, symbol, source, ambiguous=FALSE) {
  symbol <- clean_symbol(symbol)
  tab <- data.frame(feature_id=as.character(feature), gene_symbol=symbol,
                    mapping_source=source, ambiguous_mapping=as.logical(ambiguous),
                    stringsAsFactors=FALSE)
  tab <- tab[!duplicated(tab$feature_id),]
  path <- file.path(out_dir, paste0(dataset, "_gene_map.tsv"))
  write.table(tab, path, sep="\t", quote=FALSE, row.names=FALSE, na="")
  data.frame(dataset=dataset, features=nrow(tab), mapped=sum(!is.na(tab$gene_symbol)),
             mapping_fraction=mean(!is.na(tab$gene_symbol)),
             ambiguous=sum(tab$ambiguous_mapping, na.rm=TRUE), stringsAsFactors=FALSE)
}

qc <- list()

f <- read_features(file.path(matrix_dir, "GSE135803_fpkm.tsv.gz"))
s <- mapIds(org.Hs.eg.db, keys=unique(f), keytype="ENSEMBL", column="SYMBOL", multiVals="first")
qc[[length(qc)+1]] <- write_map("GSE135803", f, unname(s[f]), "org.Hs.eg.db:ENSEMBL")

f <- read_features(file.path(matrix_dir, "GSE287755_raw_counts.tsv.gz"))
s <- mapIds(org.Hs.eg.db, keys=unique(f), keytype="ENSEMBL", column="SYMBOL", multiVals="first")
qc[[length(qc)+1]] <- write_map("GSE287755", f, unname(s[f]), "org.Hs.eg.db:ENSEMBL")

f <- read_features(file.path(matrix_dir, "GSE89184_raw_counts.tsv.gz"))
s <- mapIds(org.Mm.eg.db, keys=unique(f), keytype="REFSEQ", column="SYMBOL", multiVals="first")
qc[[length(qc)+1]] <- write_map("GSE89184", f, unname(s[f]), "org.Mm.eg.db:REFSEQ")

p <- read_platform(file.path(raw_meta, "GSE61535_family.soft.gz"))
entries <- strsplit(as.character(p$gene_assignment), " /// ", fixed=TRUE)
symbols <- lapply(entries, function(z) {
  y <- unique(vapply(strsplit(z, " // ", fixed=TRUE), function(w) if(length(w)>=2) trimws(w[2]) else NA_character_, character(1)))
  y <- y[!is.na(y) & y != "" & y != "---"]
  y
})
first_symbol <- vapply(symbols, function(z) if(length(z)) z[1] else NA_character_, character(1))
ambiguous <- vapply(symbols, length, integer(1)) > 1
qc[[length(qc)+1]] <- write_map("GSE61535", p$ID, first_symbol, "GPL6244:gene_assignment", ambiguous)

p <- read_platform(file.path(raw_meta, "GSE26490_family.soft.gz"))
qc[[length(qc)+1]] <- write_map("GSE26490", p$ID, p$Symbol, "GPL7275:Symbol")

p <- read_platform(file.path(raw_meta, "GSE60085_family.soft.gz"))
qc[[length(qc)+1]] <- write_map("GSE60085", p$ID, p$GENE_SYMBOL, "GPL4134:GENE_SYMBOL")

qc <- do.call(rbind, qc)
write.table(qc, file.path(out_dir, "stage3_gene_mapping_qc.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

stats <- list(
  datasets_mapped=nrow(qc),
  minimum_mapping_fraction=min(qc$mapping_fraction),
  mapping_pass=all(qc$mapping_fraction >= 0.20),
  holdout_accessed=FALSE,
  gate3_gene_mapping=if(all(qc$mapping_fraction >= 0.20)) "PASS" else "FAIL"
)
json <- jsonlite::toJSON(stats, auto_unbox=TRUE, pretty=TRUE)
writeLines(json, file.path(out_dir, "stage3_gene_mapping_statistics.json"), useBytes=TRUE)
writeLines(json, file.path(root, "checkpoints", "stage3_gene_mapping.done"), useBytes=TRUE)
cat(json, "\n")
if (!stats$mapping_pass) quit(status=2)
