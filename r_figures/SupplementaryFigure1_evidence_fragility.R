# Supplementary Figure 1. Leave-one-source-out evidence fragility
ROOT <- Sys.getenv("LEGIONELLA_ROOT")
if (!nzchar(ROOT)) ROOT <- Sys.getenv("LEGIONELLA_ROOT")
if (!nzchar(ROOT)) stop("Set LEGIONELLA_ROOT before running this script.")
if (!dir.exists(ROOT)) stop("LEGIONELLA_ROOT does not exist: ", ROOT)
setwd(file.path(ROOT, "15_manuscript/r_figures"))
source(file.path(ROOT, "15_manuscript/r_figures/_theme.R"))
SRC <- file.path(ROOT, "13_postfreeze_interpretation")
suppressMessages({library(dplyr); library(tidyr); library(readr); library(patchwork)})

D <- read_tsv(file.path(SRC, "a2_evidence_fragility_matrix.tsv"))
datasets <- c("GSE135803", "GSE26490", "GSE287755", "GSE60085", "GSE61535", "GSE89184")

# parse leave1_dataset strings: "GSE135803->UNSTABLE; GSE61535:UP/3; ..."
parse_l1 <- function(s) {
  toks <- trimws(strsplit(s, ";")[[1]])
  out <- setNames(rep("NOT_IN", length(datasets)), datasets)
  for (t in toks) {
    if (t == "" || is.na(t)) next
    ds <- sub(":.*|->.*", "", t)
    if (ds %in% datasets) out[ds] <- ifelse(grepl("->UNSTABLE", t), "FLIPPED", "STABLE")
  }
  out
}
mat <- t(sapply(D$leave1_dataset, parse_l1))
rownames(mat) <- D$candidate
# human-only column
mat <- cbind(mat, human_only = ifelse(D$human_only_stable, "STABLE", "UNSTABLE"))
mat <- mat[c("RAB1A", setdiff(rownames(mat), "RAB1A")), ]

M <- as.data.frame(mat) %>% mutate(gene = rownames(mat)) %>%
  pivot_longer(-gene, names_to = "source", values_to = "state") %>%
  mutate(state = factor(state, levels = c("STABLE", "FLIPPED", "UNSTABLE", "NOT_IN")))
M$source <- factor(M$source, levels = colnames(mat))

p <- ggplot(M, aes(source, gene, fill = state)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  scale_fill_manual(values = c("STABLE" = "#2E8B57", "FLIPPED" = "#F4D03F",
                               "UNSTABLE" = "#C1121F", "NOT_IN" = "#B0BEC5"),
                    name = "leave-one-out") +
  scale_x_discrete(labels = c("GSE135803" = "GSE135803", "GSE26490" = "GSE26490",
                              "GSE287755" = "GSE287755", "GSE60085" = "GSE60085",
                              "GSE61535" = "GSE61535", "GSE89184" = "GSE89184",
                              "human_only" = "Human-only")) +
  labs(x = "Leave-one-dataset-out source", y = NULL,
       title = "Direction preserved in all removals; only the stability flag is threshold-sensitive") +
  theme_pub() + theme(axis.text.y = element_text(size = 9))

save_pub_pdf(p, "SupplementaryFigure1_evidence_fragility", w = 9, h = 6.5)
