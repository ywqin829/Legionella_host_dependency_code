# Figure 1. Data and evidence architecture of the Legionella host-dependency atlas
# Source data: 11_figures/source_data/Figure1*.tsv
ROOT <- Sys.getenv("LEGIONELLA_ROOT")
if (!nzchar(ROOT)) stop("Set LEGIONELLA_ROOT before running this script.")
if (!dir.exists(ROOT)) stop("LEGIONELLA_ROOT does not exist: ", ROOT)
setwd(file.path(ROOT, "15_manuscript/r_figures"))
source(file.path(ROOT, "15_manuscript/r_figures/_theme.R"))
SRC <- file.path(ROOT, "11_figures/source_data")
suppressMessages({library(dplyr); library(readr); library(patchwork)})

A <- read_tsv(file.path(SRC, "Figure1A_discovery_datasets.tsv")) %>%
  filter(analysis_route != "HOLDOUT_LOCKED")
if (anyNA(A$organism) || !all(A$organism %in% c("Homo sapiens", "Mus musculus"))) {
  stop("Figure 1A organism values must be complete and limited to Homo sapiens or Mus musculus")
}
A$organism <- factor(A$organism, levels = c("Homo sapiens", "Mus musculus"))
B <- read_tsv(file.path(SRC, "Figure1B_ppi_classes.tsv"))
C <- read_tsv(file.path(SRC, "Figure1C_crispr_tiers.tsv"))
D <- read_tsv(file.path(SRC, "Figure1D_evidence_funnel.tsv"))

# Panel A: discovery datasets - samples by organism
pA <- ggplot(A, aes(x = reorder(accession, n_samples), y = n_samples, fill = organism)) +
  geom_col(colour = "black", linewidth = 0.3, width = 0.7) +
  scale_fill_manual(values = c("Homo sapiens" = "#1B4965", "Mus musculus" = "#C9A227")) +
  coord_flip() +
  labs(x = NULL, y = "Samples", fill = NULL,
       title = "A. Discovery datasets") +
  theme_pub() + theme(legend.position = "top")

# Panel B: PPI evidence classes
pB <- ggplot(B, aes(x = evidence_class, y = pairs, fill = evidence_class)) +
  geom_col(colour = "black", linewidth = 0.3, width = 0.6) +
  geom_text(aes(label = pairs), vjust = -0.25, size = 3.5) +
  scale_fill_manual(values = c(A = "#C2185B", B = "#8E24AA", C = "#B0BEC5")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(x = "Evidence class", y = "Effector-host pairs", fill = NULL,
       title = "B. Physical interaction") +
  theme_pub() + theme(legend.position = "none")

# Panel C: CRISPR tiers
pC <- ggplot(C, aes(x = tier, y = genes, fill = tier)) +
  geom_col(colour = "black", linewidth = 0.3, width = 0.6) +
  geom_text(aes(label = genes), vjust = -0.25, size = 3.5) +
  scale_fill_manual(values = c("#1565C0", "#42A5F5", "#90CAF9")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(x = NULL, y = "Genes", fill = NULL,
       title = "C. Replicated CRISPR tiers") +
  theme_pub() + theme(legend.position = "none", axis.text.x = element_text(size = 8))

# Panel D: evidence funnel (short labels, no newlines)
D <- D %>% mutate(stage_short = c("DE calls", "Meta-convergent\ngenes", "Mechanism", "Candidates", "Tier 1")) %>%
  mutate(stage = factor(stage_short, levels = stage_short))
pD <- ggplot(D, aes(x = stage, y = n)) +
  geom_col(fill = "#2E8B57", colour = "black", linewidth = 0.3, width = 0.7) +
  geom_text(aes(label = format(n, big.mark = ",")), vjust = -0.25, size = 3.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(x = NULL, y = "Genes / candidates", title = "D. Selective evidence convergence") +
  theme_pub() + theme(axis.text.x = element_text(size = 9, angle = 25, hjust = 1))

fig <- (pA + pB + pC + pD + plot_layout(ncol = 2)) & panel_labels_topleft()
save_pub_pdf(fig, "Figure1_data_architecture", w = 10, h = 8.5)
