# Figure 2. Cross-cohort, mechanism-resolved and spatial evidence
ROOT <- Sys.getenv("LEGIONELLA_ROOT")
if (!nzchar(ROOT)) stop("Set LEGIONELLA_ROOT before running this script.")
if (!dir.exists(ROOT)) stop("LEGIONELLA_ROOT does not exist: ", ROOT)
setwd(file.path(ROOT, "15_manuscript/r_figures"))
source(file.path(ROOT, "15_manuscript/r_figures/_theme.R"))
SRC <- file.path(ROOT, "11_figures/source_data")
suppressMessages({library(dplyr); library(tidyr); library(readr); library(patchwork)})

A <- read_tsv(file.path(SRC, "Figure2A_discovery_logFC_matrix.tsv"))
B <- read_tsv(file.path(SRC, "Figure2B_primary_discovery_effects.tsv"))
C <- read_tsv(file.path(SRC, "Figure2C_mechanism_logFC_matrix.tsv"))
D <- read_tsv(file.path(SRC, "Figure2D_lcv_spatial_matrix.tsv"))

# ---- curate rows: 4 primary + key tier2 (keep panels readable) ----
key <- c("RAB1A", "RAB10", "ACTR3", "VAPA", "GNA13", "HNRNPA2B1", "HSPA8", "PRRC2C",
         "ACBD3", "TPM3", "RAB7A", "RAB8A")
A <- A %>% filter(human_symbol %in% key)
C <- C %>% filter(human_symbol %in% key)
D <- D %>% filter(human_symbol %in% key)

# ---- Panel A: discovery log2FC heatmap ----
A_long <- A %>% pivot_longer(-human_symbol, names_to = "dataset", values_to = "logFC")
pA <- ggplot(A_long, aes(dataset, human_symbol, fill = logFC)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C1121F", na.value = "#E8E8E8",
                       midpoint = 0, name = "log2FC") +
  labs(x = NULL, y = NULL, title = "A. Discovery-cohort log2FC") +
  theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
                      axis.text.y = element_text(size = 10))

# ---- Panel B: primary discovery effects (per-study, significance) ----
pB <- ggplot(B, aes(dataset, gene_symbol)) +
  geom_tile(colour = "white", linewidth = 0.5, aes(fill = logFC)) +
  geom_point(aes(alpha = study_fdr_hit), size = 3, shape = 15) +
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C1121F", midpoint = 0, name = "log2FC") +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.15), name = "study FDR<0.05") +
  labs(x = NULL, y = NULL, title = "B. Primary candidates per study") +
  theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
                      axis.text.y = element_text(size = 10))

# ---- Panel C: mechanism route logFC heatmap ----
C_long <- C %>% pivot_longer(-human_symbol, names_to = "route", values_to = "logFC")
pC <- ggplot(C_long, aes(route, human_symbol, fill = logFC)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C1121F", na.value = "#E8E8E8",
                       midpoint = 0, name = "log2FC") +
  labs(x = NULL, y = NULL, title = "C. Mechanism routes") +
  theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
                      axis.text.y = element_text(size = 10))

# ---- Panel D: LCV-associated detection ----
D_long <- D %>% pivot_longer(-human_symbol, names_to = "source", values_to = "detected")
pD <- ggplot(D_long, aes(source, human_symbol, fill = factor(detected))) +
  geom_tile(colour = "white", linewidth = 0.5) +
  scale_fill_manual(values = c("0" = "#F5F5F5", "1" = "red"), name = "detected") +
  labs(x = NULL, y = NULL, title = "D. LCV-associated proteome") +
  theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
                      axis.text.y = element_text(size = 10))

# taller canvas: 2 rows x 2 cols, each panel gets real area
fig <- (pA + pB + pC + pD + plot_layout(ncol = 2)) & panel_labels_topleft()
save_pub_pdf(fig, "Figure2_cross_cohort_mechanism", w = 11, h = 8.5)
