# Figure 4. Locked cross-context transcriptional corroboration in E-MTAB-16655 (PMA-THP-1)
ROOT <- Sys.getenv("LEGIONELLA_ROOT")
if (!nzchar(ROOT)) ROOT <- Sys.getenv("LEGIONELLA_ROOT")
if (!nzchar(ROOT)) stop("Set LEGIONELLA_ROOT before running this script.")
if (!dir.exists(ROOT)) stop("LEGIONELLA_ROOT does not exist: ", ROOT)
setwd(file.path(ROOT, "15_manuscript/r_figures"))
source(file.path(ROOT, "15_manuscript/r_figures/_theme.R"))
SRC <- file.path(ROOT, "11_figures/source_data")
suppressMessages({library(dplyr); library(readr); library(patchwork)})

A <- read_tsv(file.path(SRC, "Figure4A_holdout_candidate_volcano.tsv"))
B <- read_tsv(file.path(SRC, "Figure4B_primary_holdout_effects.tsv"))
C <- read_tsv(file.path(SRC, "Figure4C_holdout_set_tests.tsv"))
D <- read_tsv(file.path(SRC, "Figure4D_effector_boundary_tests.tsv"))

# ---- Panel A: holdout candidate volcano ----
A <- A %>% mutate(
  tier_group = case_when(final_tier == "TIER1_EXTERNALLY_VALIDATED" ~ "Tier 1 (RAB1A)",
                         final_tier == "REJECT_HOLDOUT_DIRECTIONAL_CONTRADICTION" ~ "Rejected",
                         grepl("TIER2", final_tier) ~ "Tier 2",
                         TRUE ~ "Tier 3/other"),
  highlight = human_symbol %in% c("RAB1A", "RAB10", "TPM3", "ACTR3", "VAPA"))
pA <- ggplot(A, aes(logFC, minuslog10FDR, colour = tier_group)) +
  geom_point(aes(shape = highlight), size = 1.6, alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  ggrepel::geom_text_repel(data = filter(A, highlight),
                           aes(label = human_symbol), size = 3.2, max.overlaps = 20, show.legend = FALSE) +
  scale_colour_manual(values = c("Tier 1 (RAB1A)" = "#8C2D2D", "Tier 2" = "#C65D3A",
                                 "Tier 3/other" = "#D89B5B", "Rejected" = "#6B3E2E")) +
  scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 16), guide = "none") +
  labs(x = "log2FC (WT vs uninfected)", y = "-log10 FDR", colour = NULL,
       title = "A. Frozen panel in holdout") +
  theme_pub()

# ---- Panel B: primary candidates with 95% CI ----
pB <- ggplot(B, aes(x = reorder(human_symbol, logFC), y = logFC, colour = final_tier)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15, linewidth = 0.6) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c("TIER1_EXTERNALLY_VALIDATED" = "#8C2D2D",
                                 "TIER2_MULTIMODAL" = "#C65D3A")) +
  labs(x = NULL, y = "log2FC (WT vs uninfected)", colour = NULL,
       title = "B. Primary candidates (95% CI)") +
  theme_pub() + theme(legend.position = "top")

# ---- Panel C: set-test FDR (compact display labels) ----
C <- C %>% mutate(short_label = case_when(
  set_name == "PRIMARY4_DIRECTIONAL" ~ "WTvsUI: primary4",
  set_name == "PANEL_DIRECTIONAL" ~ "WTvsUI: panel54",
  set_name == "DOT_ICM_PREDEFINED" ~ "WTvsdotA: Dot/Icm",
  grepl("^[A-Z]+_vs_WT$", contrast_id) ~ paste0(gsub("_vs_WT", "", contrast_id), " vs WT"),
  TRUE ~ set_name))
pC <- ggplot(C, aes(x = reorder(short_label, minuslog10FDR), y = minuslog10FDR,
                    fill = passes_frozen_endpoint)) +
  geom_col(colour = "black", linewidth = 0.3, width = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#C45A3C", "FALSE" = "#F0C9A8"), name = "passed 5% FDR") +
  labs(x = NULL, y = "-log10 empirical FDR", title = "C. Pre-specified set tests (15)") +
  theme_pub()

# ---- Panel D: effector-boundary tests ----
D <- D %>% mutate(contrast_id = gsub("_vs_WT", "", contrast_id))
pD <- ggplot(D, aes(x = reorder(contrast_id, score), y = score, fill = passes_frozen_endpoint)) +
  geom_col(colour = "black", linewidth = 0.3, width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#C45A3C", "FALSE" = "#F0C9A8"), name = "passed") +
  labs(x = NULL, y = "Mean absolute signed Z", title = "D. Multi-effector boundary tests") +
  theme_pub()

fig <- (pA | pB) / (pC | pD)
save_pub_pdf(fig, "Figure4_holdout_validation", w = 12, h = 9)
