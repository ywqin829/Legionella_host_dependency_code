# =====================================================================
# Revision 2 / WP8: Supplementary Figure 7 (discovery-only candidate-rule
# sensitivity) rendered to match the manuscript legend A-D.
#
# The revision-2 WP4 renderer produced a two-panel summary. The legend in the
# working draft declares four panels: (A) panel size, overlap and Jaccard for
# the five distinct local variants; (B) per-gene retention frequency across
# those five variants; (C) panel, RAB and curated trafficking counts as the
# required number of discovery-evidence domains rises from three to five;
# (D) Tier 1 counts for tier-only perturbations that leave the 56-gene panel
# unchanged.
#
# Every value is read from the archived WP4 outputs. No statistic is
# recomputed, no panel is re-selected, and the holdout is not touched.
#
# Run from this directory:  Rscript suppfigure7_rule_sensitivity.R
# =====================================================================

.get_script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1]))))
  getwd()
}

SCRIPT_DIR <- .get_script_dir()
ROOT <- Sys.getenv("LEGIONELLA_ROOT", "C:/Users/yw_qi/Legionella_host_dependency")
FROZEN_FIG <- file.path(ROOT, "FINAL_SUBMISSION_2026-10-04", "figure R script")
REV2 <- file.path(ROOT, "revision2_2026-11-23")
OUT <- file.path(REV2, "figures")
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)

source(file.path(FROZEN_FIG, "_theme.R"))
suppressMessages({library(dplyr); library(readr); library(tidyr); library(patchwork)})

save_png <- function(p, file, w, h, dpi = 300) {
  ggsave(file.path(OUT, paste0(file, ".png")), p, width = w, height = h,
         dpi = dpi, bg = "white")
  message("wrote ", file, ".png")
}
save_pdf <- function(p, file, w, h) {
  ggsave(file.path(OUT, paste0(file, ".pdf")), p, width = w, height = h,
         device = cairo_pdf)
  message("wrote ", file, ".pdf")
}

tab <- read_tsv(file.path(REV2, "REVISION2_WP4_RULE_SENSITIVITY_TABLE.tsv"),
                show_col_types = FALSE)
ret <- read_tsv(file.path(REV2, "REVISION2_VARIANT_GENE_RETENTION.tsv"),
                show_col_types = FALSE)

# ------------------------------------------------------------------ #
# Panel A: panel size and Jaccard for the five distinct local variants
# ------------------------------------------------------------------ #
local_ids <- c("A1_admit_F3", "A2_fitness_does_not_block",
               "A4_T2_orthogonal_no_spatial",
               "A6_T2_orthogonal_require_response_and_mechanism",
               "A7_drop_T2_functional_response_branch")
local_lab <- c(A1_admit_F3 = "A1  admit F3 evidence",
               A2_fitness_does_not_block = "A2  fitness does not block",
               A4_T2_orthogonal_no_spatial = "A4  no spatial domain",
               A6_T2_orthogonal_require_response_and_mechanism = "A6  require response + mechanism",
               A7_drop_T2_functional_response_branch = "A7 (= B1)  drop functional/response branch")

loc <- tab %>%
  filter(variant %in% local_ids) %>%
  mutate(label = factor(local_lab[variant], levels = rev(unname(local_lab))),
         jaccard = jaccard_vs_original,
         removed = removed_n,
         added = added_n)

ref <- tab %>% filter(variant == "ORIGINAL")

pA <- ggplot(loc, aes(panel_n, label)) +
  geom_col(fill = "#D9A066", colour = "black", linewidth = 0.3, width = 0.66) +
  geom_vline(xintercept = ref$panel_n, linetype = 2, colour = "#8C2D2D") +
  geom_text(aes(x = 104, label = sprintf("%d genes; J = %.2f", panel_n, jaccard)),
            hjust = 1, size = 3.0) +
  scale_x_continuous(limits = c(0, 105), expand = expansion(mult = c(0, 0))) +
  labs(x = "Genes in the panel", y = NULL,
       title = "A. Local rule variants: size and overlap",
       subtitle = sprintf("Dashed line = frozen panel (%d genes); J = Jaccard vs frozen",
                          ref$panel_n)) +
  theme_pub(base_size = 10) +
  theme(axis.text.y = element_text(size = 8.5),
        plot.subtitle = element_text(size = 8))

# ------------------------------------------------------------------ #
# Panel B: per-gene retention across the five distinct local panels
# ------------------------------------------------------------------ #
freq <- ret %>%
  filter(n_local_panel_effect == 5) %>%
  count(n_retained_local, name = "genes") %>%
  complete(n_retained_local = 3:5, fill = list(genes = 0)) %>%
  mutate(n_retained_local = factor(n_retained_local, levels = 3:5))

pB <- ggplot(freq, aes(n_retained_local, genes)) +
  geom_col(fill = "#C45A3C", colour = "black", linewidth = 0.3, width = 0.62) +
  geom_text(aes(label = sprintf("%d / %d", genes, sum(freq$genes))),
            vjust = -0.35, size = 3.2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Local panels retaining the gene (of 5)", y = "Frozen panel genes (n)",
       title = "B. Per-gene retention across local variants") +
  theme_pub(base_size = 10)

# ------------------------------------------------------------------ #
# Panel C: structural discovery-domain series (3 -> 5 domains)
# ------------------------------------------------------------------ #
dom <- tab %>%
  filter(variant %in% c("D3_panel_min_discovery_domains",
                        "D4_panel_min_discovery_domains",
                        "D5_panel_min_discovery_domains")) %>%
  mutate(domains = factor(c(">= 3", ">= 4", ">= 5"), levels = c(">= 3", ">= 4", ">= 5"))) %>%
  select(domains, panel_n, rab_gtpase_retained, trafficking_retained) %>%
  pivot_longer(c(panel_n, rab_gtpase_retained, trafficking_retained),
               names_to = "series", values_to = "count") %>%
  mutate(series = factor(recode(series,
                                panel_n = "Panel genes",
                                rab_gtpase_retained = "RAB GTPases",
                                trafficking_retained = "Trafficking genes"),
                         levels = c("Panel genes", "RAB GTPases", "Trafficking genes")))

pC <- ggplot(dom, aes(domains, count, fill = series)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.7,
           colour = "black", linewidth = 0.3) +
  geom_text(aes(label = count), position = position_dodge(width = 0.78),
            vjust = -0.35, size = 3.0) +
  scale_fill_manual(values = c("Panel genes" = "#B0B0B0",
                               "RAB GTPases" = "#C45A3C",
                               "Trafficking genes" = "#8C2D2D"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(x = "Required discovery-evidence domains", y = "Genes retained",
       title = "C. Structural domain-threshold series") +
  theme_pub(base_size = 10) +
  theme(legend.position = "top",
        plot.caption = element_text(hjust = 0, size = 7, colour = "black")) +
  labs(caption = "The >=5-domain boundary admits one gene and is shown as a structural endpoint, not a stable result.")

# ------------------------------------------------------------------ #
# Panel D: Tier 1 counts under tier-only perturbations
# ------------------------------------------------------------------ #
tier_lab <- c(A3_T1_orthogonal_response_only = "A3  Tier 1: response only",
              A5_T1_require_response_and_mechanism = "A5  Tier 1: response + mechanism",
              A8_T1_require_four_domains = "A8  Tier 1: >= 4 domains",
              B2_drop_T3_single_route = "B2  drop T3 single-route")

tier <- tab %>%
  filter(measured_class == "tier_only") %>%
  mutate(label = factor(tier_lab[variant], levels = rev(unname(tier_lab))))

pD <- ggplot(tier, aes(tier1_n, label)) +
  geom_col(fill = "#7B6D8D", colour = "black", linewidth = 0.3, width = 0.62) +
  geom_vline(xintercept = ref$tier1_n, linetype = 2, colour = "#8C2D2D") +
  geom_text(aes(label = sprintf("Tier 1 = %d; Tier 2 = %d; panel = %d",
                                tier1_n, tier2_n, panel_n)),
            hjust = -0.08, size = 2.8) +
  scale_x_continuous(limits = c(0, 6.6), expand = expansion(mult = c(0, 0))) +
  labs(x = "Tier 1 genes", y = NULL,
       title = "D. Tier-only perturbations (panel unchanged)",
       subtitle = sprintf("Dashed line = frozen Tier 1 (%d genes); all panels keep %d genes",
                          ref$tier1_n, ref$panel_n)) +
  theme_pub(base_size = 10) +
  theme(axis.text.y = element_text(size = 8.5),
        plot.subtitle = element_text(size = 8))

figS7 <- (pA | pB) / (pC | pD) + plot_layout(heights = c(1, 1))
figS7 <- figS7 & panel_labels_topleft()
save_pdf(figS7, "SupplementaryFigure7_rule_sensitivity", w = 15.5, h = 10.4)
save_png(figS7, "SupplementaryFigure7_rule_sensitivity", w = 15.5, h = 10.4)

write_tsv(loc %>% select(variant, label, panel_n, jaccard, added, removed),
          file.path(OUT, "SuppFigure7A_source_data.tsv"))
write_tsv(freq, file.path(OUT, "SuppFigure7B_source_data.tsv"))
write_tsv(dom, file.path(OUT, "SuppFigure7C_source_data.tsv"))
write_tsv(tier %>% mutate(label = tier_lab[variant]) %>%
            select(variant, label, tier1_n, tier2_n, panel_n),
          file.path(OUT, "SuppFigure7D_source_data.tsv"))

message("done")
