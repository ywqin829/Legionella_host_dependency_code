# =====================================================================
# Revision 2 / WP8: split the nine-panel Figure 5 into
#   (1) a six-panel main figure that keeps the association estimates,
#       the ascertainment control, the mouse-layer stratification and the
#       abundance dose-response; and
#   (2) a three-panel supplementary figure holding the correlation-robust
#       diagnostics.
#
# The frozen 2026-10-04 scripts are NOT edited. This script re-reads the same
# archived figure source data and re-draws the same panels under new letters.
# Panel content is unchanged; only grouping, letters and layout differ.
#
# Run from this directory:  Rscript figure5_split.R
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
SRC <- file.path(FROZEN_FIG, "source_data")
OUT <- file.path(ROOT, "revision2_2026-11-23", "figures")
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)

if (!dir.exists(SRC)) stop("frozen source_data not found: ", SRC)
source(file.path(FROZEN_FIG, "_theme.R"))
suppressMessages({
  library(dplyr); library(readr); library(tidyr); library(patchwork); library(jsonlite)
})

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

coupling <- read_tsv(file.path(SRC, "a1_evidence_coupling_matched.tsv"), show_col_types = FALSE)
joint_rate <- read_tsv(file.path(SRC, "g_ppi_lcv_exact_rate_ci.tsv"), show_col_types = FALSE)
balance <- read_tsv(file.path(SRC, "g_matching_balance_by_covariate.tsv"), show_col_types = FALSE)
sensitivity <- read_tsv(file.path(SRC, "g_matching_sensitivity.tsv"), show_col_types = FALSE)

# ---------------------------------------------------------------- #
# MAIN FIGURE 5
# ---------------------------------------------------------------- #

# ---- Panel A: matched fold enrichment ----
coupling <- coupling %>% mutate(
  log2fold = log2(fold_vs_matched),
  star = case_when(
    perm_p < 0.001 ~ "***",
    perm_p < 0.01 ~ "**",
    perm_p < 0.05 ~ "*",
    TRUE ~ "ns"
  ),
  short = case_when(
    grepl("PPI\\(A/B\\) -> CRISPR hit", label) ~ "PPI(A/B) -> CRISPR hit",
    grepl("PPI\\(A/B\\) -> clean", label) ~ "PPI(A/B) -> filtered CRISPR evidence",
    grepl("PPI\\(any\\) -> clean", label) ~ "PPI(any) -> filtered CRISPR evidence",
    grepl("LCV -> CRISPR hit", label) ~ "LCV -> CRISPR hit",
    grepl("LCV -> clean", label) ~ "LCV -> filtered CRISPR evidence",
    grepl("PPI\\(A/B\\) -> LCV", label) ~ "PPI(A/B) -> LCV",
    grepl("CRISPR hit -> response", label) ~ "CRISPR hit -> response",
    grepl("clean CRISPR -> response", label) ~ "filtered CRISPR evidence -> response",
    grepl("response -> clean", label) ~ "response -> filtered CRISPR evidence",
    grepl("LCV -> response", label) ~ "LCV -> response",
    grepl("PPI\\(A/B\\) -> response", label) ~ "PPI(A/B) -> response",
    TRUE ~ label
  )) %>% mutate(short = factor(short, levels = rev(short)))

pA <- ggplot(coupling, aes(log2fold, short, fill = log2fold >= 1.25)) +
  geom_col(colour = "black", linewidth = 0.3, width = 0.7) +
  geom_text(aes(label = star), hjust = -0.3, size = 3.5) +
  geom_vline(xintercept = 0, colour = "black", linewidth = 0.5) +
  scale_x_continuous(expand = expansion(mult = c(0.12, 0.20))) +
  scale_fill_manual(values = c("TRUE" = "#C45A3C", "FALSE" = "#D9A066"), guide = "none") +
  labs(x = "log2 fold-enrichment vs matched background", y = NULL,
       title = "A. Evidence coupling") +
  theme_pub(base_size = 11) +
  theme(axis.text.y = element_text(size = 9), plot.margin = margin(5.5, 16, 5.5, 5.5))

# ---- Panel B: exact uncertainty for the joint evidence rate ----
joint_rate$group <- factor(joint_rate$group, levels = c("Universe", "PPI(A/B)", "LCV", "PPI(A/B)+LCV"))
pB <- ggplot(joint_rate, aes(group, rate, fill = group)) +
  geom_col(colour = "black", linewidth = 0.3, width = 0.7) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.16, linewidth = 0.6) +
  geom_text(aes(label = paste0(successes, "/", n)), vjust = -0.55, size = 3.4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  scale_fill_manual(values = c("Universe" = "#B07D62", "PPI(A/B)" = "#D9A066",
                               "LCV" = "#C45A3C", "PPI(A/B)+LCV" = "#8C2D2D"), guide = "none") +
  labs(x = NULL, y = "P(filtered CRISPR evidence | evidence)",
       title = "B. Joint evidence, exact 95% CI") +
  theme_pub(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 9))

# ---- Panel C: covariate balance ----
balance <- balance %>% mutate(
  covariate_label = recode(covariate,
    `_deg` = "Network degree",
    `_pmids` = "Publication count",
    `_chronos` = "U937 Chronos effect",
    `_common` = "Common-essentiality",
    `_ui` = "Baseline expression",
    `_nstud` = "Expression-study count",
    .default = covariate
  ),
  covariate_label = factor(covariate_label, levels = rev(unique(covariate_label)))
)
bal_long <- balance %>%
  select(covariate_label, abs_smd_before, abs_smd_after) %>%
  pivot_longer(c(abs_smd_before, abs_smd_after), names_to = "stage", values_to = "smd") %>%
  mutate(stage = factor(recode(stage,
                               abs_smd_before = "Before matching",
                               abs_smd_after = "After matching"),
                        levels = c("Before matching", "After matching")))

pC <- ggplot(bal_long, aes(y = covariate_label)) +
  geom_segment(data = balance,
               aes(x = abs_smd_before, xend = abs_smd_after, yend = covariate_label),
               colour = "#B07D62", linewidth = 0.8) +
  geom_point(aes(x = smd, colour = stage), size = 2.7) +
  geom_vline(xintercept = 0.1, linetype = 2, colour = "#555555") +
  scale_colour_manual(values = c("Before matching" = "#B07D62",
                                 "After matching" = "#C45A3C"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Absolute standardized mean difference", y = NULL,
       title = "C. Matching-balance audit",
       subtitle = "Dashed line = 0.10") +
  theme_pub(base_size = 11) +
  theme(legend.position = "top", axis.text.y = element_text(size = 9))

# ---- Panel D: alternate PPI estimator ----
overlap <- sensitivity %>% filter(method == "overlap weighting")
pD <- ggplot(overlap, aes(x = "Overlap-weighted\nmarginal OR", y = estimate)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.12, linewidth = 0.9,
                colour = "#C45A3C") +
  geom_point(size = 3.4, colour = "#8C2D2D") +
  geom_hline(yintercept = 1, linetype = 2, colour = "#555555") +
  scale_y_log10(breaks = c(0.5, 1, 2, 5, 10),
                expand = expansion(mult = c(0.12, 0.12))) +
  labs(x = NULL, y = "Marginal odds ratio (log scale)",
       title = "D. Alternative PPI estimator",
       subtitle = "OR 2.53, 95% CI 0.56-11.51, P = 0.230") +
  theme_pub(base_size = 11) +
  theme(axis.text.x = element_text(size = 9))

# ---- Panel E: mouse-layer detection stratification ----
smD <- jsonlite::fromJSON(file.path(SRC, "p0b_detectability_adjusted_summary.json"))
thr <- tibble(
  threshold = factor(c("Full detection",
                       "High-confidence detection\n(sequence coverage >= 38.6%)"),
                     levels = c("Full detection",
                                "High-confidence detection\n(sequence coverage >= 38.6%)")),
  fold    = c(smD$full_mouse_exposure$fold, smD$highconf_mouse_exposure$fold),
  p       = c(smD$full_mouse_exposure$fisher_p, smD$highconf_mouse_exposure$fisher_p),
  n       = c(smD$full_mouse_exposure$n_query, smD$highconf_mouse_exposure$n_query),
  query   = c(smD$full_mouse_exposure$rate_query, smD$highconf_mouse_exposure$rate_query),
  matched = c(smD$full_mouse_exposure$rate_matched, smD$highconf_mouse_exposure$rate_matched)
)
longE <- thr %>%
  pivot_longer(c(query, matched), names_to = "grp", values_to = "rate") %>%
  mutate(grp = factor(recode(grp, query = "CRISPR-hit query", matched = "Confounder-matched"),
                      levels = c("CRISPR-hit query", "Confounder-matched")))
pE <- ggplot(longE, aes(threshold, rate, fill = grp)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62,
           colour = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%", rate * 100)),
            position = position_dodge(width = 0.7), vjust = -0.4, size = 3.0) +
  geom_text(data = thr, aes(x = threshold, y = Inf, label = sprintf("fold = %.2fx", fold)),
            vjust = 1.2, size = 3.4, fontface = "bold", inherit.aes = FALSE) +
  scale_fill_manual(values = c("CRISPR-hit query" = "#C45A3C",
                               "Confounder-matched" = "#B07D62")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.34)),
                     labels = function(x) sprintf("%.0f%%", x * 100)) +
  labs(x = NULL, y = "P(filtered CRISPR evidence)", fill = NULL,
       title = "E. Coupling under high-confidence detection",
       caption = sprintf("Fisher P = %.1e (full detection, n = %d) and %.1e (high-confidence, n = %d).",
                         thr$p[1], thr$n[1], thr$p[2], thr$n[2])) +
  theme_pub(base_size = 11) +
  theme(legend.position = "top", axis.text.x = element_text(size = 9),
        plot.caption = element_text(hjust = 0, size = 7, colour = "black"))

# ---- Panel F: abundance dose-response ----
dr <- read_tsv(file.path(SRC, "p0b_abundance_dose_response.tsv"), show_col_types = FALSE)
dr <- dr %>% mutate(
  level_f = factor(level_label,
                   levels = c("Not detected", "LFQ tertile 1 (low)",
                              "LFQ tertile 2 (mid)", "LFQ tertile 3 (high)")),
  pct = crispr_clean_rate * 100)
pF <- ggplot(dr, aes(level_f, pct, fill = factor(level))) +
  geom_col(colour = "black", linewidth = 0.3, width = 0.68) +
  geom_text(aes(label = sprintf("%.1f%%\nn = %s", pct, format(n, big.mark = ","))),
            vjust = -0.35, size = 3.0) +
  scale_fill_manual(values = c("0" = "#B0B0B0", "1" = "#D9A066",
                               "2" = "#C45A3C", "3" = "#8C2D2D"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.28))) +
  labs(x = NULL, y = "CRISPR-clean rate (%)",
       title = "F. Abundance dose-response",
       caption = sprintf("Ordinal exposure (not detected = 0, LFQ tertiles = 1-3);\ntrend OR = %.2f per level, 95%% CI %.2f-%.2f, P = %.0e.",
                         smD$dose_response_trend_OR_per_level,
                         smD$dose_response_trend_CI[1], smD$dose_response_trend_CI[2],
                         smD$dose_response_trend_p)) +
  theme_pub(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 9),
        plot.caption = element_text(hjust = 0, size = 7, colour = "black"))

fig5 <- (pA + pB + pC + pD + pE + pF + plot_layout(ncol = 3)) & panel_labels_topleft()
save_pdf(fig5, "Figure5_evidence_coupling_main", w = 17, h = 10.6)
save_png(fig5, "Figure5_evidence_coupling_main", w = 17, h = 10.6)
write_tsv(coupling %>% select(label, fold_vs_matched, perm_p, log2fold),
          file.path(OUT, "Figure5A_source_data.tsv"))
write_tsv(joint_rate, file.path(OUT, "Figure5B_source_data.tsv"))
write_tsv(balance, file.path(OUT, "Figure5C_source_data.tsv"))
write_tsv(overlap, file.path(OUT, "Figure5D_source_data.tsv"))
write_tsv(longE, file.path(OUT, "Figure5E_source_data.tsv"))
write_tsv(dr, file.path(OUT, "Figure5F_source_data.tsv"))

# ---------------------------------------------------------------- #
# SUPPLEMENTARY FIGURE 8 (correlation-robust diagnostics)
# ---------------------------------------------------------------- #

# ---- Panel A: correlation-robust joint model ----
gl  <- read_tsv(file.path(SRC, "corr_robust_logistic.tsv"), show_col_types = FALSE)
lblE <- c(lcv = "LCV-associated detection", ppi = "Physical interaction",
          response = "Transcriptional response", nstudy = "Expression-study count",
          chronos = "U937 Chronos effect", common = "Common-essentiality")
gl2 <- gl %>% filter(term %in% names(lblE)) %>%
  mutate(term_label = factor(lblE[term], levels = rev(unname(lblE))))
longA <- bind_rows(
  gl2 %>% transmute(term_label, estimator = "Ordinary SE",
                    OR, lo = ci_lo_ordinary, hi = ci_hi_ordinary, p = p_ordinary),
  gl2 %>% transmute(term_label, estimator = "Module-clustered SE",
                    OR, lo = ci_lo_cluster, hi = ci_hi_cluster, p = p_cluster)
) %>% mutate(estimator = factor(estimator, levels = c("Ordinary SE", "Module-clustered SE")))
sA <- ggplot(longA, aes(x = OR, y = term_label, colour = estimator)) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y",
                width = 0.24, linewidth = 0.8, position = position_dodge(width = 0.55)) +
  geom_point(size = 2.8, position = position_dodge(width = 0.55)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "#555555") +
  scale_x_log10(expand = expansion(mult = c(0.08, 0.14))) +
  scale_colour_manual(values = c("Ordinary SE" = "#D9A066",
                                 "Module-clustered SE" = "#8C2D2D")) +
  labs(x = "Odds ratio (log scale)", y = NULL, colour = NULL,
       title = "A. Correlation-robust joint model") +
  theme_pub(base_size = 11) +
  theme(legend.position = "top", axis.text.y = element_text(size = 9))

# ---- Panel B: module-decimation matched null ----
sm_ <- jsonlite::fromJSON(file.path(SRC, "corr_robust_summary.json"))
dec <- tibble(
  median = sm_$B_decimation_fold_median,
  p05    = sm_$B_decimation_fold_p05,
  p95    = sm_$B_decimation_fold_p95,
  full   = sm_$B_full_universe_matched$fold,
  nvalid = sm_$B_decimation_n_valid,
  frac   = sm_$B_decimation_frac_fisher_p_lt_0.05
)
sB <- ggplot(dec) +
  geom_segment(aes(x = p05, xend = p95, y = 1, yend = 1),
               colour = "#4C72B0", linewidth = 1.8) +
  geom_point(aes(x = median, y = 1), size = 4.2, colour = "#4C72B0") +
  geom_vline(aes(xintercept = full), linetype = 2, colour = "#C45A3C", linewidth = 0.8) +
  geom_text(aes(x = median, y = 1.30, label = sprintf("median = %.1fx", median)),
            size = 3.3, fontface = "bold", colour = "#2F4A6D") +
  geom_text(aes(x = full, y = 0.72, label = sprintf("full universe = %.1fx", full)),
            size = 3.0, colour = "#C45A3C", hjust = -0.05) +
  scale_y_continuous(limits = c(0.55, 1.45)) +
  scale_x_continuous(expand = expansion(mult = c(0.14, 0.16))) +
  labs(x = "Matched fold enrichment (LCV -> filtered CRISPR)", y = NULL,
       title = "B. Module-decimation null",
       caption = sprintf("%d/%d decimation replicates remain significant (Fisher P < 0.05).\nBlue band = 5th-95th percentile of 200 replicates; point = median.",
                         round(dec$frac * dec$nvalid), dec$nvalid)) +
  theme_pub(base_size = 11) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.caption = element_text(hjust = 0, size = 7, colour = "black"))

# ---- Panel C: within-detection abundance gradient ----
gr  <- read_tsv(file.path(SRC, "p0b_detectability_gradient.tsv"), show_col_types = FALSE)
glab <- c(log10_LFQ = "log10 LFQ intensity", seq_coverage_pct = "Sequence coverage (%)",
          "_deg" = "Network degree", "_pmids" = "Publication count",
          "_chronos" = "U937 Chronos effect", "_common_ess" = "Common-essentiality",
          "_ui_expr" = "Baseline expression")
gr2 <- gr %>%
  filter(term != "const", term %in% names(glab)) %>%
  mutate(label = factor(glab[term], levels = rev(unname(glab))),
         sig = p < 0.05)
sC <- ggplot(gr2, aes(x = OR, y = label, colour = sig)) +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), orientation = "y",
                width = 0.24, linewidth = 0.8) +
  geom_point(size = 2.8) +
  geom_vline(xintercept = 1, linetype = 2, colour = "#555555") +
  scale_x_log10(expand = expansion(mult = c(0.10, 0.14))) +
  scale_colour_manual(values = c(`TRUE` = "#C45A3C", `FALSE` = "#9E9E9E"),
                      labels = c(`TRUE` = "P < 0.05", `FALSE` = "not significant")) +
  labs(x = "Odds ratio per unit (log scale)", y = NULL, colour = NULL,
       title = "C. Within-detection abundance gradient",
       caption = "Logistic model of CRISPR-hit status within the 2,072 mouse-detected proteins;\nabundance and coverage are modelled jointly with the audited covariates.") +
  theme_pub(base_size = 11) +
  theme(legend.position = "top", axis.text.y = element_text(size = 9),
        plot.caption = element_text(hjust = 0, size = 7, colour = "black"))

figS8 <- (sA + sB + sC + plot_layout(ncol = 3)) & panel_labels_topleft()
save_pdf(figS8, "SupplementaryFigure8_coupling_diagnostics", w = 17, h = 5.4)
save_png(figS8, "SupplementaryFigure8_coupling_diagnostics", w = 17, h = 5.4)
write_tsv(longA, file.path(OUT, "SuppFigure8A_source_data.tsv"))
write_tsv(dec, file.path(OUT, "SuppFigure8B_source_data.tsv"))
write_tsv(gr2, file.path(OUT, "SuppFigure8C_source_data.tsv"))

message("done")
