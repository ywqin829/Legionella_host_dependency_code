# Figure 5. Evidence modalities are non-equivalent; post-freeze audits retain
# transparent uncertainty and show where association is estimator-sensitive.
#
# Second-round editorial round (2026-10-04): the two new robustness analyses
# (correlation-robust coupling, mouse-layer detectability) are merged here as
# panels E-H, and the abundance dose-response is added as panel I. Panels A-D
# are unchanged from the previous version. Layout is 3 x 3.
ROOT <- Sys.getenv("LEGIONELLA_ROOT")
if (!nzchar(ROOT)) stop("Set LEGIONELLA_ROOT before running this script.")
if (!dir.exists(ROOT)) stop("LEGIONELLA_ROOT does not exist: ", ROOT)
setwd(file.path(ROOT, "15_manuscript/r_figures"))
source(file.path(ROOT, "15_manuscript/r_figures/_theme.R"))
suppressMessages({library(dplyr); library(readr); library(tidyr); library(patchwork); library(jsonlite)})
SRC <- file.path(ROOT, "13_postfreeze_interpretation")

coupling <- read_tsv(file.path(SRC, "a1_evidence_coupling_matched.tsv"), show_col_types = FALSE)
joint_rate <- read_tsv(file.path(SRC, "g_ppi_lcv_exact_rate_ci.tsv"), show_col_types = FALSE)
balance <- read_tsv(file.path(SRC, "g_matching_balance_by_covariate.tsv"), show_col_types = FALSE)
sensitivity <- read_tsv(file.path(SRC, "g_matching_sensitivity.tsv"), show_col_types = FALSE)

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
       title = "B. Joint evidence with exact 95% CI") +
  theme_pub(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 9))

# ---- Panel C: covariate balance with publication-ready labels ----
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
pC <- ggplot(balance, aes(y = covariate_label)) +
  geom_segment(aes(x = abs_smd_before, xend = abs_smd_after, yend = covariate_label),
               colour = "#B07D62", linewidth = 0.8) +
  geom_point(aes(x = abs_smd_before), colour = "#B07D62", size = 2.7) +
  geom_point(aes(x = abs_smd_after), colour = "#C45A3C", size = 2.7) +
  geom_vline(xintercept = 0.1, linetype = 2, colour = "#555555") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Absolute standardized mean difference", y = NULL,
       title = "C. Matching-balance audit",
       subtitle = "Before / after matching; dashed line = 0.10") +
  theme_pub(base_size = 11) +
  theme(axis.text.y = element_text(size = 9))

# ---- Panel D: alternate estimator preserves uncertainty and weak evidence ----
overlap <- sensitivity %>% filter(method == "overlap weighting")
pD <- ggplot(overlap, aes(x = estimate, y = method)) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high, ymin = method, ymax = method),
                 orientation = "y", width = 0.18, linewidth = 0.8,
                 colour = "#C45A3C") +
  geom_point(size = 3.2, colour = "#8C2D2D") +
  geom_vline(xintercept = 1, linetype = 2, colour = "#555555") +
  scale_x_log10(expand = expansion(mult = c(0.08, 0.12))) +
  labs(x = "Marginal odds ratio (log scale)", y = "Estimator",
       title = "D. Alternative PPI-to-dependency estimator",
       subtitle = "Overlap weighting; OR 2.53, P = 0.230") +
  theme_pub(base_size = 11) +
  theme(axis.text.y = element_text(size = 9), axis.ticks.y = element_line())

# ---- Panel E: correlation-robust joint model (ordinary vs module-clustered SE) ----
gl  <- read_tsv(file.path(SRC, "corr_robust_logistic.tsv"), show_col_types = FALSE)
lblE <- c(lcv = "LCV-associated detection", ppi = "Physical interaction",
          response = "Transcriptional response", nstudy = "Expression-study count",
          chronos = "U937 Chronos effect", common = "Common-essentiality")
gl2 <- gl %>% filter(term %in% names(lblE)) %>%
  mutate(term_label = factor(lblE[term], levels = rev(unname(lblE))))
longE <- bind_rows(
  gl2 %>% transmute(term_label, estimator = "Ordinary SE",
                    OR, lo = ci_lo_ordinary, hi = ci_hi_ordinary, p = p_ordinary),
  gl2 %>% transmute(term_label, estimator = "Module-clustered SE",
                    OR, lo = ci_lo_cluster, hi = ci_hi_cluster, p = p_cluster)
) %>% mutate(estimator = factor(estimator, levels = c("Ordinary SE", "Module-clustered SE")))
pE <- ggplot(longE, aes(x = OR, y = term_label, colour = estimator)) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y",
                width = 0.24, linewidth = 0.8, position = position_dodge(width = 0.55)) +
  geom_point(size = 2.8, position = position_dodge(width = 0.55)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "#555555") +
  scale_x_log10(expand = expansion(mult = c(0.08, 0.14))) +
  scale_colour_manual(values = c("Ordinary SE" = "#D9A066",
                                 "Module-clustered SE" = "#8C2D2D")) +
  labs(x = "Odds ratio (log scale)", y = NULL, colour = NULL,
       title = "E. Correlation-robust joint model") +
  theme_pub(base_size = 11) +
  theme(legend.position = "top", axis.text.y = element_text(size = 9))

# ---- Panel F: module-decimation matched null ----
sm_ <- jsonlite::fromJSON(file.path(SRC, "corr_robust_summary.json"))
dec <- tibble(
  median = sm_$B_decimation_fold_median,
  p05    = sm_$B_decimation_fold_p05,
  p95    = sm_$B_decimation_fold_p95,
  full   = sm_$B_full_universe_matched$fold,
  nvalid = sm_$B_decimation_n_valid,
  frac   = sm_$B_decimation_frac_fisher_p_lt_0.05
)
pF <- ggplot(dec) +
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
       title = "F. Module-decimation null",
       caption = sprintf("%d/%d decimation replicates remain significant (Fisher P < 0.05).\nBlue band = 5th-95th percentile of 200 replicates; point = median.",
                         round(dec$frac * dec$nvalid), dec$nvalid)) +
  theme_pub(base_size = 11) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.caption = element_text(hjust = 0, size = 7, colour = "black"))

# ---- Panel G: coupling under full vs high-confidence detection ----
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
longG <- thr %>%
  pivot_longer(c(query, matched), names_to = "grp", values_to = "rate") %>%
  mutate(grp = factor(recode(grp, query = "CRISPR-hit query", matched = "Confounder-matched"),
                      levels = c("CRISPR-hit query", "Confounder-matched")))
pG <- ggplot(longG, aes(threshold, rate, fill = grp)) +
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
       title = "G. Coupling under high-confidence detection",
       caption = sprintf("Fisher P = %.1e (full detection, n = %d) and %.1e (high-confidence, n = %d).",
                         thr$p[1], thr$n[1], thr$p[2], thr$n[2])) +
  theme_pub(base_size = 11) +
  theme(legend.position = "top", axis.text.x = element_text(size = 9),
        plot.caption = element_text(hjust = 0, size = 7, colour = "black"))

# ---- Panel H: residual abundance gradient within mouse detections ----
gr  <- read_tsv(file.path(SRC, "p0b_detectability_gradient.tsv"), show_col_types = FALSE)
glab <- c(log10_LFQ = "log10 LFQ intensity", seq_coverage_pct = "Sequence coverage (%)",
          "_deg" = "Network degree", "_pmids" = "Publication count",
          "_chronos" = "U937 Chronos effect", "_common_ess" = "Common-essentiality",
          "_ui_expr" = "Baseline expression")
gr2 <- gr %>%
  filter(term != "const", term %in% names(glab)) %>%
  mutate(label = factor(glab[term], levels = rev(unname(glab))),
         sig = p < 0.05)
pH <- ggplot(gr2, aes(x = OR, y = label, colour = sig)) +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), orientation = "y",
                width = 0.24, linewidth = 0.8) +
  geom_point(size = 2.8) +
  geom_vline(xintercept = 1, linetype = 2, colour = "#555555") +
  scale_x_log10(expand = expansion(mult = c(0.10, 0.14))) +
  scale_colour_manual(values = c(`TRUE` = "#C45A3C", `FALSE` = "#9E9E9E"),
                      labels = c(`TRUE` = "P < 0.05", `FALSE` = "not significant")) +
  labs(x = "Odds ratio per unit (log scale)", y = NULL, colour = NULL,
       title = "H. Within-detection abundance gradient",
       caption = "Logistic model of CRISPR-hit status within the 2,072 mouse-detected proteins;\nabundance and coverage are modelled jointly with the audited covariates.") +
  theme_pub(base_size = 11) +
  theme(legend.position = "top", axis.text.y = element_text(size = 9),
        plot.caption = element_text(hjust = 0, size = 7, colour = "black"))

# ---- Panel I: ordinal abundance dose-response of CRISPR-hit status ----
dr <- read_tsv(file.path(SRC, "p0b_abundance_dose_response.tsv"), show_col_types = FALSE)
dr <- dr %>% mutate(
  level_f = factor(level_label,
                   levels = c("Not detected", "LFQ tertile 1 (low)",
                              "LFQ tertile 2 (mid)", "LFQ tertile 3 (high)")),
  pct = crispr_clean_rate * 100)
pI <- ggplot(dr, aes(level_f, pct, fill = factor(level))) +
  geom_col(colour = "black", linewidth = 0.3, width = 0.68) +
  geom_text(aes(label = sprintf("%.1f%%\nn = %s", pct, format(n, big.mark = ","))),
            vjust = -0.35, size = 3.0) +
  scale_fill_manual(values = c("0" = "#B0B0B0", "1" = "#D9A066",
                               "2" = "#C45A3C", "3" = "#8C2D2D"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.28))) +
  labs(x = NULL, y = "CRISPR-clean rate (%)",
       title = "I. Abundance dose-response",
       caption = sprintf("Ordinal exposure (not detected = 0, LFQ tertiles = 1-3);\ntrend OR = %.2f per level, 95%% CI %.2f-%.2f, P = %.0e.",
                         smD$dose_response_trend_OR_per_level,
                         smD$dose_response_trend_CI[1], smD$dose_response_trend_CI[2],
                         smD$dose_response_trend_p)) +
  theme_pub(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 9),
        plot.caption = element_text(hjust = 0, size = 7, colour = "black"))

fig <- (pA + pB + pC + pD + pE + pF + pG + pH + pI + plot_layout(ncol = 3)) & panel_labels_topleft()
save_pub_pdf(fig, "Figure5_evidence_coupling", w = 17, h = 15)
