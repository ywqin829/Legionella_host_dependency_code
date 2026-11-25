# =====================================================================
# Figure 5 - print-size redraw (6.75 in wide).
# Panel content and source data are unchanged from figure5_split.R.
# Only canvas size, font sizes, in-panel text and layout differ.
# =====================================================================
if (!exists("PRINT_W")) {
  .dir <- (function() {
    a <- commandArgs(trailingOnly = FALSE)
    f <- grep("^--file=", a, value = TRUE)
    if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1]))))
    getwd()
  })()
  source(file.path(.dir, "_print_common.R"))
}

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
    grepl("PPI\\(A/B\\) -> clean", label) ~ "PPI(A/B) -> filt. CRISPR",
    grepl("PPI\\(any\\) -> clean", label) ~ "PPI(any) -> filt. CRISPR",
    grepl("LCV -> CRISPR hit", label) ~ "LCV -> CRISPR hit",
    grepl("LCV -> clean", label) ~ "LCV -> filt. CRISPR",
    grepl("PPI\\(A/B\\) -> LCV", label) ~ "PPI(A/B) -> LCV",
    grepl("CRISPR hit -> response", label) ~ "CRISPR hit -> response",
    grepl("clean CRISPR -> response", label) ~ "filt. CRISPR -> response",
    grepl("response -> clean", label) ~ "response -> filt. CRISPR",
    grepl("LCV -> response", label) ~ "LCV -> response",
    grepl("PPI\\(A/B\\) -> response", label) ~ "PPI(A/B) -> response",
    TRUE ~ label
  )
) %>% mutate(short = factor(short, levels = rev(unique(short))))

pA <- ggplot(coupling, aes(log2fold, short, fill = log2fold >= 1.25)) +
  geom_col(colour = "black", linewidth = 0.25, width = 0.68) +
  geom_text(aes(label = star), hjust = -0.25, size = 2.2) +
  geom_vline(xintercept = 0, linewidth = 0.35) +
  scale_x_continuous(expand = expansion(mult = c(0.10, 0.18))) +
  scale_fill_manual(values = c("TRUE" = "#C45A3C", "FALSE" = "#D9A066"), guide = "none") +
  labs(x = "log2 fold-enrichment vs matched background", y = NULL,
       title = "A. Evidence coupling") +
  theme_print(7) +
  theme(axis.text.y = element_text(size = 6.2),
        plot.margin = margin(2, 10, 7, 2))

# ---- Panel B: exact uncertainty for the joint evidence rate ----
joint_rate$group <- factor(joint_rate$group,
                           levels = c("Universe", "PPI(A/B)", "LCV", "PPI(A/B)+LCV"))
pB <- ggplot(joint_rate, aes(group, rate, fill = group)) +
  geom_col(colour = "black", linewidth = 0.25, width = 0.66) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.12, linewidth = 0.45) +
  geom_text(aes(label = paste0(successes, "/", n)), vjust = -0.5, size = 2.1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  scale_fill_manual(values = c("Universe" = "#B07D62", "PPI(A/B)" = "#D9A066",
                               "LCV" = "#C45A3C", "PPI(A/B)+LCV" = "#8C2D2D"),
                    guide = "none") +
  labs(x = NULL, y = "P(filtered CRISPR evidence | evidence)",
       title = "B. Joint evidence (exact 95% CI)") +
  theme_print(7) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 6.0))

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
               colour = "#B07D62", linewidth = 0.5) +
  geom_point(aes(x = smd, colour = stage), size = 1.9) +
  geom_vline(xintercept = 0.1, linetype = 2, colour = "#555555", linewidth = 0.35) +
  scale_colour_manual(values = c("Before matching" = "#B07D62",
                                 "After matching" = "#C45A3C"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Absolute standardized mean difference", y = NULL,
       title = "C. Matching-balance audit") +
  theme_print(7) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.margin = margin(0, 0, 0, 0),
        axis.text.y = element_text(size = 6.3))

# ---- Panel D: alternate PPI estimator ----
overlap <- sensitivity %>% filter(method == "overlap weighting")
pD <- ggplot(overlap, aes(x = "Overlap-weighted\nmarginal OR", y = estimate)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.10, linewidth = 0.7,
                colour = "#C45A3C") +
  geom_point(size = 2.4, colour = "#8C2D2D") +
  geom_hline(yintercept = 1, linetype = 2, colour = "#555555", linewidth = 0.35) +
  scale_y_log10(breaks = c(0.5, 1, 2, 5, 10),
                expand = expansion(mult = c(0.12, 0.12))) +
  labs(x = NULL, y = "Marginal odds ratio (log scale)",
       title = "D. Alternative PPI estimator") +
  theme_print(7) + theme(axis.text.x = element_text(size = 6.2))

# ---- Panel E: mouse-layer detection stratification ----
smD <- jsonlite::fromJSON(file.path(SRC, "p0b_detectability_adjusted_summary.json"))
thr <- tibble(
  threshold = factor(c("Full detection",
                       "High-confidence\n(coverage >= 38.6%)"),
                     levels = c("Full detection",
                                "High-confidence\n(coverage >= 38.6%)")),
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
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           colour = "black", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.1f%%", rate * 100)),
            position = position_dodge(width = 0.7), vjust = -0.35, size = 2.0) +
  geom_text(data = thr, aes(x = threshold, y = Inf, label = sprintf("fold = %.2fx", fold)),
            vjust = 1.2, size = 2.2, fontface = "bold", inherit.aes = FALSE) +
  scale_fill_manual(values = c("CRISPR-hit query" = "#C45A3C",
                               "Confounder-matched" = "#B07D62")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.34)),
                     labels = function(x) sprintf("%.0f%%", x * 100)) +
  labs(x = NULL, y = "P(filtered CRISPR evidence)", fill = NULL,
       title = "E. High-confidence detection") +
  theme_print(7) +
  theme(legend.position = "top", axis.text.x = element_text(size = 5.8, angle = 15, hjust = 1))

# ---- Panel F: abundance dose-response ----
dr <- read_tsv(file.path(SRC, "p0b_abundance_dose_response.tsv"), show_col_types = FALSE)
dr <- dr %>% mutate(
  level_f = factor(level_label,
                   levels = c("Not detected", "LFQ tertile 1 (low)",
                              "LFQ tertile 2 (mid)", "LFQ tertile 3 (high)")),
  pct = crispr_clean_rate * 100)
pF <- ggplot(dr, aes(level_f, pct, fill = factor(level))) +
  geom_col(colour = "black", linewidth = 0.25, width = 0.66) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.35, size = 2.0) +
  scale_fill_manual(values = c("0" = "#B0B0B0", "1" = "#D9A066",
                               "2" = "#C45A3C", "3" = "#8C2D2D"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.26))) +
  labs(x = NULL, y = "CRISPR-clean rate (%)",
       title = "F. Abundance dose-response") +
  theme_print(7) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 5.8))

fig5 <- (pA | pB) / ((pC | pD) + plot_layout(widths = c(1.45, 1))) / (pE | pF) +
  plot_layout(heights = c(1.25, 1, 1))
fig5 <- fig5 & panel_labels_topleft()
save_print(fig5, "Figure5_evidence_coupling_print", w = 6.75, h = 7.6)
