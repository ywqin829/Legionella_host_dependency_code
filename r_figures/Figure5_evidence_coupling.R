# Figure 5. Evidence modalities are non-equivalent; post-freeze audits retain
# transparent uncertainty and show where association is estimator-sensitive.
ROOT <- Sys.getenv("LEGIONELLA_ROOT")
if (!nzchar(ROOT)) ROOT <- Sys.getenv("LEGIONELLA_ROOT")
if (!nzchar(ROOT)) stop("Set LEGIONELLA_ROOT before running this script.")
if (!dir.exists(ROOT)) stop("LEGIONELLA_ROOT does not exist: ", ROOT)
setwd(file.path(ROOT, "15_manuscript/r_figures"))
source(file.path(ROOT, "15_manuscript/r_figures/_theme.R"))
suppressMessages({library(dplyr); library(readr); library(patchwork)})
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

fig <- (pA | pB) / (pC | pD)
save_pub_pdf(fig, "Figure5_evidence_coupling", w = 13, h = 9.5)
