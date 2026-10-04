library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)

ROOT <- Sys.getenv("LEGIONELLA_ROOT")
if (!nzchar(ROOT)) stop("Set LEGIONELLA_ROOT before running this script.")
if (!dir.exists(ROOT)) stop("LEGIONELLA_ROOT does not exist: ", ROOT)
source(file.path(ROOT, "15_manuscript/r_figures/_theme.R"))
suppressMessages(library(patchwork))
OUT <- file.path(ROOT, "15_manuscript", "r_figures")
DIR <- file.path(ROOT, "13_postfreeze_interpretation")

enr <- read_tsv(file.path(DIR, "a4_panel_enrichment.tsv"), show_col_types = FALSE) %>%
  mutate(category = recode(label,
      "RAB_GTPase_family | panel | atlas" = "RAB GTPase (n=61)",
      "RAB_GTPase_family | panel | crispr_tested" = "RAB GTPase (n=61)",
      "RAB_associated_machinery | panel | atlas" = "RAB-associated (n=77)",
      "RAB_associated_machinery | panel | crispr_tested" = "RAB-associated (n=77)",
      "trafficking_universe | panel | atlas" = "Trafficking (n=275)",
      "trafficking_universe | panel | crispr_tested" = "Trafficking (n=275)"),
    category = factor(category, levels = c("RAB GTPase (n=61)", "RAB-associated (n=77)", "Trafficking (n=275)")),
    universe = case_when(grepl("crispr_tested", label) ~ "CRISPR-tested",
                         TRUE ~ "Integrated atlas"))

null <- jsonlite::fromJSON(file.path(DIR, "f7_selection_aware_enrichment.json"))
null_df <- tibble(category = c("RAB GTPase (n=61)", "Trafficking (n=275)"),
                  observed = c(null$obs_rab_in_panel, null$obs_trafficking_in_panel),
                  median = c(null$null_rab_median, null$null_traf_median),
                  p95 = c(null$null_rab_p95, null$null_traf_p95),
                  maximum = c(null$null_rab_max, null$null_traf_max),
                  empirical_p = c(null$bonferroni_adjusted_p_rab, null$bonferroni_adjusted_p_trafficking))

pA <- ggplot(enr, aes(x = universe, y = fold, fill = category)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62, colour = "black") +
  geom_text(aes(label = sprintf("%.1fx", fold)), position = position_dodge(width = 0.72),
            vjust = -0.35, size = 3.0, fontface = "bold") +
  scale_fill_manual(values = c("RAB GTPase (n=61)" = "#C44E52",
                               "RAB-associated (n=77)" = "#DD8452",
                               "Trafficking (n=275)" = "#4C72B0")) +
  labs(x = NULL, y = "Fold enrichment", fill = NULL,
       title = "A. Background-universe sensitivity") +
  theme_pub() +
  theme(legend.position = "top", axis.text.x = element_text(angle = 15, hjust = 1))

pB <- ggplot(null_df, aes(x = category, y = observed)) +
  geom_segment(aes(xend = category, y = median, yend = p95), linewidth = 1.3, colour = "#4C72B0") +
  geom_point(aes(y = median), size = 3, colour = "#4C72B0") +
  geom_point(aes(y = observed), size = 4, shape = 21, fill = "#C44E52", colour = "black") +
  geom_text(aes(y = observed, label = paste0("observed=", observed)), vjust = -0.8, size = 3.0) +
  labs(x = NULL, y = "Number of genes in 56-gene panel",
       title = "B. Selection-aware matched null",
       caption = "Blue point = null median; vertical line = null median to 95th percentile; red point = observed panel count.\nBonferroni-adjusted empirical P = 2.0 x 10^-4 for each prespecified enrichment test.") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 15, hjust = 1),
        plot.caption = element_text(hjust = 0, size = 7, colour = "black"))

fig <- (pA | pB) & panel_labels_topleft()
save_pub_pdf(fig, "SupplementaryFigure6_ascertainment_sensitivity", w = 9, h = 5.2)

write_tsv(enr, file.path(DIR, "SupplementaryFigure6_universe_enrichment_source.tsv"))
write_tsv(null_df, file.path(DIR, "SupplementaryFigure6_selection_null_source.tsv"))
