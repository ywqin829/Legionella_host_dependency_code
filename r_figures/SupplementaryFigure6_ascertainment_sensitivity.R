library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)

ROOT <- Sys.getenv("LEGIONELLA_ROOT")
if (!nzchar(ROOT)) ROOT <- Sys.getenv("LEGIONELLA_ROOT")
if (!nzchar(ROOT)) stop("Set LEGIONELLA_ROOT before running this script.")
if (!dir.exists(ROOT)) stop("LEGIONELLA_ROOT does not exist: ", ROOT)
OUT <- file.path(ROOT, "15_manuscript", "r_figures")
DIR <- file.path(ROOT, "13_postfreeze_interpretation")

enr <- read_tsv(file.path(DIR, "a4_panel_enrichment.tsv"), show_col_types = FALSE) %>%
  mutate(category = if_else(grepl("RAB", label), "RAB family", "Trafficking universe"),
         universe = case_when(grepl("crispr_tested", label) ~ "CRISPR-tested",
                              TRUE ~ "Integrated atlas"))

null <- jsonlite::fromJSON(file.path(DIR, "f7_selection_aware_enrichment.json"))
null_df <- tibble(category = c("RAB family", "Trafficking universe"),
                  observed = c(null$obs_rab_in_panel, null$obs_trafficking_in_panel),
                  median = c(null$null_rab_median, null$null_traf_median),
                  p95 = c(null$null_rab_p95, null$null_traf_p95),
                  maximum = c(null$null_rab_max, null$null_traf_max),
                  empirical_p = c(null$bonferroni_adjusted_p_rab, null$bonferroni_adjusted_p_trafficking))

p1 <- ggplot(enr, aes(x = universe, y = fold, fill = category)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62, colour = "black") +
  geom_text(aes(label = sprintf("%.1fx", fold)), position = position_dodge(width = 0.72),
            vjust = -0.35, size = 3.2, fontface = "bold") +
  scale_fill_manual(values = c("RAB family" = "#C44E52", "Trafficking universe" = "#4C72B0")) +
  labs(x = NULL, y = "Fold enrichment", fill = NULL) +
  theme_classic(base_size = 10) + theme(legend.position = "top", axis.text.x = element_text(angle = 15, hjust = 1))

p2 <- ggplot(null_df, aes(x = category, y = observed)) +
  geom_segment(aes(xend = category, y = median, yend = p95), linewidth = 1.3, colour = "#4C72B0") +
  geom_point(aes(y = median), size = 3, colour = "#4C72B0") +
  geom_point(aes(y = observed), size = 4, shape = 21, fill = "#C44E52", colour = "black") +
  geom_text(aes(y = observed, label = paste0("observed=", observed)), vjust = -0.8, size = 3.0) +
  labs(x = NULL, y = "Number of genes in 54-gene panel",
       caption = "Blue point = null median; vertical line = null median to 95th percentile; red point = observed panel count.\nBonferroni-adjusted empirical P = 2.0 × 10^-4 for each prespecified enrichment test.") +
  theme_classic(base_size = 10) + theme(axis.text.x = element_text(angle = 15, hjust = 1), plot.caption = element_text(hjust = 0, size = 8))

pdf(file.path(OUT, "SupplementaryFigure6_ascertainment_sensitivity.pdf"), width = 8.5, height = 4.6, useDingbats = FALSE)
gridExtra::grid.arrange(p1, p2, ncol = 2, top = "Background-universe and ascertainment-matched sensitivity")
dev.off()

write_tsv(enr, file.path(DIR, "SupplementaryFigure6_universe_enrichment_source.tsv"))
write_tsv(null_df, file.path(DIR, "SupplementaryFigure6_selection_null_source.tsv"))
