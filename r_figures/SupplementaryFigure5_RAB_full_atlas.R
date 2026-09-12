# Supplementary Figure 5. Full RAB GTPase family evidence atlas (n=61)
ROOT <- Sys.getenv("LEGIONELLA_ROOT")
if (!nzchar(ROOT)) ROOT <- Sys.getenv("LEGIONELLA_ROOT")
if (!nzchar(ROOT)) stop("Set LEGIONELLA_ROOT before running this script.")
if (!dir.exists(ROOT)) stop("LEGIONELLA_ROOT does not exist: ", ROOT)
setwd(file.path(ROOT, "15_manuscript/r_figures"))
source(file.path(ROOT, "15_manuscript/r_figures/_theme.R"))
SRC <- file.path(ROOT, "13_postfreeze_interpretation")
suppressMessages({library(dplyr); library(tidyr); library(readr)})

TR <- read_tsv(file.path(SRC, "a3_trafficking_gene_evidence_matrix.tsv"))
rab <- TR %>% filter(family == "RAB") %>%
  mutate(n_ev = as.integer(ppi_AB) + as.integer(apms_high) + as.integer(crispr_hit) +
           as.integer(crispr_clean) + as.integer(stable_cross_cohort) + as.integer(lcv_spatial)) %>%
  arrange(desc(n_ev)) %>% mutate(gene = factor(gene, levels = rev(gene)))
tracks <- c("ppi_AB", "apms_high", "crispr_hit", "crispr_clean", "stable_cross_cohort", "lcv_spatial")
tlab <- c("PPI(A/B)", "AP-MS", "CRISPR", "Filtered", "Response", "LCV")

rab_long <- rab %>% select(gene, all_of(tracks)) %>%
  pivot_longer(-gene, names_to = "track", values_to = "val") %>%
  mutate(val = as.integer(val))
p <- ggplot(rab_long, aes(track, gene, fill = factor(val))) +
  geom_tile(colour = "white", linewidth = 0.5) +
  scale_fill_manual(values = c("0" = "#F5F5F5", "1" = "#D81B60"), name = "Evidence present") +
  scale_x_discrete(labels = tlab) +
  labs(x = "Evidence modality", y = "RAB GTPase gene", title = "Full RAB GTPase family evidence atlas (n = 61)") +
  theme_pub(base_size = 11) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                      axis.text.y = element_text(size = 9))

save_pub_pdf(p, "SupplementaryFigure5_RAB_full_atlas", w = 8, h = 12)
