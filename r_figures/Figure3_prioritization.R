# Figure 3. Multimodal prioritization and internal positive-control recovery
ROOT <- Sys.getenv("LEGIONELLA_ROOT")
if (!nzchar(ROOT)) stop("Set LEGIONELLA_ROOT before running this script.")
if (!dir.exists(ROOT)) stop("LEGIONELLA_ROOT does not exist: ", ROOT)
setwd(file.path(ROOT, "15_manuscript/r_figures"))
source(file.path(ROOT, "15_manuscript/r_figures/_theme.R"))
SRC <- file.path(ROOT, "11_figures/source_data")
suppressMessages({library(dplyr); library(tidyr); library(readr); library(patchwork)})

# The candidate evidence matrix is provided in Table 1/source data rather than
# repeated here. The main figure therefore has two panels, labelled A and B.
tier_counts <- read_tsv(file.path(SRC, "Figure3A_final_tier_counts.tsv"), show_col_types = FALSE)
positive_controls <- read_tsv(file.path(SRC, "Figure3B_positive_control_recovery.tsv"), show_col_types = FALSE) %>%
  rename(gene = human_symbol)

# ---- Panel A: final tier counts ----
tier_counts <- tier_counts %>%
  mutate(
    tier_label = case_when(
      grepl("TIER1", final_tier) ~ "Tier 1",
      grepl("TIER2", final_tier) ~ "Tier 2",
      grepl("TIER3_CONTEXT", final_tier) ~ "Tier 3 context",
      grepl("FITNESS", final_tier) ~ "Tier 3 fitness",
      TRUE ~ "Rejected"
    ),
    tier_label = factor(tier_label, levels = rev(unique(tier_label)))
  )

# Reserve right-side space for the numeric label above the longest horizontal bar.
count_limit <- max(tier_counts$n) + max(3, ceiling(max(tier_counts$n) * 0.15))

pA <- ggplot(tier_counts, aes(x = tier_label, y = n, fill = tier_label)) +
  geom_col(colour = "black", linewidth = 0.3, width = 0.7) +
  geom_text(aes(label = n), hjust = -0.15, size = 3.5) +
  coord_flip(clip = "off") +
  scale_y_continuous(
    limits = c(0, count_limit),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_fill_manual(values = c(
    "Tier 1" = "#C1121F",
    "Tier 2" = "#E07A1F",
    "Tier 3 context" = "#F4D03F",
    "Tier 3 fitness" = "#8D99AE",
    "Rejected" = "#333333"
  )) +
  labs(x = NULL, y = "Candidates", fill = NULL, title = "A. Final strict tiers") +
  theme_pub() +
  theme(
    legend.position = "none",
    plot.margin = margin(t = 5.5, r = 22, b = 5.5, l = 5.5)
  )

# ---- Panel B: positive-control recovery ----
positive_controls_long <- positive_controls %>%
  pivot_longer(-gene, names_to = "domain", values_to = "recovered")

pB <- ggplot(positive_controls_long, aes(domain, gene, fill = factor(recovered))) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = c("0" = "#F5F5F5", "1" = "#2E8B57"), name = "recovered") +
  labs(x = NULL, y = NULL, title = "B. Positive-control recovery (9/9)") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

fig <- (pA | pB) & panel_labels_topleft()
save_pub_pdf(fig, "Figure3_prioritization", w = 9, h = 4.8)
