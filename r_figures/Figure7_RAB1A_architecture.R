# Figure 7. RAB1A: the sole Tier-1 candidate
# A: direct effector manipulation (literature-confirmed edges)
# B: host trafficking program evidence matrix (final tiers)
ROOT <- Sys.getenv("LEGIONELLA_ROOT")
if (!nzchar(ROOT)) stop("Set LEGIONELLA_ROOT before running this script.")
if (!dir.exists(ROOT)) stop("LEGIONELLA_ROOT does not exist: ", ROOT)
setwd(file.path(ROOT, "15_manuscript/r_figures"))
source(file.path(ROOT, "15_manuscript/r_figures/_theme.R"))
SRC <- file.path(ROOT, "13_postfreeze_interpretation")
suppressMessages({library(dplyr); library(tidyr); library(readr); library(patchwork)})

edges <- read_tsv(file.path(SRC, "fig7_effector_rab1a_edges.tsv"), show_col_types = FALSE)
host_evidence <- read_tsv(file.path(SRC, "fig7_host_trafficking_evidence_finaltier.tsv"), show_col_types = FALSE)

# ---- Panel A: effector -> RAB1A network (manual layout) ----
effectors <- edges %>% mutate(
  x = 0.18,
  y = seq(0.86, 0.14, length.out = n()),
  role_short = gsub(" \\|.*", "", role)
)

pA <- ggplot() +
  geom_segment(
    data = effectors,
    aes(x = x + 0.05, xend = 0.70, y = y, yend = 0.50),
    colour = "#D35400", linewidth = 0.8, alpha = 0.8
  ) +
  geom_point(data = effectors, aes(x, y), size = 5.5, colour = "#D35400") +
  # Long effector names are placed outside the nodes so every character remains visible.
  geom_text(
    data = effectors,
    aes(x = x - 0.035, y = y, label = effector),
    size = 3.5, hjust = 1, colour = "#6B2E0D", fontface = "bold"
  ) +
  geom_text(
    data = effectors,
    aes(x = (x + 0.70) / 2, y = (y + 0.50) / 2 + 0.025, label = role_short),
    size = 3.0, colour = "#7B3F00"
  ) +
  annotate("point", x = 0.70, y = 0.50, size = 18, colour = "#C1121F") +
  annotate("text", x = 0.70, y = 0.535, label = "RAB1A", size = 3.2,
           colour = "white", fontface = "bold") +
  annotate("text", x = 0.70, y = 0.465, label = "TIER 1", size = 2.9,
           colour = "white", fontface = "bold") +
  coord_cartesian(xlim = c(-0.03, 0.86), ylim = c(0.04, 0.96), clip = "off") +
  labs(title = "A. Direct effector manipulation of RAB1A") +
  theme_pub(base_size = 11) +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(),
    axis.title = element_blank(), axis.line = element_blank(),
    plot.margin = margin(t = 5.5, r = 8, b = 5.5, l = 24)
  )

# ---- Panel B: compact host trafficking evidence matrix ----
host_evidence <- host_evidence %>% mutate(
  tier_label = case_when(
    grepl("TIER1", tier) ~ "T1",
    grepl("TIER2", tier) ~ "T2",
    grepl("TIER3", tier) ~ "T3",
    TRUE ~ "Other"
  ),
  gene_label = paste0(gene, " (", tier_label, ")"),
  gene_label = factor(gene_label, levels = rev(unique(gene_label)))
)
tracks <- c("PPI(A/B)", "AP-MS", "CRISPR", "LCV", "Response")
host_long <- host_evidence %>%
  select(gene_label, all_of(tracks)) %>%
  pivot_longer(all_of(tracks), names_to = "track", values_to = "val") %>%
  mutate(val = as.integer(val))

pB <- ggplot(host_long, aes(track, gene_label, fill = factor(val))) +
  geom_tile(colour = "white", linewidth = 0.5) +
  scale_fill_manual(values = c("0" = "#F5F5F5", "1" = "#2E8B57"), name = "Evidence") +
  labs(x = NULL, y = NULL, title = "B. Host trafficking evidence") +
  theme_pub(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    legend.position = "top"
  )

fig <- ((pA | pB) + plot_layout(widths = c(1.9, 1.0))) & panel_labels_topleft()
save_pub_pdf(fig, "Figure7_RAB1A_architecture", w = 12, h = 6.5)
