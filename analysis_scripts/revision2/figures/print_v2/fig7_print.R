# =====================================================================
# Figure 7 - print-size redraw (6.75 in wide).
# Main figure keeps panels A and B with the same panel letters.
# Panel B drops the 7 "background" rows; those rows are shown in
# Supplementary Figure 9B so the tiered evidence matrix stays legible.
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
    aes(x = x + 0.05, xend = 0.66, y = y, yend = 0.50),
    colour = "#D35400", linewidth = 0.55, alpha = 0.8
  ) +
  geom_point(data = effectors, aes(x, y), size = 2.6, colour = "#D35400") +
  geom_text(
    data = effectors,
    aes(x = x - 0.035, y = y, label = effector),
    size = 2.35, hjust = 1, colour = "#6B2E0D", fontface = "bold"
  ) +
  geom_label(
    data = effectors,
    aes(x = (x + 0.66) / 2, y = (y + 0.50) / 2 + 0.055, label = role_short),
    size = 1.85, colour = "#7B3F00", hjust = 0.5,
    fill = "white", label.size = 0, label.padding = unit(0.08, "lines")
  ) +
  annotate("point", x = 0.66, y = 0.50, size = 8, colour = "#C1121F") +
  annotate("text", x = 0.715, y = 0.50, label = "RAB1A\nTier 1", size = 2.4,
           colour = "#8C1B12", fontface = "bold", hjust = 0, lineheight = 0.9) +
  coord_cartesian(xlim = c(-0.10, 0.95), ylim = c(0.04, 0.96), clip = "off") +
  labs(title = "A. Direct effector manipulation of RAB1A") +
  theme_print(7) +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(),
    axis.title = element_blank(), axis.line = element_blank(),
    plot.margin = margin(2, 8, 2, 2)
  )

# ---- Panel B: tiered host trafficking evidence matrix ----
host_tiered <- host_evidence %>%
  filter(tier != "background") %>%
  mutate(
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
host_long <- host_tiered %>%
  select(gene_label, all_of(tracks)) %>%
  pivot_longer(all_of(tracks), names_to = "track", values_to = "val") %>%
  mutate(val = as.integer(val), track = factor(track, levels = tracks))

pB <- ggplot(host_long, aes(track, gene_label, fill = factor(val))) +
  geom_tile(colour = "white", linewidth = 0.25) +
  scale_fill_manual(values = c("0" = "#F5F5F5", "1" = "#2E8B57"), name = "Evidence") +
  labs(x = NULL, y = NULL, title = "B. Host trafficking evidence (tiered genes)") +
  theme_print(7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.9),
        axis.text.y = element_text(size = 5.9))

fig7 <- (pA | pB) + plot_layout(widths = c(1.35, 1.0))
fig7 <- fig7 & panel_labels_topleft()
save_print(fig7, "Figure7_RAB1A_architecture_print", w = 6.75, h = 4.8)
