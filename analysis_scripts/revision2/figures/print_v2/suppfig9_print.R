# =====================================================================
# Supplementary Figure 9 - print-size redraw (6.75 in wide).
# A: leave-out influence (formerly Figure 6E)
# B: background rows of the host trafficking matrix (formerly in Figure 7B)
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

LO <- read_tsv(file.path(SRC, "g_panel_leaveout.tsv"), show_col_types = FALSE)

# ---- Panel A: leave-out influence ----
lo <- LO %>%
  mutate(label = factor(recode(scenario,
      full_54_gene_panel = "Full 54-gene panel",
      minus_RAB1A = "Minus RAB1A",
      minus_RAB10 = "Minus RAB10",
      minus_RAB1A_RAB10 = "Minus RAB1A + RAB10"),
    levels = c("Full 54-gene panel", "Minus RAB1A", "Minus RAB10",
               "Minus RAB1A + RAB10")),
    p_lab = ifelse(empirical_p < 0.001, "p < 0.001", sprintf("p = %.3f", empirical_p)))

pA <- ggplot(lo, aes(observed_mean_oriented_z, label)) +
  geom_segment(aes(x = null_p95, xend = observed_mean_oriented_z, yend = label),
               colour = "#8D99AE", linewidth = 0.45, linetype = 2) +
  geom_point(size = 1.9, colour = "#1B4965") +
  geom_point(aes(x = null_p95), shape = 124, size = 3.0, colour = "#C1121F") +
  geom_text(aes(label = sprintf("%.2f", observed_mean_oriented_z)),
            hjust = -0.4, size = 2.1) +
  scale_x_continuous(limits = c(0.95, 4.25), breaks = c(1, 2, 3, 4),
                     expand = expansion(mult = c(0.02, 0.05))) +
  labs(x = "Mean oriented z of panel genes (|: null 95th percentile)", y = NULL,
       title = "A. Leave-out influence") +
  theme_print(7) +
  theme(axis.text.y = element_text(size = 6.0))

# ---- Panel B: background rows of the trafficking evidence matrix ----
host_bg <- read_tsv(file.path(SRC, "fig7_host_trafficking_evidence_finaltier.tsv"),
                    show_col_types = FALSE) %>%
  filter(tier == "background") %>%
  mutate(gene_label = factor(gene, levels = rev(gene)))
tracks <- c("PPI(A/B)", "AP-MS", "CRISPR", "LCV", "Response")
bg_long <- host_bg %>%
  select(gene_label, all_of(tracks)) %>%
  pivot_longer(all_of(tracks), names_to = "track", values_to = "val") %>%
  mutate(val = as.integer(val), track = factor(track, levels = tracks))

pB <- ggplot(bg_long, aes(track, gene_label, fill = factor(val))) +
  geom_tile(colour = "white", linewidth = 0.25) +
  scale_fill_manual(values = c("0" = "#F5F5F5", "1" = "#2E8B57"), name = "Evidence") +
  labs(x = NULL, y = NULL, title = "B. Background rows of the trafficking matrix") +
  theme_print(7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.8),
        axis.text.y = element_text(size = 5.8))

fig9 <- (pA | pB) + plot_layout(widths = c(1.25, 1.0))
fig9 <- fig9 & panel_labels_topleft()
save_print(fig9, "SupplementaryFigure9_print", w = 6.75, h = 3.2)
