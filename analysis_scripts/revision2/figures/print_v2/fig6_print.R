# =====================================================================
# Figure 6 - print-size redraw (6.75 in wide).
# Main figure keeps panels A-D. Panel E (leave-out influence) moves to
# Supplementary Figure 9A. Panel letters A-D are unchanged, so Results text
# only needs one pointer edit (Fig. 6E -> Supplementary Fig. 9A).
# Panel B drops in-cell numbers; exact fractions remain in the source data.
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

TR <- read_tsv(file.path(SRC, "a3_trafficking_gene_evidence_matrix.tsv"), show_col_types = FALSE)
FAM <- read_tsv(file.path(SRC, "a3_family_evidence_saturation.tsv"), show_col_types = FALSE)
EN <- read_tsv(file.path(SRC, "a4_panel_enrichment.tsv"), show_col_types = FALSE)
CT <- read_tsv(file.path(SRC, "a4_panel_trafficking_per_contrast.tsv"), show_col_types = FALSE)

# ---- Panel A: RAB family evidence atlas ----
rab <- TR %>% filter(family == "RAB") %>%
  mutate(n_ev = rowSums(across(c(ppi_AB, apms_high, crispr_hit, crispr_clean,
                                 stable_cross_cohort, lcv_spatial))))
rab_main <- rab %>% filter(n_ev >= 3 | in_frozen_panel) %>%
  arrange(desc(n_ev)) %>% mutate(gene = factor(gene, levels = rev(gene)))
tracks <- c("ppi_AB", "apms_high", "crispr_hit", "crispr_clean",
            "stable_cross_cohort", "lcv_spatial")
tlab <- c("PPI(A/B)", "AP-MS", "CRISPR", "Filt.\nCRISPR",
          "Cross-cohort\nresponse", "LCV")
rab_long <- rab_main %>% select(gene, all_of(tracks)) %>%
  pivot_longer(-gene, names_to = "track", values_to = "val") %>%
  mutate(val = as.integer(val), track = factor(track, levels = tracks))

pA <- ggplot(rab_long, aes(track, gene, fill = factor(val))) +
  geom_tile(colour = "white", linewidth = 0.25) +
  scale_fill_manual(values = c("0" = "#F5F5F5", "1" = "#D81B60"), name = "Evidence") +
  scale_x_discrete(labels = tlab) +
  labs(x = NULL, y = NULL, title = "A. RAB family evidence atlas") +
  theme_print(7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.8),
        axis.text.y = element_text(size = 6.3),
        legend.position = "right")

# ---- Panel B: family evidence saturation (no in-cell numbers) ----
fam_short <- c("Arf_effectors_GEF" = "Arf GEF", "RAB_effector_GEF" = "RAB GEF",
               "GARP_EARP" = "GARP/EARP", "COG_tether" = "COG",
               "ER_contact" = "ER contact", "ER_retention" = "ER retention",
               "RAB_GAP" = "RAB GAP")
FAM_long <- FAM %>% rename(family = 1) %>%
  mutate(family = ifelse(family %in% names(fam_short), fam_short[family], family)) %>%
  pivot_longer(-family, names_to = "track", values_to = "frac") %>%
  mutate(track = factor(track, levels = c("ppi_AB", "crispr_hit", "crispr_clean",
                                          "stable_cross_cohort", "lcv_spatial")))

pB <- ggplot(FAM_long, aes(track, family, fill = frac)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  scale_fill_gradient(low = "white", high = "#1B4965", limits = c(0, 1),
                      name = "Fraction") +
  scale_x_discrete(labels = c(ppi_AB = "PPI(A/B)", crispr_hit = "CRISPR",
                              crispr_clean = "Filt.\nCRISPR",
                              stable_cross_cohort = "Response", lcv_spatial = "LCV")) +
  labs(x = NULL, y = NULL, title = "B. Family evidence saturation") +
  theme_print(7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.6),
        axis.text.y = element_text(size = 5.6),
        legend.key.height = unit(1.4, "lines"),
        legend.key.width = unit(0.9, "lines"),
        plot.margin = margin(8, 3, 2, 2))

# ---- Panel C: panel enrichment ----
ENm <- EN %>% filter(grepl("atlas$", label)) %>%
  mutate(group = recode(label,
      "RAB_GTPase_family | panel | atlas" = "RAB GTPase\n(n=61)",
      "RAB_associated_machinery | panel | atlas" = "RAB-associated\n(n=77)",
      "trafficking_universe | panel | atlas" = "Trafficking\n(n=275)"),
    group = factor(group, levels = c("RAB GTPase\n(n=61)",
                                     "RAB-associated\n(n=77)",
                                     "Trafficking\n(n=275)")))

pC <- ggplot(ENm, aes(group, fold, fill = group)) +
  geom_col(colour = "black", linewidth = 0.25, width = 0.62) +
  geom_text(aes(label = paste0(round(fold, 1), "x")), vjust = -0.5, size = 2.1) +
  scale_fill_manual(values = c("RAB GTPase\n(n=61)" = "#C1121F",
                               "RAB-associated\n(n=77)" = "#E07A5F",
                               "Trafficking\n(n=275)" = "#2E86AB"), guide = "none") +
  labs(x = NULL, y = "Fold enrichment", title = "C. Panel enrichment") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  theme_print(7) +
  theme(axis.text.x = element_text(size = 5.4))

# ---- Panel D: response by trafficking membership ----
CT <- CT %>% mutate(passed = factor(passed, levels = c(TRUE, FALSE),
                                    labels = c("Passing", "Non-passing")))
CTmean <- CT %>%
  summarise(mean_abs_z_traf = mean(mean_abs_z_traf),
            mean_abs_z_nontraf = mean(mean_abs_z_nontraf))
pD <- ggplot(CT) +
  geom_segment(aes(x = 0, xend = 1, y = mean_abs_z_traf, yend = mean_abs_z_nontraf,
                   colour = passed), alpha = 0.30, linewidth = 0.35) +
  geom_segment(data = CTmean,
               aes(x = 0, xend = 1, y = mean_abs_z_traf, yend = mean_abs_z_nontraf),
               colour = "black", linewidth = 0.9, linetype = 2,
               arrow = arrow(length = unit(0.10, "cm"), type = "closed")) +
  geom_point(aes(x = 0, y = mean_abs_z_traf, colour = passed), size = 1.6) +
  geom_point(aes(x = 1, y = mean_abs_z_nontraf, colour = passed), size = 1.6) +
  scale_colour_manual(values = c("Passing" = "#2E86AB",
                                 "Non-passing" = "#8D99AE"), name = NULL) +
  scale_x_continuous(breaks = c(0, 1),
                     labels = c("Trafficking\nsubset", "Non-trafficking\nsubset")) +
  labs(x = NULL, y = "Mean |signed z| of panel genes",
       title = "D. Response by trafficking membership") +
  theme_print(7) +
  theme(axis.text.x = element_text(size = 6.2), legend.position = "right")

fig6 <- pA / ((pB | pC) + plot_layout(widths = c(1.25, 1))) / pD +
  plot_layout(heights = c(1.45, 1.3, 1))
fig6 <- fig6 & panel_labels_topleft()
save_print(fig6, "Figure6_trafficking_atlas_print", w = 6.75, h = 8.4)
