# Figure 6. A restricted Rab/ER-Golgi core remains after influence sensitivity.
ROOT <- Sys.getenv("LEGIONELLA_ROOT")
if (!nzchar(ROOT)) stop("Set LEGIONELLA_ROOT before running this script.")
if (!dir.exists(ROOT)) stop("LEGIONELLA_ROOT does not exist: ", ROOT)
setwd(file.path(ROOT, "15_manuscript/r_figures"))
source(file.path(ROOT, "15_manuscript/r_figures/_theme.R"))
suppressMessages({library(dplyr); library(tidyr); library(readr); library(ggplot2); library(patchwork)})
SRC <- file.path(ROOT, "13_postfreeze_interpretation")

TR <- read_tsv(file.path(SRC, "a3_trafficking_gene_evidence_matrix.tsv"), show_col_types=FALSE)
FAM <- read_tsv(file.path(SRC, "a3_family_evidence_saturation.tsv"), show_col_types=FALSE)
EN <- read_tsv(file.path(SRC, "a4_panel_enrichment.tsv"), show_col_types=FALSE)
CT <- read_tsv(file.path(SRC, "a4_panel_trafficking_per_contrast.tsv"), show_col_types=FALSE)
LO <- read_tsv(file.path(SRC, "g_panel_leaveout.tsv"), show_col_types=FALSE)

rab <- TR %>% filter(family == "RAB") %>%
  mutate(n_ev=rowSums(select(., ppi_AB, apms_high, crispr_hit, crispr_clean,
                             stable_cross_cohort, lcv_spatial)))
rab_main <- rab %>% filter(n_ev >= 3 | in_frozen_panel) %>%
  arrange(desc(n_ev)) %>% mutate(gene=factor(gene, levels=rev(gene)))
tracks <- c("ppi_AB", "apms_high", "crispr_hit", "crispr_clean", "stable_cross_cohort", "lcv_spatial")
tlab <- c("PPI(A/B)", "AP-MS", "CRISPR", "Filtered\nCRISPR", "Cross-cohort\nresponse", "LCV")
rab_long <- rab_main %>% select(gene, all_of(tracks)) %>%
  pivot_longer(-gene, names_to="track", values_to="val") %>% mutate(val=as.integer(val))
pA <- ggplot(rab_long, aes(track, gene, fill=factor(val))) +
  geom_tile(colour="white", linewidth=0.4) +
  scale_fill_manual(values=c("0"="#F5F5F5", "1"="#D81B60"), name="evidence") +
  scale_x_discrete(labels=tlab) + labs(x=NULL, y=NULL, title="A. RAB family evidence atlas") +
  theme_pub(base_size=11) + theme(axis.text.x=element_text(angle=45, hjust=1, size=10), axis.text.y=element_text(size=10))

FAM_long <- FAM %>% rename(family=1) %>% pivot_longer(-family, names_to="track", values_to="frac")
FAM_long$track <- factor(FAM_long$track, levels=c("ppi_AB", "crispr_hit", "crispr_clean", "stable_cross_cohort", "lcv_spatial"))
pB <- ggplot(FAM_long, aes(track, family, fill=frac)) +
  geom_tile(colour="white", linewidth=0.4) +
  geom_text(aes(label=ifelse(frac > 0, sprintf("%.2f", frac), "")), size=3) +
  scale_fill_gradient(low="white", high="#1B4965", limits=c(0,1), name="fraction") +
  scale_x_discrete(labels=c(ppi_AB="PPI(A/B)", crispr_hit="CRISPR", crispr_clean="Filtered",
                            stable_cross_cohort="Response", lcv_spatial="LCV")) +
  labs(x=NULL, y=NULL, title="B. Family evidence saturation") +
  theme_pub(base_size=11) + theme(axis.text.x=element_text(angle=45, hjust=1, size=10), axis.text.y=element_text(size=10))

ENm <- EN %>% filter(grepl("atlas$", label)) %>%
  mutate(group = recode(label,
      "RAB_GTPase_family | panel | atlas" = "RAB GTPase\n(n=61)",
      "RAB_associated_machinery | panel | atlas" = "RAB-associated\n(n=77)",
      "trafficking_universe | panel | atlas" = "Trafficking\n(n=275)"),
    group = factor(group, levels=c("RAB GTPase\n(n=61)", "RAB-associated\n(n=77)", "Trafficking\n(n=275)")))
pC <- ggplot(ENm, aes(group, fold, fill=group)) +
  geom_col(colour="black", linewidth=0.3, width=0.6) +
  geom_text(aes(label=paste0(round(fold,1), "x")), vjust=-0.5, size=3.3) +
  scale_fill_manual(values=c("RAB GTPase\n(n=61)"="#C1121F", "RAB-associated\n(n=77)"="#E07A5F", "Trafficking\n(n=275)"="#2E86AB"), guide="none") +
  labs(x=NULL, y="Fold enrichment in 56-gene panel", title="C. Panel enrichment",
       subtitle="selection-aware null:\nadjusted P ~ 2e-4") +
  scale_y_continuous(expand=expansion(mult=c(0, 0.16))) +
  theme_pub(base_size=11) + theme(axis.text.x=element_text(size=9), axis.text.y=element_text(size=10))

CT <- CT %>% mutate(passed=factor(passed, levels=c(TRUE,FALSE), labels=c("Passing","Non-passing")))
pD <- ggplot(CT) +
  geom_segment(aes(x=0, xend=1, y=mean_abs_z_traf, yend=mean_abs_z_nontraf, colour=passed), alpha=0.5, linewidth=0.7) +
  geom_point(aes(x=0, y=mean_abs_z_traf, colour=passed), size=2.5) +
  geom_point(aes(x=1, y=mean_abs_z_nontraf, colour=passed), size=2.5) +
  scale_colour_manual(values=c("Passing"="#2E86AB", "Non-passing"="#8D99AE"), name=NULL) +
  scale_x_continuous(breaks=c(0,1), labels=c("Trafficking\nsubset", "Non-trafficking\nsubset")) +
  labs(x=NULL, y="Mean |signed z| of panel genes", title="D. Response by\ntrafficking membership",
       subtitle="12 deletion contrasts;\nmatched-null 0/12 enriched") +
  theme_pub(base_size=11) + theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10))

LO <- LO %>% mutate(
  label = recode(scenario,
    full_54_gene_panel = "Full 54-gene panel",
    minus_RAB1A = "-RAB1A",
    minus_RAB10 = "-RAB10",
    minus_RAB1A_RAB10 = "-RAB1A/-RAB10"
  ),
  label = factor(label, levels=rev(unique(label)))
)
pE <- ggplot(LO, aes(y=label)) +
  geom_segment(aes(x=null_p95, xend=observed_mean_oriented_z, yend=label), colour="#2E86AB", linewidth=0.9) +
  geom_point(aes(x=null_p95), colour="#8D99AE", size=2.3) +
  geom_point(aes(x=observed_mean_oriented_z), colour="#1B4965", size=2.8) +
  labs(x="Mean direction-oriented signed z\n(grey: null 95th percentile; blue: observed)", y=NULL,
       title="E. Panel signal persists\nafter Rab leave-out",
       subtitle="10,000 gene-label permutations;\nempirical P <= 1e-4 in every scenario") +
  theme_pub(base_size=11) +
  theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10))

fig <- (((pA | pB) / (pC | pD | pE)) +
  plot_layout(heights=c(1.15, 1))) & panel_labels_topleft()
save_pub_pdf(fig, "Figure6_trafficking_atlas", w=13, h=9.2)
