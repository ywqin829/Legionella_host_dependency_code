# Supplementary Figure 2. Effector-deletion context of E-MTAB-16655 mutants
ROOT <- Sys.getenv("LEGIONELLA_ROOT")
if (!nzchar(ROOT)) ROOT <- Sys.getenv("LEGIONELLA_ROOT")
if (!nzchar(ROOT)) stop("Set LEGIONELLA_ROOT before running this script.")
if (!dir.exists(ROOT)) stop("LEGIONELLA_ROOT does not exist: ", ROOT)
setwd(file.path(ROOT, "15_manuscript/r_figures"))
source(file.path(ROOT, "15_manuscript/r_figures/_theme.R"))
SRC <- file.path(ROOT, "13_postfreeze_interpretation")
suppressMessages({library(dplyr); library(tidyr); library(readr)})

D <- read_tsv(file.path(SRC, "c1_region_redundancy_x_holdout.tsv"))
# region -> effectors (split ;)
rows <- list()
eff_cols <- c()
for (i in 1:nrow(D)) {
  effs <- strsplit(D$effectors_deleted[i], ";")[[1]]
  effs <- effs[nzchar(effs)]
  rows[[i]] <- setNames(rep(0, 0), character(0))
  for (e in effs) {
    e <- trimws(e)
    if (!e %in% eff_cols) eff_cols <- c(eff_cols, e)
  }
}
# order columns by module clusters: RAB/isolated | ARF1/RNA | SidE/SidJ | Rab1-AMP | Rab10/SidC
order_eff <- c("cegC1", "legA7", "lem4", "lem5", "legL3",
               "ralF", "legC4", "lirA", "lirB", "lirE", "lirF", "pieE", "lem16", "lem17",
               "sdeA", "sdeB", "sdeC", "sde", "sidJ", "legK2", "legS2", "lem19",
               "sidM", "sidD", "ceg30", "legA14", "legA15",
               "sidC", "sdcA", "sdeD", "ceg32", "lem26", "lem27")
eff_cols <- order_eff[order_eff %in% eff_cols]

M <- matrix(0, nrow = nrow(D), ncol = length(eff_cols), dimnames = list(D$mutant, eff_cols))
for (i in 1:nrow(D)) {
  for (e in strsplit(D$effectors_deleted[i], ";")[[1]]) {
    e <- trimws(e)
    if (e %in% eff_cols) M[i, e] <- 1
  }
}
M <- M[c("A", "B", "C", "D", "E", "F", "H", "I"), , drop = FALSE]
status <- c(A = "PASS", B = "PASS", C = "PASS", D = "PASS", E = "FAIL", F = "FAIL",
            H = "in FGHI", I = "FAIL")

M_long <- as.data.frame(M) %>% mutate(region = rownames(M)) %>%
  pivot_longer(-region, names_to = "effector", values_to = "val")
M_long$status <- status[M_long$region]
M_long$region <- factor(M_long$region, levels = rev(rownames(M)))

p <- ggplot(M_long, aes(effector, region, fill = factor(val))) +
  geom_tile(colour = "white", linewidth = 0.5) +
  scale_fill_manual(values = c("0" = "#F5F5F5", "1" = "#2E8B57"), name = "Effector deleted") +
  labs(x = "Representative Dot/Icm effectors", y = "Deletion region",
       title = "Deletion regions vs representative Dot/Icm effectors (text colour = holdout status)",
       subtitle = "Delta-F deletes the SidE family AND SidJ; Delta-H deletes SidM+SidD; Delta-I deletes SidC/SdcA") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
        axis.text.y = element_text(size = 10, colour = ifelse(M_long$status[!duplicated(M_long$region)] == "PASS", "#2E8B57",
                                                              ifelse(M_long$status[!duplicated(M_long$region)] == "FAIL", "#C1121F", "#666666"))))

save_pub_pdf(p, "SupplementaryFigure2_effector_context", w = 11, h = 5.0)
