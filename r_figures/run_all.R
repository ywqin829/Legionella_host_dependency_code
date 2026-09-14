# Convenience runner. Detailed plotting code remains in the individual figure scripts.
ROOT <- Sys.getenv("LEGIONELLA_ROOT")
if (!nzchar(ROOT)) stop("Set LEGIONELLA_ROOT before running this script.")
if (!dir.exists(ROOT)) stop("LEGIONELLA_ROOT does not exist: ", ROOT)
FIG_DIR <- file.path(ROOT, "15_manuscript/r_figures")
setwd(FIG_DIR)
figs <- c(
  "Figure1_data_architecture", "Figure2_cross_cohort_mechanism",
  "Figure3_prioritization", "Figure4_holdout_validation",
  "Figure5_evidence_coupling", "Figure6_trafficking_atlas",
  "Figure7_RAB1A_architecture", "SupplementaryFigure1_evidence_fragility",
  "SupplementaryFigure2_effector_context", "SupplementaryFigure3_species_composition",
  "SupplementaryFigure4_depmap_coessentiality", "SupplementaryFigure5_RAB_full_atlas", "SupplementaryFigure6_ascertainment_sensitivity"
)
for (f in figs) {
  cat("== running", f, "==\n")
  source(file.path(FIG_DIR, paste0(f, ".R")), local = FALSE)
}
cat("All", length(figs), "figure scripts completed in:", FIG_DIR, "\n")