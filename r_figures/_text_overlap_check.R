# Detect text-overlap in the 7 figure PDFs via pdfplumber (Python) - token-efficient "visual" QA
ROOT <- Sys.getenv("LEGIONELLA_ROOT")
if (!nzchar(ROOT)) stop("Set LEGIONELLA_ROOT before running this script.")
if (!dir.exists(ROOT)) stop("LEGIONELLA_ROOT does not exist: ", ROOT)
setwd(file.path(ROOT, "15_manuscript/r_figures"))
files <- c("Figure1_data_architecture", "Figure2_cross_cohort_mechanism", "Figure3_prioritization",
           "Figure4_holdout_validation", "Figure5_evidence_coupling",
           "Figure6_trafficking_atlas", "Figure7_RAB1A_architecture")
writeLines(files, "_figlist.txt")
