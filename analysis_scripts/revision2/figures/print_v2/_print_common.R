# =====================================================================
# Common setup for print-size redraws of Figures 5-7 and Supplementary
# Figure 9. Target width = 6.75 in (ASM full-page width). The frozen
# 2026-10-04 scripts and outputs are never modified; everything is written
# to revision2_2026-11-23/figures_v2_print/.
# =====================================================================

ROOT <- Sys.getenv("LEGIONELLA_ROOT", "C:/Users/yw_qi/Legionella_host_dependency")
FROZEN_FIG <- file.path(ROOT, "FINAL_SUBMISSION_2026-10-04", "figure R script")
SRC <- file.path(FROZEN_FIG, "source_data")
OUT <- file.path(ROOT, "revision2_2026-11-23", "figures_v2_print")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

source(file.path(FROZEN_FIG, "_theme.R"))
suppressMessages({
  library(ggplot2); library(patchwork); library(dplyr)
  library(readr); library(tidyr); library(jsonlite)
})

PRINT_W <- 6.75
BASE <- 7

theme_print <- function(base = BASE) {
  theme_pub(base_size = base) +
    theme(
      plot.title = element_text(size = base + 1.2, face = "bold", hjust = 0,
                                margin = margin(b = 2)),
      axis.title = element_text(size = base + 0.3, face = "bold"),
      axis.text = element_text(size = base - 0.3, colour = "black"),
      legend.title = element_text(size = base - 0.2, face = "bold"),
      legend.text = element_text(size = base - 0.5),
      legend.key.size = unit(0.8, "lines"),
      plot.margin = margin(2, 3, 2, 2)
    )
}

save_print <- function(p, file, w = PRINT_W, h) {
  ggsave(file.path(OUT, paste0(file, ".pdf")), p, width = w, height = h,
         device = cairo_pdf)
  ggsave(file.path(OUT, paste0(file, ".png")), p, width = w, height = h,
         dpi = 300, bg = "white")
  message("wrote ", file, " (", w, " x ", h, " in)")
}
