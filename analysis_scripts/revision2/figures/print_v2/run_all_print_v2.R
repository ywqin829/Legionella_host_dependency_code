# Run all print-size redraws in this directory.
.dir <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1]))))
  getwd()
})()
source(file.path(.dir, "_print_common.R"))
for (f in c("fig5_print.R", "fig6_print.R", "fig7_print.R", "suppfig9_print.R")) {
  message("== ", f)
  source(file.path(.dir, f))
}
message("All print figures written to ", OUT)
