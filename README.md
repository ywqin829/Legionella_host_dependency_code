# Legionella host-dependency atlas: reproducibility code

This repository contains the analysis and figure-generation code for the Legionella host-dependency atlas. It contains no author names, email addresses, or machine-specific absolute paths.

## Configuration

Set the project root before running scripts:

```r
Sys.setenv(LEGIONELLA_ROOT = "/path/to/Legionella_host_dependency")
```

The project root should contain the numbered data and analysis directories (`01_data_raw`, `08_evidence_atlas`, `09_holdout_validation`, `13_postfreeze_interpretation`, and `15_manuscript`).

## Figure reproduction

Run `r_figures/run_all.R` after setting `LEGIONELLA_ROOT`. The scripts write PDFs to `<LEGIONELLA_ROOT>/15_manuscript/r_figures`.

## Data and environment

Large raw data files are not included in this code-only repository. They should be obtained from the public accessions listed in the manuscript and data-availability documentation. R package requirements are listed in `R_DEPENDENCIES.txt`; Python scripts require Python 3 and the packages imported by each script.

## Scope

The joint models and post-freeze sensitivity analyses are explanatory and do not alter the frozen candidate tiers or primary endpoints.
