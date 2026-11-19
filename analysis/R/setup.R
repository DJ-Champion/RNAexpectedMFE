# Shared setup for the analysis scripts. Run them from the project root:
#   Rscript analysis/run_all.R
# or open RNAexpected.Rproj and source the numbered scripts in order.

if (!file.exists("DESCRIPTION") || !dir.exists("data-raw")) {
  stop("Run the analysis scripts from the RNAexpected project root.")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

# Package code (model definitions, fitted coefficients, expected_mfe()), so the
# pipeline uses exactly what the package exports.
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)

# Analysis-only helpers.
for (f in list.files("analysis/R", pattern = "\\.R$", full.names = TRUE)) {
  if (basename(f) != "setup.R") source(f)
}

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/diagnostics", recursive = TRUE, showWarnings = FALSE)
dir.create("results/models", recursive = TRUE, showWarnings = FALSE)
