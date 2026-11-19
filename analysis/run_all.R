# Rebuild every result from data-raw/. Run from the project root:
#   Rscript analysis/run_all.R

scripts <- c("analysis/01_fit_models.R",
             "analysis/02_compare_models.R",
             "analysis/03_figures.R")

for (s in scripts) {
  message("== ", s)
  source(s, local = new.env())
}

writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
