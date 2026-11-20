# Key figures, using the recommended model (C) unless noted.
#
# Outputs:
#   results/figures/mfe_vs_length.{pdf,png}     observed mean MFE by GC, model overlaid
#   results/figures/mfe_vs_gc.{pdf,png}         observed mean MFE by length, model overlaid
#   results/figures/error_vs_length.{pdf,png}   relative error of A-C along length, with
#                                               human transcript region lengths
#   results/tables/rel_error_by_length.csv      data behind error_vs_length
#
# Figure legends are drafted in docs/methods.md rather than printed on the
# figures. plot_mfe_vs_length() and plot_mfe_vs_gc() also take log = TRUE for
# a log-magnitude MFE axis (diagnostic; not saved by default).

source("analysis/R/setup.R")

data <- load_mfe_data()
fits <- readRDS("results/models/fits.rds")
regions <- read_csv("data-raw/human_region_summary.csv", show_col_types = FALSE)
cells <- cell_means(data, model = "C")
rel <- rel_error_by_length(fits$full, data)

write_csv(rel, "results/tables/rel_error_by_length.csv")

save_figure(plot_mfe_vs_length(cells), "mfe_vs_length", width = 9, height = 6)
save_figure(plot_mfe_vs_gc(cells), "mfe_vs_gc", width = 9, height = 6)
save_figure(plot_error_vs_length(rel, regions), "error_vs_length", width = 10, height = 7)
