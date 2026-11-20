# Key figures. The MFE-vs-length and MFE-vs-GC figures are drawn once per model
# in `figure_models` (default: the recommended model C), with that model's
# fitted coefficients in the subtitle.
#
# Outputs:
#   results/figures/mfe_vs_length_model_<M>.{pdf,png}  observed mean MFE by GC, model M overlaid
#   results/figures/mfe_vs_gc_model_<M>.{pdf,png}      observed mean MFE by length, model M overlaid
#   results/figures/error_vs_length.{pdf,png}   relative error of A-C along length, with
#                                               human transcript region lengths
#   results/tables/rel_error_by_length.csv      data behind error_vs_length
#
# Figure legends are drafted in docs/methods.md rather than printed on the
# figures. plot_mfe_vs_length() and plot_mfe_vs_gc() also take log = TRUE for
# a log-magnitude MFE axis (diagnostic; not saved by default).

source("analysis/R/setup.R")

# Any of names(mfe_models), e.g. c("A", "B", "C") for one version each.
figure_models <- "C"

data <- load_mfe_data()
fits <- readRDS("results/models/fits.rds")
regions <- read_csv("data-raw/human_region_summary.csv", show_col_types = FALSE)
cells <- cell_means(data)
rel <- rel_error_by_length(fits$full, data)

write_csv(rel, "results/tables/rel_error_by_length.csv")

for (m in figure_models) {
  save_figure(plot_mfe_vs_length(cells, model = m),
              paste0("mfe_vs_length_model_", m), width = 9, height = 6.5)
  save_figure(plot_mfe_vs_gc(cells, model = m),
              paste0("mfe_vs_gc_model_", m), width = 9, height = 6.5)
}
save_figure(plot_error_vs_length(rel, regions), "error_vs_length", width = 10, height = 7)
