# Key figures: observed mean MFE with the recommended model (C) overlaid.
#
# Outputs:
#   results/figures/mfe_vs_length.{pdf,png}
#   results/figures/mfe_vs_gc.{pdf,png}

source("analysis/R/setup.R")

data <- load_mfe_data()
cells <- cell_means(data, model = "C")

save_figure(plot_mfe_vs_length(cells), "mfe_vs_length")
save_figure(plot_mfe_vs_gc(cells), "mfe_vs_gc")
