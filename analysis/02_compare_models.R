# Compare models A, B and C on identical data.
#
# Outputs:
#   results/tables/model_comparison.csv   overall fit, information criteria,
#                                         held-out error (replicates 9-10)
#   results/tables/rmse_by_length.csv     error by length band vs noise floor
#   results/figures/rmse_by_length.{pdf,png}
#   results/figures/diagnostics/model_{A,B,C}.{pdf,png}

source("analysis/R/setup.R")

data <- load_mfe_data()
fits <- readRDS("results/models/fits.rds")
test <- filter(data, iteration > 8)

comparison <- imap_dfr(fits$full, function(fit, name) {
  cells <- cell_means(data, model = name)
  fit_metrics(fit, data) |>
    mutate(
      model          = name,
      description    = mfe_models[[name]]$description,
      formula        = deparse1(mfe_models[[name]]$formula),
      cell_mean_rmse = sqrt(mean(cells$residual^2)),
      heldout_rmse   = sqrt(mean((test$MFE - predict_capped(fits$train[[name]], test))^2)),
      .before = 1
    )
}) |>
  mutate(delta_aic = aic - min(aic), delta_bic = bic - min(bic))

by_length <- rmse_by_length(fits$full, data)

write_csv(comparison, "results/tables/model_comparison.csv")
write_csv(by_length, "results/tables/rmse_by_length.csv")

save_figure(plot_rmse_by_length(by_length), "rmse_by_length", width = 9, height = 4.5)

# Diagnostics share one residual scale so the three models compare directly.
cells <- map(set_names(names(mfe_models)), ~ cell_means(data, model = .x))
limit <- max(map_dbl(cells, ~ quantile(abs(.x$residual), 0.995)))
for (name in names(cells)) {
  save_figure(plot_diagnostics(cells[[name]], name, c(-limit, limit)),
              paste0("model_", name), width = 12, height = 8,
              dir = "results/figures/diagnostics")
}

print(select(comparison, model, n_params, r2, rmse, weighted_rmse,
             cell_mean_rmse, heldout_rmse, delta_aic), width = Inf)
print(by_length, n = Inf)
