# Goodness-of-fit summaries. Errors are computed on capped predictions, i.e.
# exactly what expected_mfe() returns.
#
# Unweighted RMSE is reported but favours unweighted fits by construction, so
# models are compared mainly on weighted RMSE and on relative error within
# length bands.

length_bands <- c(0, 50, 100, 150, 250, 500, 1000, 2000, 4000)

# GC range used for relative-error summaries: where most real transcript
# regions sit, and away from the 0 kcal/mol floor at very low GC.
rel_error_gc <- c(30, 70)

fit_metrics <- function(fit, data) {
  pred <- predict_capped(fit, data)
  res <- data$MFE - pred
  w <- fit_weights(data)
  n <- length(res)
  k <- length(coef(fit))
  r2 <- 1 - sum(res^2) / sum((data$MFE - mean(data$MFE))^2)
  tibble(
    n_params      = k,
    r2            = r2,
    adj_r2        = 1 - (1 - r2) * (n - 1) / (n - k),
    rmse          = sqrt(mean(res^2)),
    weighted_rmse = sqrt(sum(w * res^2) / sum(w)),
    mae           = mean(abs(res)),
    bias          = mean(res),
    aic           = AIC(fit),
    bic           = BIC(fit)
  )
}

# Error of the model against the mean of the 10 replicates in each
# (length, target GC) design cell: the systematic part of the error, with
# replicate noise averaged out.
cell_mean_rmse <- function(data, pred) {
  data |>
    mutate(res = MFE - pred) |>
    group_by(length, target_GC) |>
    summarise(res = mean(res), .groups = "drop") |>
    with(sqrt(mean(res^2)))
}

# Pooled within-cell SD of MFE across replicates: roughly the error any
# length+GC model must make on an individual sequence.
replicate_sd <- function(data) {
  data |>
    group_by(length, target_GC) |>
    summarise(v = var(MFE), .groups = "drop") |>
    with(sqrt(mean(v)))
}

# Signed relative error (%) of each design-cell mean, for cells in the
# rel_error_gc range. Negative = model predicts too little structure.
cell_rel_error <- function(fit, data) {
  data |>
    filter(target_GC >= rel_error_gc[1], target_GC <= rel_error_gc[2]) |>
    mutate(pred = predict_capped(fit, pick(everything()))) |>
    group_by(length, target_GC) |>
    summarise(obs = mean(MFE), pred = mean(pred), .groups = "drop") |>
    filter(abs(obs) >= 1) |>
    mutate(rel_error = 100 * (obs - pred) / abs(obs))
}

rel_error_by_length <- function(fits, data) {
  imap_dfr(fits, function(fit, name) {
    cell_rel_error(fit, data) |>
      group_by(length) |>
      summarise(rel_error = median(rel_error), .groups = "drop") |>
      mutate(model = name, .before = 1)
  })
}

rmse_by_length <- function(fits, data) {
  data <- mutate(data, band = cut(length, length_bands, dig.lab = 5))
  floor <- data |>
    group_by(band) |>
    group_modify(~ tibble(replicate_sd = replicate_sd(.x)))
  imap_dfr(fits, function(fit, name) {
    rel <- cell_rel_error(fit, data) |>
      mutate(band = cut(length, length_bands, dig.lab = 5)) |>
      group_by(band) |>
      summarise(median_rel_error = median(rel_error),
                median_abs_rel_error = median(abs(rel_error)), .groups = "drop")
    data |>
      mutate(pred = predict_capped(fit, data)) |>
      group_by(band) |>
      summarise(rmse = sqrt(mean((MFE - pred)^2)),
                cell_mean_rmse = cell_mean_rmse(pick(everything()), pred),
                .groups = "drop") |>
      left_join(rel, by = "band") |>
      mutate(model = name, .before = 1)
  }) |>
    left_join(floor, by = "band") |>
    mutate(band = as.character(band))
}
