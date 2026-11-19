# Goodness-of-fit summaries. Errors are computed on capped predictions, i.e.
# exactly what expected_mfe() returns.

length_bands <- c(0, 100, 500, 1000, 2000, 4000)

fit_metrics <- function(fit, data) {
  pred <- predict_capped(fit, data)
  res <- data$MFE - pred
  n <- length(res)
  k <- length(coef(fit))
  r2 <- 1 - sum(res^2) / sum((data$MFE - mean(data$MFE))^2)
  tibble(
    n_params = k,
    r2       = r2,
    adj_r2   = 1 - (1 - r2) * (n - 1) / (n - k),
    rmse     = sqrt(mean(res^2)),
    mae      = mean(abs(res)),
    bias     = mean(res),
    aic      = AIC(fit),
    bic      = BIC(fit)
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

rmse_by_length <- function(fits, data) {
  data <- mutate(data, band = cut(length, length_bands, dig.lab = 5))
  floor <- data |>
    group_by(band) |>
    group_modify(~ tibble(replicate_sd = replicate_sd(.x)))
  imap_dfr(fits, function(fit, name) {
    data |>
      mutate(pred = predict_capped(fit, data)) |>
      group_by(band) |>
      summarise(rmse = sqrt(mean((MFE - pred)^2)),
                cell_mean_rmse = cell_mean_rmse(pick(everything()), pred),
                .groups = "drop") |>
      mutate(model = name, .before = 1)
  }) |>
    left_join(floor, by = "band") |>
    mutate(band = as.character(band))
}
