# Fit models A, B and C (R/models.R) to the full synthetic dataset, and again
# on replicates 1-8 only for the held-out check in 02_compare_models.R.
#
# Outputs:
#   R/coefficients.R                    coefficients used by expected_mfe()
#   results/tables/model_coefficients.csv
#   results/models/fits.rds             nls objects (git-ignored)

source("analysis/R/setup.R")

data <- load_mfe_data()
message("Loaded ", nrow(data), " sequences")

fits <- fit_all_models(data)
fits_train <- fit_all_models(filter(data, iteration <= 8))

saveRDS(list(full = fits, train = fits_train), "results/models/fits.rds")
write_csv(coefficient_table(fits), "results/tables/model_coefficients.csv")
write_coefficients_r(fits, nrow(data))

print(map(fits, coef))
