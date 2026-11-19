# Candidate expected-MFE models.
#
# This file is the single source of truth for the model definitions: the
# pipeline (analysis/01_fit_models.R) fits these formulas, and expected_mfe()
# evaluates the same right-hand sides with the fitted coefficients in
# R/coefficients.R.
#
# Variables: GC_content is the GC fraction in [0, 1]; length is in nt;
# MFE is the RNAfold minimum free energy in kcal/mol.

mfe_models <- list(
  A = list(
    formula     = MFE ~ (a * GC_content^2 + b) * length + c,
    start       = list(a = -0.8, b = -0.1, c = 13),
    description = "Quadratic GC slope, constant intercept"
  ),
  B = list(
    formula     = MFE ~ (a * GC_content^2 + b) * length + (c * GC_content^2 + d),
    start       = list(a = -0.8, b = -0.1, c = 10, d = 10),
    description = "Quadratic GC slope, quadratic GC intercept"
  ),
  C = list(
    formula     = MFE ~ (a * GC_content^b + c) * length + d,
    start       = list(a = -0.84, b = 2.35, c = -0.15, d = 13.5),
    description = "Power-law GC slope, constant intercept"
  )
)

# Range of the synthetic training data. Predictions outside it are
# extrapolations.
fitted_range <- list(length = c(20, 4000), gc = c(0, 1))
