#' Expected minimum free energy from GC content and length
#'
#' Returns the MFE (kcal/mol) that RNAfold gives, on average, for a random
#' sequence of the given GC content and length. Comparing an observed MFE
#' with this value shows whether a sequence is more or less structured than
#' its composition alone would predict.
#'
#' @param gc Numeric vector of GC content as a fraction in [0, 1].
#' @param length Numeric vector of sequence lengths (nt). Must be the same
#'   length as `gc`, or either may be length 1.
#' @param model Which fitted model to use: "C" (default, power-law GC term),
#'   "A" or "B". The models are compared in docs/methods.md.
#' @return Numeric vector of expected MFE values, capped at 0.
#' @export
expected_mfe <- function(gc, length, model = "C") {
  model <- match.arg(model, names(mfe_models))
  if (!is.numeric(gc) || !is.numeric(length)) {
    stop("`gc` and `length` must be numeric.", call. = FALSE)
  }
  n_gc <- base::length(gc)
  n_len <- base::length(length)
  if (n_gc != n_len && n_gc != 1 && n_len != 1) {
    stop("`gc` and `length` must be the same length, or one of them length 1.",
         call. = FALSE)
  }
  if (any(gc < 0 | gc > 1, na.rm = TRUE)) {
    stop("`gc` must be a fraction in [0, 1] (got range ",
         paste(format(range(gc, na.rm = TRUE)), collapse = " to "),
         "). If you have percentages, divide by 100.", call. = FALSE)
  }
  if (any(length <= 0, na.rm = TRUE)) {
    stop("`length` must be positive.", call. = FALSE)
  }

  outside <- length < fitted_range$length[1] | length > fitted_range$length[2]
  n_outside <- sum(outside, na.rm = TRUE)
  if (n_outside > 0) {
    message(n_outside, " of ", n_len, " lengths are outside the fitted range (",
            fitted_range$length[1], "-", fitted_range$length[2],
            " nt); those predictions are extrapolated.")
  }

  vars <- c(as.list(mfe_coefficients[[model]]),
            list(GC_content = gc, length = length))
  raw <- eval(mfe_models[[model]]$formula[[3]], vars, baseenv())
  pmin(raw, 0)
}

#' Fitted coefficients of the expected-MFE models
#'
#' @param model A model name ("A", "B" or "C"), or `NULL` for all of them.
#' @return A named numeric vector for one model, or a named list of vectors.
#' @export
expected_mfe_coefficients <- function(model = NULL) {
  if (is.null(model)) return(mfe_coefficients)
  model <- match.arg(model, names(mfe_models))
  mfe_coefficients[[model]]
}
