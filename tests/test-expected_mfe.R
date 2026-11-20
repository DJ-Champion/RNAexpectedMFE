# Plain base-R tests, run by R CMD check (no testthat dependency).
library(RNAexpected)

expect_error <- function(expr, pattern) {
  msg <- tryCatch({ expr; NULL }, error = conditionMessage)
  stopifnot(!is.null(msg), grepl(pattern, msg))
}
expect_message <- function(expr, pattern) {
  msg <- NULL
  withCallingHandlers(expr, message = function(m) {
    msg <<- conditionMessage(m)
    invokeRestart("muffleMessage")
  })
  stopifnot(!is.null(msg), grepl(pattern, msg))
}

# Model C reproduces its published coefficients (1/length-weighted fit). A
# refit that drifts from them should fail here: downstream projects
# (mRNA_stability_R, RNAstab) hard-code these values and must be updated too.
reference_c <- function(gc, length) {
  pmin((-0.8392021 * gc^2.3535788 - 0.1517445) * length + 9.3178327, 0)
}
gc <- c(0, 0.25, 0.5, 0.75, 1)
len <- c(20, 100, 1000, 2500, 4000)
stopifnot(all(abs(expected_mfe(gc, len) - reference_c(gc, len)) < 0.01))

# All three models return capped, finite values across the fitted range.
for (m in c("A", "B", "C")) {
  grid <- expand.grid(gc = seq(0, 1, 0.1), len = seq(20, 4000, 220))
  out <- expected_mfe(grid$gc, grid$len, model = m)
  stopifnot(length(out) == nrow(grid), all(is.finite(out)), all(out <= 0))
}

# Short sequences hit the 0 kcal/mol cap rather than going positive.
stopifnot(expected_mfe(0.1, 20) == 0)

# Recycling and NA handling.
stopifnot(length(expected_mfe(0.5, c(100, 200, 300))) == 3)
stopifnot(is.na(expected_mfe(NA_real_, 100)), is.na(expected_mfe(0.5, NA_real_)))
expect_error(expected_mfe(c(0.1, 0.2), c(100, 200, 300)), "same length")

# Percent instead of fraction is rejected; bad lengths and models too.
expect_error(expected_mfe(50, 1000), "fraction")
expect_error(expected_mfe(-0.1, 1000), "fraction")
expect_error(expected_mfe(0.5, 0), "positive")
expect_error(expected_mfe(0.5, 1000, model = "D"), "should be one of")

# Out-of-range lengths only produce a message.
expect_message(expected_mfe(0.5, c(10, 1000, 5000)), "2 of 3 lengths")

# Coefficient accessor.
stopifnot(identical(names(expected_mfe_coefficients()), c("A", "B", "C")),
          identical(names(expected_mfe_coefficients("B")), c("a", "b", "c", "d")))
