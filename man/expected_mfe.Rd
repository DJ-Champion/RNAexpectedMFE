\name{expected_mfe}
\alias{expected_mfe}
\title{Expected minimum free energy from GC content and length}
\usage{
expected_mfe(gc, length, model = "C")
}
\arguments{
\item{gc}{Numeric vector of GC content as a fraction in [0, 1].}
\item{length}{Numeric vector of sequence lengths (nt). Must be the same
  length as \code{gc}, or either may be length 1.}
\item{model}{Which fitted model to use: \code{"C"} (default, power-law GC
  term), \code{"A"} or \code{"B"}.}
}
\value{
Numeric vector of expected MFE values (kcal/mol), capped at 0.
}
\description{
Returns the minimum free energy (MFE) that RNAfold gives, on average, for a
random sequence of the given GC content and length. Comparing an observed MFE
with this value shows whether a sequence is more or less structured than its
composition alone would predict.
}
\details{
The models were fitted by nonlinear least squares to 406,470 synthetic
sequences (20-4000 nt, 0-100\% GC, 10 replicates per design cell) folded with
RNAfold using default settings:
\itemize{
  \item A: \eqn{(a \cdot GC^2 + b) \cdot L + c}
  \item B: \eqn{(a \cdot GC^2 + b) \cdot L + (c \cdot GC^2 + d)}
  \item C: \eqn{(a \cdot GC^b + c) \cdot L + d}
}
Model C fits best (RMSE about 19 kcal/mol overall, rising with length).
Lengths outside 20-4000 nt give a message, because those predictions are
extrapolated. A \code{gc} outside [0, 1] is an error, since it almost always
means percentages were passed.
}
\examples{
expected_mfe(gc = 0.45, length = 1500)
expected_mfe(gc = c(0.3, 0.5, 0.7), length = 1000)
expected_mfe(0.45, 1500, model = "A")
}
