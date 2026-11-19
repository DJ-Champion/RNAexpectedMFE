\name{expected_mfe_coefficients}
\alias{expected_mfe_coefficients}
\title{Fitted coefficients of the expected-MFE models}
\usage{
expected_mfe_coefficients(model = NULL)
}
\arguments{
\item{model}{A model name (\code{"A"}, \code{"B"} or \code{"C"}), or
  \code{NULL} for all of them.}
}
\value{
A named numeric vector for one model, or a named list of vectors.
}
\description{
Returns the coefficients used by \code{\link{expected_mfe}}. See that page
for the model formulas.
}
\examples{
expected_mfe_coefficients("C")
}
