% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lmSolve.R
\name{bootlmSolve}
\alias{bootlmSolve}
\title{Bootstrap minimum norm regression}
\usage{
bootlmSolve(formula, data, iter = 4000, tolerance = 0.001)
}
\arguments{
\item{formula}{formula}

\item{data}{data frame}

\item{iter}{number of iterations. defaults to 4000.}

\item{tolerance}{tolerance}
}
\value{
a vector
}
\description{
Bootstrap variant of lmSolve. Returns samples of coefficients only.
Solves with moore-penrose pseudoinversion which is very fast and tolerant
to non-positive definite model matrices.
}
\examples{
bootlmSolve()

}
