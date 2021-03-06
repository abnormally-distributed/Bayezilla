% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distributions.R
\name{qal}
\alias{qal}
\title{Asymmetric Laplace quantile function}
\usage{
qal(p, mu = 0, scale = 1, k = 0.5, lower.tail = TRUE,
  log.p = FALSE)
}
\arguments{
\item{p}{vector of probabilities.}

\item{mu}{vector of location parameter values}

\item{scale}{vector of scale parameter values}

\item{k}{vector of quantile locations}

\item{lower.tail}{if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]}

\item{log.p}{if TRUE, probabilities p are given as log(p).}
}
\value{
vector
}
\description{
Asymmetric Laplace quantile function
}
\examples{
qal()
}
