% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HS.R
\name{HS}
\alias{HS}
\title{Horseshoe}
\usage{
HS(formula, data, family = "gaussian", log_lik = FALSE, iter = 4000,
  warmup = 3000, adapt = 3000, chains = 4, thin = 2,
  method = "rjparallel", cl = makeCluster(2), ...)
}
\arguments{
\item{formula}{the model formula}

\item{data}{a data frame.}

\item{family}{one of "gaussian", "binomial", or "poisson".}

\item{log_lik}{Should the log likelihood be monitored? The default is FALSE.}

\item{iter}{How many post-warmup samples? Defaults to 10000.}

\item{warmup}{How many warmup samples? Defaults to 1000.}

\item{adapt}{How many adaptation steps? Defaults to 2000.}

\item{chains}{How many chains? Defaults to 4.}

\item{thin}{Thinning interval. Defaults to 1.}

\item{method}{Defaults to "rjparallel". For an alternative parallel option, choose "parallel" or. Otherwise, "rjags" (single core run).}

\item{cl}{Use parallel::makeCluster(# clusters) to specify clusters for the parallel methods. Defaults to two cores.}

\item{...}{Other arguments to run.jags.}
}
\value{
an rjags object
}
\description{
This is the horseshoe model described by Carvalho et al. (2010). This tends to run very quickly
even for larger data sets or larger numbers of predictors and in my experience is faster and more stable (at least
on the tested data sets!) than the same model implemetned in Stan. \cr
\cr
\cr
Model Specification: \cr 
\cr

\if{html}{\figure{Horseshoe.png}{}}
\if{latex}{\figure{Horseshoe.png}{}}
\cr
\cr
Plugin Pseudo-Variances: \cr
\cr
\if{html}{\figure{pseudovar.png}{}}
\if{latex}{\figure{pseudovar.png}{}}
\cr
}
\examples{
HS()

}
\references{
Carvalho, C. M., Polson, N. G., and Scott, J. G. (2010). The horseshoe estimator for sparse signals. Biometrika, 97(2):465–480.
}
