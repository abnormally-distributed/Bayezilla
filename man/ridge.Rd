% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RR.R
\name{ridge}
\alias{ridge}
\title{Ridge Regression}
\usage{
ridge(formula, data, family = "gaussian", lambda.prior = "dmouch",
  log_lik = FALSE, iter = 8000, warmup = 1000, adapt = 1000,
  chains = 4, thin = 1, method = "parallel", cl = makeCluster(2),
  ...)
}
\arguments{
\item{formula}{the model formula}

\item{data}{a data frame.}

\item{family}{one of "gaussian", "st" (Student-t with nu = 3), "binomial", or "poisson".}

\item{lambda.prior}{either "dmouch" (the default) or "gamma"}

\item{log_lik}{Should the log likelihood be monitored? The default is FALSE.}

\item{iter}{How many post-warmup samples? Defaults to 10000.}

\item{warmup}{How many warmup samples? Defaults to 1000.}

\item{adapt}{How many adaptation steps? Defaults to 2000.}

\item{chains}{How many chains? Defaults to 4.}

\item{thin}{Thinning interval. Defaults to 1.}

\item{method}{Defaults to "rjparallel". For an alternative parallel option, choose "parallel". Otherwise, "rjags" (single core run).}

\item{cl}{Use parallel::makeCluster(# clusters) to specify clusters for the parallel methods. Defaults to two cores.}

\item{...}{Other arguments to run.jags.}
}
\value{
a runjags object
}
\description{
The Bayesian implementation of ridge regression. Plug-in pseudovariances are used for 
the binomial and poisson likelihood functions. 
\cr
\cr
Model Specification: \cr
\cr
\if{html}{\figure{RR.png}{}}
\if{latex}{\figure{RR.png}{}}
\cr
\cr
Plugin Pseudo-Variances: \cr
\cr
\if{html}{\figure{pseudovar.png}{}}
\if{latex}{\figure{pseudovar.png}{}}
}
\examples{
ridge()
}
