% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ridgeSpike.R
\name{ridgeSpike}
\alias{ridgeSpike}
\title{Ridge Regression Stochastic Search Variable Selection (Bernoulli-Normal Mixture)}
\usage{
ridgeSpike(formula, data, family = "gaussian", log_lik = FALSE,
  iter = 10000, warmup = 1000, adapt = 2000, chains = 4,
  thin = 1, method = "parallel", cl = makeCluster(2), ...)
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

\item{method}{Defaults to "parallel". For an alternative parallel option, choose "rjparallel" or. Otherwise, "rjags" (single core run).}

\item{cl}{Use parallel::makeCluster(# clusters) to specify clusters for the parallel methods. Defaults to two cores.}

\item{...}{Other arguments to run.jags.}
}
\value{
a runjags object
}
\description{
The Bayesian implementation of ridge regression combined with
the Bernoulli-Normal mixture model for stochastic search variable selection. 
Plug-in pseudovariances are used for the binomial and poisson likelihood functions. 

\cr
In a way this is comparable to the elastic net. The elastic net is a convex combination of L1 and L2
norm penalities, while this model is a combination of L2 and L0 penalities (albeit not convex, since
the L0 norm is not convex).

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
ridgeSpike()
}
