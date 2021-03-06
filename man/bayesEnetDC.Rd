% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bayesEnetDC.R
\name{bayesEnetDC}
\alias{bayesEnetDC}
\title{Bayesian Elastic Net with additional unpenalized design covariates}
\usage{
bayesEnetDC(formula, design.formula, data, family = "gaussian",
  log_lik = FALSE, iter = 10000, warmup = 1000, adapt = 2000,
  chains = 4, thin = 1, method = "parallel", cl = makeCluster(2),
  ...)
}
\arguments{
\item{formula}{the model formula}

\item{design.formula}{formula for the design covariates.}

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
A run.jags object
}
\description{
The Bayesian elastic net described by Li and Lin (2010). This function has the further allowance for a set of covariates that are not penalized. 
For example, you may wish to include variables such as age and gender so that  the coefficients for the other variables are 
penalized while controlling for these. This is a common need in research.  For the binomial and poisson likelihoods
plug-in pseudo-variances are used. 


\cr
The model structure is given below: \cr
\cr
\cr
\if{html}{\figure{elasticNetDC.png}{}}
\if{latex}{\figure{elasticNetDC.png}{}}
\cr
\cr
\cr
Plugin Pseudo-Variances: \cr
\if{html}{\figure{pseudovar.png}{}}
\if{latex}{\figure{pseudovar.png}{}}
\cr
}
\examples{
bayesEnetDC()

}
\references{
Li, Qing; Lin, Nan. The Bayesian elastic net. Bayesian Anal. 5 (2010), no. 1, 151--170. doi:10.1214/10-BA506. https://projecteuclid.org/euclid.ba/1340369796
}
