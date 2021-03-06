% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HSreg.R
\name{HSreg}
\alias{HSreg}
\title{Regularized Horseshoe}
\usage{
HSreg(formula, data, family = "gaussian", phi = 0.5, slab_scale = 2,
  slab_df = 3, log_lik = FALSE, iter = 10000, warmup = 4000,
  adapt = 5000, chains = 4, thin = 1, method = "parallel",
  cl = makeCluster(2), ...)
}
\arguments{
\item{formula}{the model formula}

\item{data}{a data frame.}

\item{family}{one of "gaussian", "binomial", or "poisson".}

\item{phi}{your prior guess on the inclusion probability. Defaults to .50. Best way to come up with a figure is a prior guess on how many coefficients are non-zero out of the total number of predictors.}

\item{slab_scale}{the standard deviation of the "slab". Defaults to 2.}

\item{slab_df}{the degrees of freedom fo the slab. Higher degrees of freedom give increased L2-like regularization. Defaults to 3.}

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
This is the horseshoe model described by Piironen & Vehtari (2017). This tends to run very quickly
even for larger data sets or larger numbers of predictors and in my experience is faster and more stable (at least
on the tested data sets!) than the same model implemetned in Stan. If the horseshoe+ is analagous to the 
adaptive Bayesian LASSO, then this could be compared to the Bayesian Elastic Net in that it imposes a combination
of different shrinkage penalties (the elastic net being a combination of L1 and L2, and the regularized horseshoe being
a combination of sub-L1 and student-t penalties). \cr
\cr
Model Specification: \cr 
\cr

\if{html}{\figure{regularizedHorseshoe.png}{}}
\if{latex}{\figure{regularizedHorseshoe.png}{}}
\cr
\cr
Plugin Pseudo-Variances: \cr
\cr
\if{html}{\figure{pseudovar.png}{}}
\if{latex}{\figure{pseudovar.png}{}}
\cr
}
\examples{
HSreg()

}
\references{
Piironen, Juho; Vehtari, Aki. Sparsity information and regularization in the horseshoe and other shrinkage priors. Electron. J. Statist. 11 (2017), no. 2, 5018--5051. doi:10.1214/17-EJS1337SI. https://projecteuclid.org/euclid.ejs/1513306866 \cr
\cr
Carvalho, C. M., Polson, N. G., and Scott, J. G. (2010). The horseshoe estimator for sparse signals. Biometrika, 97(2):465–480. \cr
}
