% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bridgeDC.R
\name{bridgeDC}
\alias{bridgeDC}
\title{Bayesian Bridge Regression with unpenalized design covariates}
\usage{
bridgeDC(formula, design.formula, data, family = "gaussian",
  kappa = 1.4, log_lik = FALSE, iter = 10000, warmup = 1000,
  adapt = 2000, chains = 4, thin = 1, method = "parallel",
  cl = makeCluster(2), ...)
}
\arguments{
\item{formula}{the model formula}

\item{design.formula}{formula for the design covariates.}

\item{data}{a data frame.}

\item{family}{one of "gaussian", "binomial", or "poisson".}

\item{kappa}{the Lp norm you wish to utilize. Default is 1.4.}

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
The Bayesian Bridge model of Mallick & Yi (2018), but with the allowance for a set of covariates that are not penalized. 
For example, you may wish to include variables such as age and gender in all models so that 
the coefficients for the other variables are penalized while controlling for these. This
is a common need in research. Bridge regression allows you to utilize different Lp norms for the shape 
of the prior through the shape parameter kappa of the power exponential distribution (also known as generalized Gaussian). 
Norms of 1 and 2 give the Laplace and Gaussian distributions respectively (corresponding to the LASSO and Ridge Regression). 
Norms smaller than 1 are very difficult to estimate directly, but have very tall modes at zero and very long, cauchy like tails. 
Values greater than 2 become increasingly platykurtic, with the uniform distribution arising as it approaches infinity. #' Using kappa = 1 yields the New Bayesian LASSO, which is a re-parameterization of the Bayesian LASSO utilizing a scale mixture of
uniform distributions to obtain the Laplacian priors (Mallick & Yi, 2014).  \cr
\cr 
\cr
JAGS has no built in power exponential distribution, so the distribution is parameterized as a uniform-gamma mixture just as in Mallick & Yi (2018). 
The parameterization is given below. For generalized linear models plug-in pseudovariances are used. \cr
\cr
Model Specification:
\cr
\if{html}{\figure{bridgeDC.png}{}}
\if{latex}{\figure{bridgeDC.png}{}}
\cr
\cr
Plugin Pseudo-Variances: \cr
\if{html}{\figure{pseudovar.png}{}}
\if{latex}{\figure{pseudovar.png}{}}
}
\examples{
bridgeDC()

}
\references{
Mallick, H. & Yi, N. (2018) Bayesian bridge regression, Journal of Applied Statistics, 45:6, 988-1008, DOI: 10.1080/02664763.2017.1324565 \cr
\cr
Mallick, H., & Yi, N. (2014). A New Bayesian Lasso. Statistics and its interface, 7(4), 571–582. doi:10.4310/SII.2014.v7.n4.a12 \cr
}
