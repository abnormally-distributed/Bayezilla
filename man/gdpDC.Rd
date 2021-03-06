% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gdpDC.R
\name{gdpDC}
\alias{gdpDC}
\title{Generalized double pareto shrinkage prior with unpenalized design covariates}
\usage{
gdpDC(formula, design.formula, data, family = "gaussian",
  log_lik = FALSE, iter = 10000, warmup = 1000, adapt = 2000,
  chains = 4, thin = 1, method = "parallel", cl = makeCluster(2),
  ...)
}
\arguments{
\item{formula}{the model formula}

\item{design.formula}{formula for the design covariates.}

\item{data}{a data frame.}

\item{family}{one of "gaussian" (the default), "binomial", or "poisson"}

\item{log_lik}{Should the log likelihood be monitored? The default is FALSE.}

\item{iter}{How many post-warmup samples? Defaults to 10000.}

\item{warmup}{How many warmup samples? Defaults to 2500.}

\item{adapt}{How many adaptation steps? Defaults to 5000.}

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
The generalized double pareto shrinkage prior of Armagan, Dunson, & Lee (2013). For the binomial
and poisson likelihoods the same pseudovariance functions used in all other models in this package are used. 
This variant grants the allowance for a set of covariates that are not penalized. 
For example, you may wish to include variables such as age and gender in all models so that 
the coefficients for the other variables are penalized while controlling for these. This
is a common need in research. \cr\cr
\cr
This model is parameterized similarly to the Bayesian LASSO of Park & Casella (2008), Normal-Exponential-Gamma prior of
Griffin, J. E. and Brown, P. J. (2011), and adaptive Bayesian LASSO of Leng, Tran and David Nott (2018).
The key feature is that this model explicitly utilizes generalized double pareto priors through a scale mixture
of normals, while the Bayesian LASSO utilizes double exponential priors through a scale mixture of normals.
The Bayesian adaptive LASSO also utilizes double exponential just as the BLASSO, but has coefficient
specific shrinkage parameters. The NEG-BLASSO utilizes (as the name suggests) normal-exponential-gamma
priors, which behave very similarly to the GDP. Both the NEG and GDP distributions
have a peak at zero, just as the double exponential distribution, but have very long, student-t-like tails. \cr
\cr
In this model, the coefficient specific shrinkage parameters are given gamma distributions that with shape
and rate parameters (alpha and zeta, respectively) each with independent gamma(4, 8) hyperpriors, 
which are very concentrated near 1. \cr
\cr
To quote directly from Armagan, Dunson, & Lee (2013): \cr
\cr 
\cr
\emph{As α grows, the density becomes lighter tailed, more peaked and the variance becomes
smaller, while as ζ grows, the density becomes flatter and the variance increases. Hence if
we increase α, we may cause unwanted bias for large signals, though causing stronger
shrinkage for noise-like signals; if we increase ζ we may lose the ability to shrink noise-like
signals, as the density is not as pronounced around zero; and finally, if we increase α and η
at the same rate, the variance remains constant but the tails become lighter, converging to a
Laplace density in the limit. This leads to over-shrinking of coefficients that are away from
zero. As a typical default specification for the hyperparameters, one can take α = ζ = 1. This
choice leads to Cauchy-like tail behavior, which is well-known to have desirable Bayesian
robustness properties.} \cr
\cr
\cr
The reason for not just fixing the values at 1 is that I have observed that this does not always result
in much, if any, shrinkage, and that using hyperpriors results in much better sampling (better chain mixing
and less autocorrelation). Furthermore, hyperpriors allow the data to speak as to which values are best. 

\cr
Model Specification:
\cr
\cr
\if{html}{\figure{gdpDC.png}{}}
\if{latex}{\figure{gdpDC.png}{}}
\cr
\cr
The marginal probability density function for the coefficients is of the form \cr
\if{html}{\figure{genparetoPDF.png}{}}
\if{latex}{\figure{genparetoPDF.png}{}}
\cr 
Which makes the implied prior on the coefficients \cr
\if{html}{\figure{gdpMarginal.png}{}}
\if{latex}{\figure{gdpMarginal.png}{}}
\cr
}
\examples{
gdpDC()

}
\references{
Armagan, A., Dunson, D. B., & Lee, J. (2013). Generalized Double Pareto Shrinkage. Statistica Sinica, 23(1), 119–143. \cr
}
