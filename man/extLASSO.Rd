% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extLASSO.R
\name{extLASSO}
\alias{extLASSO}
\title{Extended Bayesian LASSO}
\usage{
extLASSO(formula, data, family = "gaussian", eta_prior = "classic",
  local_u = 4, top_u = 50, log_lik = FALSE, iter = 10000,
  warmup = 2000, adapt = 1000, chains = 4, thin = 3,
  method = "rjparallel", cl = makeCluster(2), ...)
}
\arguments{
\item{formula}{the model formula}

\item{data}{a data frame.}

\item{family}{one of "gaussian", "st" (Student-t with nu=3), "binomial", or "poisson".}

\item{eta_prior}{one of "classic (default)", or "gamma".}

\item{local_u}{This must be assigned a value. Default is 4.}

\item{top_u}{If using eta_prior = "classic" this must be assigned a value. Default is 50.}

\item{log_lik}{Should the log likelihood be monitored? The default is FALSE.}

\item{iter}{How many post-warmup samples? Defaults to 10000.}

\item{warmup}{How many warmup samples? Defaults to 2000.}

\item{adapt}{How many adaptation steps? Defaults to 1000.}

\item{chains}{How many chains? Defaults to 4.}

\item{thin}{Thinning interval. Defaults to 3.}

\item{method}{Defaults to "rjparallel". For an alternative parallel option, choose "parallel" or. Otherwise, "rjags" (single core run).}

\item{cl}{Use parallel::makeCluster(# clusters) to specify clusters for the parallel methods. Defaults to two cores.}

\item{...}{Other arguments to run.jags.}
}
\value{
A run.jags object
}
\description{
This is the extended Bayesian LASSO presented by
Crispin M. Mutshinda and Mikko J. Sillanpää (2010) which is an improvement on the Baysian LASSO
of Park & Casella (2008). Two variants are provided here. \cr
\cr
The first version of the model is the original specification in Mutshinda and Sillanpää (2010), labeled
"classic" in the options here. This requires you to choose upper limits for the uniform priors on both 
the top-level shrinkage hyperparameter as well as the local shrinkage parameters. These can be tuned through
model comparison if neccessary. 
\cr
\cr
Model Specification:
\cr
\if{html}{\figure{extLASSO.png}{}}
\if{latex}{\figure{extLASSO.png}{}}
\cr
\cr
The second version is the "gamma" prior. This places a gamma(0.50 , 0.20) prior on the
top level shrinkage hyperparameter. The individual shrinkage parameters are still given independent uniform(0, local_u) 
priors just as in the classic version. 
\cr
\cr
Model Specification:
\cr
\if{html}{\figure{extLASSO2.png}{}}
\if{latex}{\figure{extLASSO2.png}{}}
\cr
\cr
The author of the extended Bayesian Lasso (Sillanpää, personal communication) confirmed that gamma prior on the top
level shrinkage parameter 
does work well in some settings, which is why I opted to include the "gamma" variant.
\cr
\cr
One really nice feature of the extended Bayesian Lasso is that inclusion probabilities and Bayes Factors for the 
coefficients are easily obtained. The prior inclusion probability is given by 1/local_u, so for example
uniform(0, 2) priors on the shrinkage parameters indicate a 50% prior inclusion
probability. Common in Bayesian variable selection is to use a 20% probability if
dealing with a high dimensional problem, so for this choose local_u = 5. If you have genuine prior
information you can and should use this to guide your choice. If you are unsure, use model comparison
to select which value of u to choose. Inclusion indicators are given by a step function based on the 
marginal individual shrinkage parameter, delta = step(1 - eta). Inlcusion probabilities are then given as the number of
1s that appear in the vector of monte carlo samples out of the total number of iterations. This will appear as
the mean for each inclusion indicator in the summary. \cr
\cr
Bayes Factors for each cofficient can then be manually derived using the formula below (Mutshinda, & Sillanpää, 2010; 2012). 
\cr
\cr
\if{html}{\figure{extLASSO_BF.png}{}}
\if{latex}{\figure{extLASSO_BF.png}{}}
\cr
\cr
}
\examples{
extLASSO()

}
\references{
Li, Z.,and Sillanpää, M. J. (2012) Overview of LASSO-related penalized regression methods for quantitative trait mapping and genomic selection. Theoretical and Applied Genetics 125: 419-435. \cr
\cr
Mutshinda, C.M., & Sillanpää, M.J. (2010). Extended Bayesian LASSO for multiple quantitative trait loci mapping and unobserved phenotype prediction. Genetics, 186 3, 1067-75 . \cr 
\cr
Mutshinda, C. M., and M. J. Sillanpää (2012) A decision rule for quantitative trait locus detection under the extended Bayesian LASSO model. Genetics 192: 1483-1491. \cr
}
