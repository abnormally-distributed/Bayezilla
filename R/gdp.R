#' Generalized double pareto shrinkage prior
#'
#' @description The generalized double pareto shrinkage prior of Armagan, Dunson, & Lee (2013). Note only the Gaussian 
#' likelihood is provided because the model requires conditioning on the error variance, which GLM-families
#' do not have. \cr
#' \cr
#' This model is parameterized extremely similarly to the Bayesian LASSO of Park & Casella (2008), Normal-Exponential-Gamma prior of
#' Griffin, J. E. and Brown, P. J. (2011), and adaptive Bayesian LASSO of Leng, Tran and David Nott (2018).
#' The key feature is that this model explicitly utilizes generalized double pareto priors through a scale mixture
#' of normals, while the Bayesian LASSO utilizes double exponential priors through a scale mixture of normals.
#' The Bayesian adaptive LASSO also utilizes double exponential just as the BLASSO, but has coefficient
#' specific shrinkage parameters. The NEG-BLASSO utilizes (as the name suggests) normal-exponential-gamma
#' priors, which behave very similarly to the GDP. Both the NEG and GDP distributions
#' have a peak at zero, just as the double exponential distribution, but have very long, student-t-like tails. \cr
#' \cr
#' In this model, the coefficient specific shrinkage parameters are given gamma distributions that with shape
#' and rate parameters (alpha and zeta, respectively) each with independent gamma(0.64, 0.64) hyperpriors, which
#' have an expected value of 1 and a standard deviation of 1.25. An expected value of 1 is justified by the fact that
#' for most circumstances alpha=zeta=1 is ideal.
#' 
#' To quote directly from Armagan, Dunson, & Lee (2013): 
#' \cr 
#' \cr
#' \emph{As α grows, the density becomes lighter tailed, more peaked and the variance becomes
#' smaller, while as ζ grows, the density becomes flatter and the variance increases. Hence if
#' we increase α, we may cause unwanted bias for large signals, though causing stronger
#' shrinkage for noise-like signals; if we increase ζ we may lose the ability to shrink noise-like
#' signals, as the density is not as pronounced around zero; and finally, if we increase α and η
#' at the same rate, the variance remains constant but the tails become lighter, converging to a
#' Laplace density in the limit. This leads to over-shrinking of coefficients that are away from
#' zero. As a typical default specification for the hyperparameters, one can take α = ζ = 1. This
#' choice leads to Cauchy-like tail behavior, which is well-known to have desirable Bayesian
#' robustness properties.}
#' \cr
#' The reason for not just fixing the values at 1 is that I have observed that this does not always result
#' in much, if any, shrinkage, and that using hyperpriors results in much better sampling (better chain mixing
#' and less autocorrelation). Furthermore, hyperpriors allow the data to speak as to which values are best. 
#' 
#' \cr
#' Model Specification:
#' \cr
#' \cr
#' \if{html}{\figure{gdp.png}{}}
#' \if{latex}{\figure{gdp.png}{}}
#' \cr
#' \cr
#' The marginal probability density function for the coefficients is of the form \cr
#' \if{html}{\figure{genparetoPDF.png}{}}
#' \if{latex}{\figure{genparetoPDF.png}{}}
#' \cr 
#' Which makes the implied prior on the coefficients \cr
#' \if{html}{\figure{gpdMarginal.png}{}}
#' \if{latex}{\figure{gpdMarginal.png}{}}
#' \cr 
#'
#' @param formula the model formula
#' @param data a data frame.
#' @param alpha default is 1
#' @param zeta default is 1
#' @param log_lik Should the log likelihood be monitored? The default is FALSE.
#' @param iter How many post-warmup samples? Defaults to 10000.
#' @param warmup How many warmup samples? Defaults to 1000.
#' @param adapt How many adaptation steps? Defaults to 2000.
#' @param chains How many chains? Defaults to 4.
#' @param thin Thinning interval. Defaults to 1.
#' @param method Defaults to "parallel". For an alternative parallel option, choose "rjparallel" or. Otherwise, "rjags" (single core run).
#' @param cl Use parallel::makeCluster(# clusters) to specify clusters for the parallel methods. Defaults to two cores.
#' @param ... Other arguments to run.jags.
#'
#' @references 	Armagan, A., Dunson, D. B., & Lee, J. (2013). Generalized Double Pareto Shrinkage. Statistica Sinica, 23(1), 119–143. \cr
#'
#' @return
#' a runjags object
#' @export
#' 
#' @seealso 
#' \code{\link[Bayezilla]{adaLASSO}}
#' \code{\link[Bayezilla]{negLASSO}} 
#' \code{\link[Bayezilla]{extLASSO}} 
#' \code{\link[Bayezilla]{HS}}
#' \code{\link[Bayezilla]{HSplus}}
#' \code{\link[Bayezilla]{HSreg}}
#'
#' @examples
#' gdp()
#' 
gdp = function(formula, data, log_lik = FALSE, iter=10000, warmup=1000, adapt=2000, chains=4, thin=1, method = "parallel", cl = makeCluster(2), ...){
  
  X = model.matrix(formula, data)[,-1]
  y = model.frame(formula, data)[,1]
  
  jags_gdp = "model{
  
  tau ~ dgamma(.01, .01) 
  sigma2 <- 1/tau
  alpha ~ dgamma(0.64 , 0.64)
  zeta ~ dgamma(0.64 , 0.64)
  for (p in 1:P){
    lambda[p] ~ dgamma(alpha , zeta)
    eta[p] ~ dexp(lambda[p]^2 / 2)
    omega[p] <- 1 / (sigma2 * eta[p])
    beta[p] ~ dnorm(0, omega[p])
  }
  
  Intercept ~ dnorm(0, 1)
  
  for (i in 1:N){
    y[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
    log_lik[i] <- logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]), tau)
    ySim[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
  }
  
  sigma <- sqrt(1/tau)
  Deviance <- -2 * sum(log_lik[1:N])
}"
  
  P <- ncol(X)
  write_lines(jags_gdp, "jags_gdp.txt")
  jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X))
  monitor <- c("Intercept", "beta", "sigma", "Deviance", "alpha", "zeta", "lambda", "eta", "ySim", "log_lik")
  if (log_lik == FALSE){
    monitor = monitor[-(length(monitor))]
  }
  inits <- lapply(1:chains, function(z) list("Intercept" = 0, "beta" = rep(0, P), "alpha" = 1, "zeta" = 1, "eta" = rep(1, P), "lambda" = rep(1, P), "tau" = 1, .RNG.name= "lecuyer::RngStream", .RNG.seed= sample(1:10000, 1)))
  
  out = run.jags(model = "jags_gdp.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
  file.remove("jags_gdp.txt")
  if (!is.null(cl)) {
    parallel::stopCluster(cl = cl)
  }
  return(out)
}
