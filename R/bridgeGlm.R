#' Bridge Regression (Lq penalty) 
#'
#' @description Classical shrinkage methods such as the LASSO and ridge regression correspond to what are called Lebesgue norms, 
#' abbreviated as Lp norms, where p is a non-negative real number representing the norm. The LASSO and ridge respectively correspond
#' to L1 and L2 norms. However, Frank & Friedman (1993) consider specifying other norms within a frequentist context and called this
#' "bridge regression". Dolinar (1991) figured this out two years prior in a Bayesian paradigm in a rather obscure technical report
#' prepared for NASA. A variant of Bayesian bridge regression was also put forth in Mutshinda, C. M., and M. J. Sillanpää's (2011) \code{\link[Bayezilla]{SASP}}
#' model. 
#' 
#' \cr
#' 
#' The variant of bridge regression offered here is analagous to the \code{\link[Bayezilla]{glmBayes}} function, but with the ability to specify a norm through 
#' the kappa parameter of the power exponential distribution. Since the power exponential distribution is not implemented in JAGS it is specified here through
#' a change of variables (see model specification at end of description). Lower values of kappa result in greater probability density at zero (infinite if kappa
#' is below 1) and much heavier tails, which is useful for shrinkage. When kappa > 2, the tails of the power exponential distribution begin to contract and as kappa -> infinity, 
#' becomes a uniform distribution. A notable feature of this is that the range of possible values where prior mass exists begins to contract with the tails. 
#' A suggested use of this model would be for situations where you suspect that coefficients fall within a relatively narrow range, but do not suspect that
#' a large number of them are truly exactly zero. Hence, you may wish to utilize kappa = 4 (set as the default value here). \cr
#' \cr
#' This function is not parameterized in a way that is meant to offer efficient variable selection. You may find that the horeseshoe model offered in \code{\link[Bayezilla]{HS}} or
#' \code{\link[Bayezilla]{extLASSO}} give better shrinkage. \cr
#' \cr
#' The model: \cr
#' \cr
#' ##### top level parameters ##### \cr
#' Intercept ~ normal(0, 1) \cr
#' omega ~ gamma(0.125, 0.125) \cr
#' tau ~ gamma(.01, .01) # only for Gaussian model \cr 
#' \cr 
#' \cr
#' ##### coefficient level parameters ##### \cr
#' u_i ~ gamma(1 / kappa, 1 / omega^kappa) \cr
#' z_i <- u_i^kappa \cr
#' xi_i ~ bernoulli(.5) \cr
#' beta_i <- -1^xi_i * z_i \cr
#' \cr
#' \cr
#' !!!!! NOTE !!!!!!! \cr
#' I have an implementation of bridge regression in the Stan language in the companion package to this one,
#' BayezillaPlus. It is far more likely to produce stable results.
#' 
#' \cr
#'
#' @references 
#' Dolinar, S. (1991). Maximum-entropy probability distributions under Lp-norm constraints. The Telecommunications and Data Acquisition Report, 74-87. Retrieved from https://ntrs.nasa.gov/search.jsp?R=19910009002. \cr 
#' 
#' Frank, I., & Friedman, J. (1993). A Statistical View of Some Chemometrics Regression Tools. Technometrics, 35(2), 109-135. doi:10.2307/1269656 \cr
#' 
#' Mutshinda, C. M., and M. J. Sillanpää (2011) Bayesian shrinkage analysis of QTLs under shape-adaptive shrinkage priors, and accurate re-estimation of genetic effects. Heredity 107: 405-412. \cr
#' 
#' @param formula The model formula
#' @param data A data frame.
#' @param family One of "gaussian", "binomial", or "poisson".
#' @param kappa The norm. Defaults to 4. 
#' @param log_lik Should the log likelihood be monitored? The default is FALSE.
#' @param iter How many post-warmup samples? Defaults to 6000.
#' @param warmup How many warmup samples? Defaults to 2500.
#' @param adapt How many adaptation steps? Defaults to 5000.
#' @param chains How many chains? Defaults to 4.
#' @param thin Thinning interval. Defaults to 3.
#' @param method Defaults to "parallel". For an alternative parallel option, choose "rjparallel" or. Otherwise, "rjags" (single core run).
#' @param cl Use parallel::makeCluster(# clusters) to specify clusters for the parallel methods. Defaults to two cores.
#' @param ... Other arguments to run.jags.
#'
#' @return A run.jags object
#' @export
#'
#' @examples
#' glmBridge()
#'
glmBridge= function(formula, data, family = "gaussian", kappa = 4, log_lik = FALSE, iter= 6000, warmup = 2500, adapt = 2000, chains=4, thin = 1, method = "parallel", cl = makeCluster(3), ...){
  
  X = model.matrix(formula, data)[,-1]
  y = model.frame(formula, data)[,1]
  
  if (family == "gaussian"){
    
    jags_bridge= "model{
    
  Intercept ~ dnorm(0, 1)
  tau ~ dgamma(.1, .1)
  sigma <- sqrt(1/tau)
  lambda ~ dt(0, 1, 6) T(0, )
  
  for (j in 1:P){
    eta[j] ~ dnorm(1, 1) T(0, )
    omega[j] <- lambda*eta[j] 
    pp[j]<- 1/pow(omega[j], kappa)
    u[j] ~ dgamma(1 / kappa, pp[j])
    z[j] <- pow(u[j], kappa)
    xi[j] ~ dbern(0.5)
    beta[j]<- pow(-1, xi[j])*z[j]
  }

  for(i in 1:N){
    y[i]~dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
    ySim[i]~dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
    log_lik[i] = logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]), tau)
  }
  
  Deviance <- -2 * sum(log_lik[1:N])
}"
    
    P = ncol(X)
    monitor = c("Intercept", "beta", "sigma", "Deviance", "lambda", "eta", "ySim" ,"log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    jagsdata = list(X = X, y = y, N = length(y), P = ncol(X), "kappa" = kappa)
    inits = lapply(1:chains, function(z) list("Intercept" = 0, "ySim" = y, "xi" = sample(c(0,1), size = P, replace = TRUE), "u" = rexp(P, 1), "eta" = abs(rnorm(P)), "lambda" = 2, "tau" = 1 , .RNG.name="lecuyer::RngStream", .RNG.seed= sample(1:10000, 1)))  
    write_lines(jags_bridge, "jags_bridge.txt")
  }
  
  if (family == "binomial"){
    jags_bridge= "model{
    
  Intercept ~ dnorm(0, 1)
  
  omega_raw ~ dscaled.gamma(1, 8)
  omega <- 1 / sqrt(omega_raw)
  
  for (j in 1:P){
    pp[j]<- pow(omega, 1 / kappa) 
    u[j] ~ dgamma(1 / kappa, pp[j])
    z[j] <- pow(u[j], kappa)
    xi[j] ~ dbern(0.5)
    beta[j]<- pow(-1, xi[j])*z[j]
  }
  
  for (i in 1:N){
    logit(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
    y[i] ~ dbern(psi[i])
    log_lik[i] <- logdensity.bern(y[i], psi[i])
    ySim[i] ~ dbern(psi[i])
   }
   Deviance <- -2 * sum(log_lik[1:N])
}"
    
    P = ncol(X)
    monitor = c("Intercept", "beta", "Deviance", "omega", "ySim" ,"log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    jagsdata = list(X = X, y = y, N = length(y), P = ncol(X))
    inits = lapply(1:chains, function(z) list("Intercept" = 0, "ySim" = y, "xi" = sample(c(0,1), size = P, replace = TRUE), "u" = rexp(P, 1), "omega_raw" = 2,  .RNG.name="lecuyer::RngStream", .RNG.seed= sample(1:10000, 1)))  
    write_lines(jags_bridge, "jags_bridge.txt")
  }
  
  if (family == "poisson"){
    jags_bridge= "model{
  
  Intercept ~ dnorm(0, 1)
  omega_raw ~ dscaled.gamma(1, 8)
  omega <- 1 / sqrt(omega_raw)
  
  for (j in 1:P){
    pp[j]<- pow(omega, 1 / kappa) 
    u[j] ~ dgamma(1 / kappa, pp[j])
    z[j] <- pow(u[j], kappa)
    xi[j] ~ dbern(0.5)
    beta[j]<- pow(-1, xi[j])*z[j]
  }
  
  for (i in 1:N){
    log(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
    y[i] ~ dpois(psi[i])
    log_lik[i] <- logdensity.pois(y[i], psi[i])
    ySim[i] ~ dpois(psi[i])
   }
   Deviance <- -2 * sum(log_lik[1:N])
}"
    
    P = ncol(X)
    monitor = c("Intercept", "beta", "Deviance", "omega", "ySim" ,"log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    jagsdata = list(X = X, y = y, N = length(y), P = ncol(X))
    inits = lapply(1:chains, function(z) list("Intercept" = 0, "ySim" = y, "xi" = sample(c(0,1), size = P, replace = TRUE), "u" = rexp(P, 1), "omega_raw" = 2, .RNG.name="lecuyer::RngStream", .RNG.seed= sample(1:10000, 1)))  
    write_lines(jags_bridge, "jags_bridge.txt")
  }
  
  out = run.jags(model = "jags_bridge.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, n.chains = chains, method = method, cl = cl, summarise = FALSE, ...)
  file.remove("jags_bridge.txt")
  if (is.null(cl) == FALSE){
    parallel::stopCluster(cl = cl)
  }
  return(out)
}