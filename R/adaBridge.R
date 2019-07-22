#' Adaptive Bayesian Bridge Regression
#'
#' @description The Bayesian Bridge model of Mallick & Yi (2018) combined 
#' with the Adaptive Bayesian LASSO of Leng, Tran and David Nott (2018). 
#' Bridge regression allows you to utilize different Lp norms for the shape 
#' of the prior through the shape parameter kappa of the power exponential distribution 
#' (also known as generalized Gaussian). Norms of 1 and 2 give the Laplace and Gaussian 
#' distributions respectively  (corresponding to the LASSO and Ridge Regression). 
#' Norms smaller than 1 are very difficult to estimate directly, but have 
#' very tall modes at zero and very long, cauchy like tails. 
#' Values greater than 2 become increasingly platykurtic, with the uniform 
#' distribution arising as it approaches infinity. \cr 
#' \cr
#' The individual lambdas on each parameter defined by a gamma(sh, ra) distribution,
#' where sh and ra are shape and rate hyperparameters.  Here sh and ra are given
#' gamma(4, 8) and gamma(1, 5) priors respectively. This places the expected
#' values for the shape and rate parameters at 0.50 and 0.20 respectively, 
#' which is consistent with the gamma(0.50, 0.20) prior on lambda
#' used for most other shrinkage models in this package \cr
#' \cr
#' JAGS has no built in power exponential distribution, so the distribution is 
#' parameterized as a uniform-gamma mixture just as in Mallick & Yi (2018). 
#' The parameterization is given below. For generalized linear models plug-in 
#' pseudovariances are used. \cr
#' \cr
#' Model Specification: \cr
#' \cr
#' \if{html}{\figure{adaBridge.png}{}}
#' \if{latex}{\figure{adaBridge.png}{}}
#' \cr
#' \cr
#' Plugin Pseudo-Variances: \cr
#' \if{html}{\figure{pseudovar.png}{}}
#' \if{latex}{\figure{pseudovar.png}{}}
#' \cr
#'
#' @param formula the model formula
#' @param data a data frame.
#' @param family one of "gaussian", "binomial", or "poisson".
#' @param kappa the Lp norm you wish to utilize. Default is 1.4.
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
#' @references 
#' Leng, C., Tran, M.N., & Nott, D.J. (2014). Bayesian adaptive Lasso. arXiv:1009.2300 \cr
#' \cr
#' Mallick, H. & Yi, N. (2018) Bayesian bridge regression, Journal of Applied Statistics, 45:6, 988-1008, DOI: 10.1080/02664763.2017.1324565 \cr
#' 
#' @return
#' a runjags object
#' @export
#'
#' @examples
#' adaBridge()
#' 
adaBridge = function(formula, data, family = "gaussian", kappa = 1.4, log_lik = FALSE, iter=10000, warmup=1000, adapt=2000, chains=4, thin=1, method = "parallel", cl = makeCluster(2), ...){
  
  X = model.matrix(formula, data)[,-1]
  y = model.frame(formula, data)[,1]
  
  if (family == "gaussian"){
    
    jags_adaptive_bridge = "model{
  
  sh ~ dgamma(4, 8)
  ra ~ dgamma(1, 5)
  
  tau ~ dgamma(.01, .01) 
  sigma <- sqrt(1/tau)

  for (i in 1:P){
    lambda[i] ~ dgamma(sh , ra)
    u[i] ~ dgamma( (1/kappa) + 1  , lambda[i] )
    beta[i] ~ dunif(-1 * pow(sigma * u[i], 1/kappa), pow(sigma * u[i], 1/kappa))
  }
  
  Intercept ~ dnorm(0, 1e-10)
  
  for (i in 1:N){
    y[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
    log_lik[i] <- logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]), tau)
    ySim[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
  }
  
  Deviance <- -2 * sum(log_lik[1:N])
}"
    
    P <- ncol(X)
    write_lines(jags_adaptive_bridge, "jags_adaptive_bridge.txt")
    jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), kappa = kappa)
    monitor <- c("Intercept", "beta", "sigma", "lambda", "Deviance", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits <- lapply(1:chains, function(z) list("Intercept" = lmSolve(formula, data)[1], 
                                               "beta" = rep(0, P), 
                                               "u" = rgamma(P, (1 / kappa) + 1, 1), 
                                               "lambda" = rep(1, P), 
                                               "tau" = 1, 
                                               "ySim" = sample(y, length(y)),, 
                                               .RNG.name= "lecuyer::RngStream", 
                                               .RNG.seed= sample(1:10000, 1)))
    
    out = run.jags(model = "jags_adaptive_bridge.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
    file.remove("jags_adaptive_bridge.txt")
    if (!is.null(cl)) {
      parallel::stopCluster(cl = cl)
    }
    return(out)
  }
  
  if (family == "binomial" || family == "logistic"){
    
    jags_adaptive_bridge = "model{

  sh ~ dgamma(4, 8)
  ra ~ dgamma(1, 5)

  for (i in 1:P){
    lambda[i] ~ dgamma(sh , ra)
    u[i] ~ dgamma( (1/kappa) + 1  , lambda[i] )
    beta[i] ~ dunif(-1 * pow(sigma * u[i], 1/kappa), pow(sigma * u[i], 1/kappa))
  }
  
  Intercept ~ dnorm(0, 1e-10)
  
    for (i in 1:N){
      logit(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
      y[i] ~ dbern(psi[i])
      log_lik[i] <- logdensity.bern(y[i], psi[i])
      ySim[i] ~ dbern(psi[i])
    }
  
  Deviance <- -2 * sum(log_lik[1:N])
}"
    
    P <- ncol(X)
    write_lines(jags_adaptive_bridge, "jags_adaptive_bridge.txt")
    jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), kappa = kappa, sigma = sqrt(pow(mean(y), -1) * pow(1 - mean(y), -1)))
    monitor <- c("Intercept", "beta", "lambda", "Deviance", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits <- lapply(1:chains, function(z) list("Intercept" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "binomial", lambda = 0.025, alpha = 0, standardize = FALSE))[1,1]), 
                                               "beta" = rep(0, P), 
                                               "u" = rgamma(P, (1 / kappa) + 1, 1), 
                                               "lambda" = rep(1, P), 
                                               "ySim" = sample(y, length(y)),, 
                                               .RNG.name= "lecuyer::RngStream", 
                                               .RNG.seed= sample(1:10000, 1)))
    
    out = run.jags(model = "jags_adaptive_bridge.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
    file.remove("jags_adaptive_bridge.txt")
    if (!is.null(cl)) {
      parallel::stopCluster(cl = cl)
    }
    return(out)
  }
  
  
  if (family == "poisson"){
    
    jags_adaptive_bridge = "model{

  sh ~ dgamma(4, 8)
  ra ~ dgamma(1, 5)

  for (i in 1:P){
    lambda[i] ~ dgamma(sh , ra)
    u[i] ~ dgamma( (1/kappa) + 1  , lambda[i] )
    beta[i] ~ dunif(-1 * pow(sigma * u[i], 1/kappa), pow(sigma * u[i], 1/kappa))
  }
  
  Intercept ~ dnorm(0, 1e-10)
  
  for (i in 1:N){
    log(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
    y[i] ~ dpois(psi[i])
    log_lik[i] <- logdensity.pois(y[i], psi[i])
    ySim[i] ~ dpois(psi[i])
}
              
  Deviance <- -2 * sum(log_lik[1:N])
}"
  
  P <- ncol(X)
  write_lines(jags_adaptive_bridge, "jags_adaptive_bridge.txt")
  jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), kappa = kappa, sigma = sqrt(pow(mean(y) , -1)))
  monitor <- c("Intercept", "beta", "lambda", "Deviance", "ySim", "log_lik")
  
  if (log_lik == FALSE){
    monitor = monitor[-(length(monitor))]
  }
  
  inits <- lapply(1:chains, function(z) list("Intercept" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "poisson", lambda = 0.025, alpha = 0, standardize = FALSE))[1,1]), 
                                             "beta" = rep(0, P), 
                                             "u" = rgamma(P, (1 / kappa) + 1, 1), 
                                             "lambda" = rep(1, P), 
                                             "ySim" = sample(y, length(y)),, 
                                             .RNG.name= "lecuyer::RngStream", 
                                             .RNG.seed= sample(1:10000, 1)))
  
  out = run.jags(model = "jags_adaptive_bridge.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
  file.remove("jags_adaptive_bridge.txt")
  if (!is.null(cl)) {
    parallel::stopCluster(cl = cl)
  }
  return(out)
  }

}

