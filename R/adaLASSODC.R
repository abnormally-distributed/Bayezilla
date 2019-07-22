#' Adaptive Bayesian Lasso with Unpenalized Design Covariates
#'
#' @description The Adaptive Bayesian LASSO of Leng, Tran and David Nott (2018).This function has the further allowance for a set of covariates that are not penalized. 
#' For example, you may wish to include variables such as age and gender so that  the coefficients for the other variables are 
#' penalized while controlling for these. This is a common need in research. Basically just the Bayesian LASSO of Park & Casella (2008) but with
#' individual lambdas on each parameter defined by a gamma(sh, ra) distribution, where sh and ra are shape and rate hyperparameters. 
#' Here sh and ra are given independent gamma(4, 8) and gamma(1, 5) priors respectively. This places the expected
#' values for the shape and rate parameters at 0.50 and 0.20 respectively, which is consistent with the gamma(0.25, 0.20) prior on lambda
#' used for most other shrinkage models with additional design covariates in this package. 
#' For the binomial and poisson likelihood functions the uniform-gamma scale mixture for the
#' variant of the Bayesian LASSO is adapted for use here. 
#' 
#' For alternatives 
#' see \code{\link[Bayezilla]{negLASSO}} (which is extremely similar) or \code{\link[Bayezilla]{extLASSO}}.
#' \cr
#' \cr 
#' Model Specification: \cr
#' \cr
#' \if{html}{\figure{adaLASSODC.png}{}}
#' \if{latex}{\figure{adaLASSODC.png}{}}
#'\cr
#' \cr
#' Plugin Pseudo-Variances: \cr
#' \if{html}{\figure{pseudovar.png}{}}
#' \if{latex}{\figure{pseudovar.png}{}}
#'
#' @param formula the model formula
#' @param design.formula formula for the design covariates.
#' @param data a data frame.
#' @param family one of "gaussian", "binomial", or "poisson".
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
#' Park, T., & Casella, G. (2008). The Bayesian Lasso. Journal of the American Statistical Association, 103(482), 681-686. Retrieved from http://www.jstor.org/stable/27640090 \cr 
#' \cr
#' Mallick, H., & Yi, N. (2014). A New Bayesian Lasso. Statistics and its interface, 7(4), 571–582. doi:10.4310/SII.2014.v7.n4.a12 \cr
#' \cr
#' Leng, C., Tran, M.N., & Nott, D.J. (2014). Bayesian adaptive Lasso. arXiv:1009.2300 \cr
#' @return
#' a runjags object
#' @export
#'
#' @examples
#' adaLASSODC()
#'
adaLASSODC = function(formula, design.formula, data, family = "gaussian", log_lik = FALSE, iter=10000, warmup=1000, adapt=2000, chains=4, thin=1, method = "parallel", cl = makeCluster(2), ...){
  
  X = as.matrix(model.matrix(formula, data)[,-1])
  y = model.frame(formula, data)[,1]
  FX <- as.matrix(model.matrix(design.formula, data)[, -1])
  
  if (family == "gaussian"){
    
  jags_blasso = "model{
    
  tau ~ dgamma(.01, .01)
  sh ~ dgamma(4, 8)
  ra ~ dgamma(1, 5)
    
  for (p in 1:P){
    lambda[p] ~ dgamma(sh , ra)
    eta[p] ~ dexp(lambda[p]^2 / 2)
    omega[p] <- 1 / ( (1 / tau) * eta[p])
    beta[p] ~ dnorm(0, omega[p])
  }
  
  for (f in 1:FP){
      design_beta[f] ~ dnorm(0, 1e-200)
  }

  Intercept ~ dnorm(0, 1e-10)
  for (i in 1:N){
    mu[i] <- Intercept + sum(beta[1:P] * X[i,1:P]) + sum(design_beta[1:FP] * FX[i,1:FP])
    y[i] ~ dnorm(mu[i], tau) 
    log_lik[i] <- logdensity.norm(y[i], mu[i], tau)
    ySim[i] ~ dnorm(mu[i], tau)
  }
  sigma <- sqrt(1/tau)
  Deviance <- -2 * sum(log_lik[1:N])
}"
  
    FP <- ncol(FX)
    P <- ncol(X)
    write_lines(jags_blasso, "jags_blasso.txt")
    jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), FP = FP, FX = FX)
    monitor <- c("Intercept", "beta", "design_beta", "sigma", "sh" , "ra", "lambda", "Deviance", "ySim", "log_lik")
    
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    
    inits <- lapply(1:chains, function(z) list("Intercept" = lmSolve(formula, data)[1], 
                                               .RNG.name= "lecuyer::RngStream", 
                                               .RNG.seed= sample(1:10000, 1), 
                                               "design_beta" = coef(lm(design.formula, data))[-1],
                                               "beta" = lmSolve(formula, data)[-1], 
                                               "eta" = rep(1, P), 
                                               "sh" = .5, 
                                               "ra" = .1, 
                                               "ySim" = sample(y, length(y)),
                                               "lambda" = sample(1:3, size = P, replace =TRUE), 
                                               "tau" = 1))
    
    out = run.jags(model = "jags_blasso.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, n.chains = chains, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
    file.remove("jags_blasso.txt")
    if (!is.null(cl)) {
      parallel::stopCluster(cl = cl)
    }
    return(out)
  }
  
  
  if (family == "binomial"){
    
    jags_bridge = "model{

  sh ~ dgamma(4, 8)
  ra ~ dgamma(1, 5)
  
  for (i in 1:P){
    lambda[i] ~ dgamma(sh , ra)
    u[i] ~ dgamma( 2  , lambda[i] )
    beta[i] ~ dunif(-1 * (sigma * u[i]), sigma * u[i])
  }
  
  Intercept ~ dnorm(0, 1e-10)
  
                  
  for (f in 1:FP){
      design_beta[f] ~ dnorm(0, 1e-200)
  }
  
    for (i in 1:N){
      logit(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P]) + sum(design_beta[1:FP] * FX[i,1:FP])
      y[i] ~ dbern(psi[i])
      log_lik[i] <- logdensity.bern(y[i], psi[i])
      ySim[i] ~ dbern(psi[i])
    }
  
  Deviance <- -2 * sum(log_lik[1:N])
  
}"
    FP <- ncol(FX)
    P <- ncol(X)
    write_lines(jags_bridge, "jags_bridge.txt")
    jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), sigma = sqrt(pow(mean(y), -1) * pow(1 - mean(y), -1)), FP = FP, FX = FX)
    monitor <- c("Intercept", "beta", "design_beta", "sh" , "ra", "lambda", "Deviance", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits <- lapply(1:chains, function(z) list("Intercept" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "binomial", lambda = 0.025, alpha = .5, standardize = FALSE))[1,1]), 
                                               "beta" = rep(0, P), 
                                               "design_beta" = as.vector(coef(glm(design.formula, data, family = "binomial")))[-1],
                                               "sh" = .5, 
                                               "ra" = .10, 
                                               "u" = rgamma(P, 2, 1), 
                                               "lambda" = rep(1, P), 
                                               "ySim" = sample(y, length(y)),
                                               .RNG.name= "lecuyer::RngStream", 
                                               .RNG.seed= sample(1:10000, 1)))
    
    out = run.jags(model = "jags_bridge.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor,  n.chains = chains,  data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
    file.remove("jags_bridge.txt")
    if (!is.null(cl)) {
      parallel::stopCluster(cl = cl)
    }
    return(out)
  }

  if (family == "poisson"){
    
    jags_bridge = "model{

  sh ~ dgamma(4, 8)
  ra ~ dgamma(1, 5)

  for (i in 1:P){
    lambda[i] ~ dgamma(sh , ra)
    u[i] ~ dgamma( 2 , lambda[i] )
    beta[i] ~ dunif(-1 * (sigma * u[i]), sigma * u[i])
  }
  
  Intercept ~ dnorm(0, 1e-10)
  
                  
  for (f in 1:FP){
      design_beta[f] ~ dnorm(0, 1e-200)
  }
  
  for (i in 1:N){
    log(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P]) + sum(design_beta[1:FP] * FX[i,1:FP])
    y[i] ~ dpois(psi[i])
    log_lik[i] <- logdensity.pois(y[i], psi[i])
    ySim[i] ~ dpois(psi[i])
}
              
  Deviance <- -2 * sum(log_lik[1:N])
}"
  FP <- ncol(FX)
  P <- ncol(X)
  write_lines(jags_bridge, "jags_bridge.txt")
  jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), sigma = sqrt(pow(mean(y) , -1)), FP = FP, FX = FX)
  monitor <- c("Intercept", "beta", "design_beta", "sh" , "ra", "lambda", "Deviance", "ySim", "log_lik")
  
  if (log_lik == FALSE){
    monitor = monitor[-(length(monitor))]
  }
  
  inits <- lapply(1:chains, function(z) list("Intercept" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "poisson", lambda = 0.025, alpha = 0, standardize = FALSE))[1,1]), 
                                             "beta" = rep(0, P), 
                                             "design_beta" = as.vector(coef(glm(design.formula, data, family = "poisson")))[-1],
                                             "sh" = .5, 
                                             "ra" = .10, 
                                             "u" = rgamma(P, 2, 1), 
                                             "lambda" = rep(1, P), 
                                             "ySim" = sample(y, length(y)),
                                             .RNG.name= "lecuyer::RngStream", 
                                             .RNG.seed= sample(1:10000, 1)))
  
  out = run.jags(model = "jags_bridge.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, n.chains = chains,  data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
  file.remove("jags_bridge.txt")
  if (!is.null(cl)) {
    parallel::stopCluster(cl = cl)
  }
  return(out)
  }

}


