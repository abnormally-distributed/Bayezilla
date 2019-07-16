#' Generalized Ridge Regression 
#'
#' @description The Bayesian implementation of the generalized ridge regression estimators discussed by
#' Ishwaran & Rao (2014) and Yuzbacsi et al. (2017). This is similar to the adaptive Bayesian LASSO
#' in that it utilizes coefficient-specific shrinkage parameters. Plug-in pseudovariances are used for 
#' the binomial and poisson likelihood functions. 
#' \cr
#' \cr
#' Model Specification: \cr
#' \cr
#' \if{html}{\figure{GRR.png}{}}
#' \if{latex}{\figure{GRR.png}{}}
#' \cr
#' \cr
#' Plugin Pseudo-Variances: \cr
#' \cr
#' \if{html}{\figure{pseudovar.png}{}}
#' \if{latex}{\figure{pseudovar.png}{}}
#'
#'
#' @param formula the model formula
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
#' 
#' Ishwaran, H. & Rao, J. (2014) Geometry and properties of generalized ridge regression in high dimensions. Contemporary Mathematics , 622. doi: 10.1090/conm/622 \cr
#' \cr
#' Yuzbacsi, B., Arashi, M., & Ahmed, S.E. (2017). Shrinkage Estimation Strategies in Generalized Ridge Regression Models Under Low/High-Dimension Regime (preprint). https://arxiv.org/abs/1707.02331v1 \cr
#' \cr
#' 
#' @return
#' a runjags object
#' @export

#' @examples
#' genRidge()
genRidge = function(formula, data, family = "gaussian", log_lik = FALSE, iter=10000, warmup=1000, adapt=2000, chains=4, thin=1, method = "parallel", cl = makeCluster(2), ...){
  
  X = model.matrix(formula, data)[,-1]
  y = model.frame(formula, data)[,1]
  
  
  if (family == "gaussian"){
    
    jags_grr = "model{
  tau ~ dgamma(.01, .01) 
  sigma2 <- 1/tau

  for (p in 1:P){
    lambda[p] ~ dgamma(0.50, 0.20)
    omega[p] <- 1 / (sigma2 / lambda[p])
    beta[p] ~ dnorm(0, omega[p])
  }
  
  Intercept ~ dnorm(0, 1e-10)
  
  for (i in 1:N){
    y[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
    log_lik[i] <- logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]), tau)
    ySim[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
  }

  sigma <- sqrt(1/tau)
  Deviance <- -2 * sum(log_lik[1:N])
}"
    
    P <- ncol(X)
    write_lines(jags_grr, "jags_grr.txt")
    jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X))
    monitor <- c("Intercept", "beta", "sigma", "lambda", "Deviance", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits <- lapply(1:chains, function(z) list("Intercept" = lmSolve(formula, data)[1], 
                                               "beta" = lmSolve(formula, data)[-1], 
                                               "lambda" = rep(1, P), 
                                               "tau" = 1, 
                                               "ySim" = y, 
                                               .RNG.name= "lecuyer::RngStream", 
                                               .RNG.seed= sample(1:10000, 1)))
    
    out = run.jags(model = "jags_grr.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
    file.remove("jags_grr.txt")
    if (!is.null(cl)) {
      parallel::stopCluster(cl = cl)
    }
    return(out)
    
  }  
  
  
  if (family == "binomial" || family == "logistic"){
    
    jags_grr = "model{
    
  for (p in 1:P){
    lambda[p] ~ dgamma(0.50, 0.20)
    omega[p] <- 1 / (sigma2 / lambda[p])
    beta[p] ~ dnorm(0, omega[p])
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
    write_lines(jags_grr, "jags_grr.txt")
    jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), sigma2 = pow(mean(y), -1) * pow(1 - mean(y), -1))
    monitor <- c("Intercept", "beta", "lambda", "Deviance", "ySim", "log_lik")
    
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    
    inits <- lapply(1:chains, function(z) list("Intercept" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "binomial", lambda = 0.025, alpha = 0, standardize = FALSE))[1,1]), 
                                               "beta" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "binomial", lambda = 0.025, alpha = 0, standardize = FALSE))[-1,1]), 
                                               "lambda" = rgamma(P, 2, 1), 
                                               "ySim" = y, 
                                               .RNG.name= "lecuyer::RngStream", 
                                               .RNG.seed= sample(1:10000, 1)))
    
    out = run.jags(model = "jags_grr.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
    file.remove("jags_grr.txt")
    if (!is.null(cl)) {
      parallel::stopCluster(cl = cl)
    }
    return(out)
  }
  
  
  if (family == "poisson"){
    
    jags_grr = "model{
    
  for (p in 1:P){
    lambda[p] ~ dgamma(0.50, 0.20)
    omega[p] <- 1 / (sigma2 / lambda[p])
    beta[p] ~ dnorm(0, omega[p])
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
  write_lines(jags_grr, "jags_grr.txt")
  jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), sigma2 = pow(mean(y) , -1))
  monitor <- c("Intercept", "beta", "lambda", "Deviance", "ySim", "log_lik")
  
  if (log_lik == FALSE){
    monitor = monitor[-(length(monitor))]
  }
  
  inits <- lapply(1:chains, function(z) list("Intercept" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "poisson", lambda = 0.025, alpha = 0, standardize = FALSE))[1,1]), 
                                             "beta" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "poisson", lambda = 0.025, alpha = 0, standardize = FALSE))[-1,1]), 
                                             "lambda" = rgamma(P, 2, 1), 
                                             "ySim" = y, 
                                             .RNG.name= "lecuyer::RngStream", 
                                             .RNG.seed= sample(1:10000, 1)))
  
  out = run.jags(model = "jags_grr.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
  file.remove("jags_grr.txt")
  if (!is.null(cl)) {
    parallel::stopCluster(cl = cl)
  }
  return(out)
  }

}

