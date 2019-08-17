#' Bayesian GAMs
#'
#' @description This implements generalized additive models (GAMs) via the basis functions
#' available in the splines and splines2 packages. For example, you can use bs(), ns(),
#' cSpline(), iSpline(), mSpline(), or the base R function poly() to generate basis functions
#' for the variables intended for non-linear modeling. The basis functions are regularized
#' using a hierarchical prior where the top-level shrinkage parameter lambda is given a 
#' gamma prior, and the predictor specific shrinkage parameters are given a DuMouchel's
#' prior. Note that predictor specific means for each variable a basis function represents,
#' i.e., if variable x1 is represented by a basis function with 6 degrees of freedom then
#' all 6 components share a single shrinkage parameter, akin to the Group LASSO. This is standard 
#' procedure for GAMs.
#' 
#' @param formula the model formula
#' @param data a data frame
#' @param family one of "gaussian", "st" (Student-t with nu=3), "binomial", or "poisson".
#' @param df degrees of freedom for prior.
#' @param s The desired prior scale. Defaults to 1. Is automatically squared within the model so
#' select a number here on the standard deviation scale.
#' @param log_lik Should the log likelihood be monitored? The default is FALSE.
#' @param iter How many post-warmup samples? Defaults to 4000.
#' @param warmup How many warmup samples? Defaults to 1000.
#' @param adapt How many adaptation steps? Defaults to 1000.
#' @param chains How many chains? Defaults to 4.
#' @param thin Thinning interval. Defaults to 1.
#' @param method Defaults to "rjparallel". Otherwise, "rjags" (single core run).
#' @param cl Use parallel::makeCluster(# clusters) to specify clusters for the parallel methods. Defaults to two cores.
#' @param ... Other arguments to run.jags.
#'
#' @return A run.jags object
#' @export
#'
#' @examples
#' gamBayes()
#'
#' 
gamBayes = function(formula, data, family = "gaussian", log_lik = FALSE, iter= 4000, warmup=1000, adapt=1000, chains=4, thin=1, method = "rjparallel", cl = makeCluster(2), ...){
  
  X = model.matrix(formula, data)
  y = model.frame(formula, data)[,1]
  idx = attr(X, "assign")[-1]
  a = sapply(unique(idx), function(x) which(idx == x))
  start = sapply(a, function(x) min(x))
  end = sapply(a, function(x) max(x))
  X = as.matrix(X[,-1])
  
  if (method == "parallel"){
    message("method switching to rjparallel to enable use of DuMouchel's prior")
    method <- "rjparallel"
  }
  
  if (family == "gaussian"){
    
    jags_glm = "model{
              tau ~ dgamma(0.01, 0.01) 
              
              lambda ~ dgamma(1, 1)
              
              for (s in 1:S){
                eta[s] ~ dmouch(lambda)
              }
              
              for (p in 1:P){
                beta[p] ~ dnorm(0, eta[idx[p]])
              }

              Intercept ~ dnorm(0, 1e-10)

              for (s in 1:S){
                    fX[1:N, s] <- beta[start[s]:end[s]] %*% t(X[1:N, start[s]:end[s]])
              }

              for (i in 1:N){
                 y[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 log_lik[i] <- logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 ySim[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
              }
              
              sigma <- sqrt(1/tau)
              Deviance <- -2 * sum(log_lik[1:N])
          }"
    
    P = ncol(X)
    write_lines(jags_glm, "jags_glm.txt")
    jagsdata = list(X = X, y = y,  N = length(y), P = ncol(X), S = max(idx), idx = idx, start = start, end = end)
    monitor = c("Intercept", "beta", "sigma", "eta" , "lambda", "fX", "Deviance", "ySim" ,"log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = lmSolve(formula, data)[1], 
                                              "beta" = lmSolve(formula, data)[-1], 
                                              "tau" = 1, 
                                              "ySim" = sample(y, length(y)), 
                                              "eta" = rep(2, max(idx)), "lambda" = 2,
                                              .RNG.name= "lecuyer::RngStream", 
                                              .RNG.seed = sample(1:10000, 1)))
  }
  
  
  
  if (family == "st"){
    
    jags_glm = "model{
              
              tau ~ dmouch(1)
              
              lambda ~ dgamma(1, 1)
              
              for (s in 1:S){
                eta[s] ~ dmouch(lambda)
              }
              
              for (p in 1:P){
                beta[p] ~ dnorm(0, eta[idx[p]])
              }

              Intercept ~ dnorm(0, 1e-10)

              for (s in 1:S){
                    fX[1:N, s] <- beta[start[s]:end[s]] %*% t(X[1:N, start[s]:end[s]])
              }
              
              for (i in 1:N){
                 mu[i] <- Intercept + sum(beta[1:P] * X[i,1:P])
                 y[i] ~ dt(mu[i], tau, 3)
                 log_lik[i] <- logdensity.t(y[i], mu[i], tau, 3)
                 ySim[i] ~ dt(mu[i], tau, 3)
              }
              
              sigma <- sqrt(1/tau)
              Deviance <- -2 * sum(log_lik[1:N])
          }"
    
    P = ncol(X)
    write_lines(jags_glm, "jags_glm.txt")
    jagsdata = list(X = X, y = y,  N = length(y), P = ncol(X), S = max(idx), idx = idx, start = start, end = end)
    monitor = c("Intercept", "beta", "sigma", "eta" , "lambda", "fX", "Deviance", "ySim" ,"log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = lmSolve(formula, data)[1], 
                                              "beta" = lmSolve(formula, data)[-1], 
                                              "tau" = 1, 
                                              "ySim" = sample(y, length(y)), 
                                              "eta" = rep(2, max(idx)), "lambda" = 2,
                                              .RNG.name= "lecuyer::RngStream", 
                                              .RNG.seed = sample(1:10000, 1)))
  }
  
  
  
  if (family == "binomial"){
    
    jags_glm = "model{
              
              lambda ~ dgamma(1, 1)
              
              for (s in 1:S){
                eta[s] ~ dmouch(lambda)
              }
              
              for (p in 1:P){
                beta[p] ~ dnorm(0, eta[idx[p]])
              }

              Intercept ~ dnorm(0, 1e-10)

              for (s in 1:S){
                    fX[1:N, s] <- beta[start[s]:end[s]] %*% t(X[1:N, start[s]:end[s]])
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
    write_lines(jags_glm, "jags_glm.txt")
    jagsdata = list(X = X, y = y,  N = length(y), P = ncol(X), S = max(idx), idx = idx, start = start, end = end)
    monitor = c("Intercept", "beta", "sigma", "eta" , "lambda", "fX", "Deviance", "ySim" ,"log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = coef(glm(formula, data, family = "binomial"))[1], 
                                              "beta" = coef(glm(formula, data, family = "binomial"))[-1],
                                              "ySim" = sample(y, length(y)), 
                                              "eta" = rep(2, max(idx)), "lambda" = 2,
                                              .RNG.name= "lecuyer::RngStream", 
                                              .RNG.seed= sample(1:10000, 1)))
  }
  
  
  if (family == "poisson"){
    
    jags_glm = "model{

              lambda ~ dgamma(1, 1)

              for (s in 1:S){
                eta[s] ~ dmouch(lambda)
              }
              
              for (p in 1:P){
                beta[p] ~ dnorm(0, eta[idx[p]])
              }

              Intercept ~ dnorm(0, 1e-10)

              for (s in 1:S){
                    fX[1:N, s] <- beta[start[s]:end[s]] %*% t(X[1:N, start[s]:end[s]])
              }
              
              for (i in 1:N){
                 log(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                 y[i] ~ dpois(psi[i])
                 log_lik[i] <- logdensity.pois(y[i], psi[i])
                 ySim[i] ~ dpois(psi[i])
              }

              Deviance <- -2 * sum(log_lik[1:N])
          }"
    
    write_lines(jags_glm, "jags_glm.txt")
    P = ncol(X)
    jagsdata = list(X = X, y = y,  N = length(y), P = ncol(X), S = max(idx), idx = idx, start = start, end = end)
    monitor = c("Intercept", "beta", "sigma", "eta" , "lambda", "fX", "Deviance", "ySim" ,"log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = coef(glm(formula, data, family = "poisson"))[1], 
                                              "ySim" = sample(y, length(y)),
                                              "eta" = rep(2, max(idx)), "lambda" = 2,
                                              "beta" = coef(glm(formula, data, family = "poisson"))[-1], 
                                              .RNG.name= "lecuyer::RngStream", 
                                              .RNG.seed= sample(1:10000, 1)))
  }
  
  
  out = run.jags(model = "jags_glm.txt", modules = c("glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, n.chains = chains, summarise = FALSE,...)
  
  if (is.null(cl) == FALSE){
    parallel::stopCluster(cl = cl)
  }
  
  file.remove("jags_glm.txt")
  return(out)
}
