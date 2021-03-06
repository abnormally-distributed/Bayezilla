#' Ridge Regression Stochastic Search Variable Selection (Bernoulli-Normal Mixture)
#'
#' @description The Bayesian implementation of ridge regression combined with
#' the Bernoulli-Normal mixture model for stochastic search variable selection. 
#' Plug-in pseudovariances are used for the binomial and poisson likelihood functions. 
#' 
#' \cr
#' In a way this is comparable to the elastic net. The elastic net is a convex combination of L1 and L2
#' norm penalities, while this model is a combination of L2 and L0 penalities (albeit not convex, since
#' the L0 norm is not convex).
#' 
#' \cr
#' \cr
#' Model Specification: \cr
#' \cr
#' \if{html}{\figure{RR.png}{}}
#' \if{latex}{\figure{RR.png}{}}
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
#' @return
#' a runjags object
#' @export

#' @examples
#' ridgeSpike()
ridgeSpike = function(formula, data, family = "gaussian", log_lik = FALSE, iter=10000, warmup=1000, adapt=2000, chains=4, thin=1, method = "parallel", cl = makeCluster(2), ...){
  
  X = model.matrix(formula, data)[,-1]
  y = model.frame(formula, data)[,1]
  
  
  if (family == "gaussian"){
    
    jags_ridgeSpike = "model{
  tau ~ dgamma(.01, .01) 
  sigma2 <- 1/tau
  phi ~ dbeta(1, 1)
  
  lambda ~ dgamma(0.50 , 0.20)

  for (p in 1:P){
    delta[p] ~ dbern(phi)
    omega[p] <- 1 / (sigma2 / lambda)
    theta[p] ~ dnorm(0, omega[p])
    beta[p] <- delta[p] * theta[p]
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
    write_lines(jags_ridgeSpike, "jags_ridgeSpike.txt")
    jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X))
    monitor <- c("Intercept", "beta", "sigma", "lambda", "phi", "Deviance","delta", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    
    inits <- lapply(1:chains, function(z) list("Intercept" = lmSolve(formula, data)[1], 
                                               "theta" = lmSolve(formula, data)[-1], 
                                               "phi" = rbeta(1, 2, 2),
                                               "delta" = rep(0, P),
                                               "lambda" = 2, 
                                               "tau" = 1, 
                                               "ySim" = sample(y, length(y)),
                                               .RNG.name= "lecuyer::RngStream", 
                                               .RNG.seed= sample(1:10000, 1)))
    
    out = run.jags(model = "jags_ridgeSpike.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
    file.remove("jags_ridgeSpike.txt")
    if (!is.null(cl)) {
      parallel::stopCluster(cl = cl)
    }
    return(out)
    
  }  
  
  if (family == "binomial"){
    
    jags_ridgeSpike = "model{
    
  lambda ~ dgamma(0.50 , 0.20)
  phi ~ dbeta(1, 1)
  
  for (p in 1:P){
    omega[p] <- 1 / (sigma2 / lambda)
    theta[p] ~ dnorm(0, omega[p])
    delta[p] ~ dbern(phi)
    beta[p] <- delta[p] * theta[p]
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
    write_lines(jags_ridgeSpike, "jags_ridgeSpike.txt")
    jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), sigma2 = pow(mean(y), -1) * pow(1 - mean(y), -1))
    monitor <- c("Intercept", "beta", "lambda", "phi", "Deviance", "delta", "ySim", "log_lik")
    
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    
    inits <- lapply(1:chains, function(z) list("Intercept" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "binomial", lambda = 0.025, alpha = 0, standardize = FALSE))[1,1]), 
                                               "theta" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "binomial", lambda = 0.025, alpha = 0, standardize = FALSE))[-1,1]), 
                                               "lambda" = 2, 
                                               "phi" = rbeta(1, 2, 2),
                                               "delta" = rep(0, P),
                                               "ySim" = sample(y, length(y)),
                                               .RNG.name= "lecuyer::RngStream", 
                                               .RNG.seed= sample(1:10000, 1)))
    
    out = run.jags(model = "jags_ridgeSpike.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
    file.remove("jags_ridgeSpike.txt")
    if (!is.null(cl)) {
      parallel::stopCluster(cl = cl)
    }
    return(out)
  }
  
  
  if (family == "poisson"){
    
    jags_ridgeSpike = "model{
    
  lambda ~ dgamma(0.50 , 0.20)
  phi ~ dbeta(1, 1)

  for (p in 1:P){
    delta[p] ~ dbern(phi)
    omega[p] <- 1 / (sigma2 / lambda)
    theta[p] ~ dnorm(0, omega[p])
    beta[p] <- delta[p] * theta[p]
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
  write_lines(jags_ridgeSpike, "jags_ridgeSpike.txt")
  jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), sigma2 = pow(mean(y) , -1))
  monitor <- c("Intercept", "beta", "lambda", "phi", "Deviance","delta", "ySim", "log_lik")
  
  if (log_lik == FALSE){
    monitor = monitor[-(length(monitor))]
  }
  
  inits <- lapply(1:chains, function(z) list("Intercept" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "poisson", lambda = 0.025, alpha = 0, standardize = FALSE))[1,1]), 
                                             "theta" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "poisson", lambda = 0.025, alpha = 0, standardize = FALSE))[-1,1]), 
                                             "lambda" = 2, 
                                             "delta" = rep(0, P),
                                             "phi" = rbeta(1, 2, 2),
                                             "ySim" = sample(y, length(y)),
                                             .RNG.name= "lecuyer::RngStream", 
                                             .RNG.seed= sample(1:10000, 1)))
  
  out = run.jags(model = "jags_ridgeSpike.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
  file.remove("jags_ridgeSpike.txt")
  if (!is.null(cl)) {
    parallel::stopCluster(cl = cl)
  }
  return(out)
  }

}

