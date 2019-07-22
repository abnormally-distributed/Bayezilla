#' Bayesian Elastic Net
#'
#' @description The Bayesian elastic net described by Li and Lin (2010). For the binomial and poisson likelihoods
#' plug-in pseudo-variances are used. 
#' 
#' \cr
#' The model structure is given below: \cr
#' \cr
#' \cr
#' \if{html}{\figure{elasticNet.png}{}}
#' \if{latex}{\figure{elasticNet.png}{}}
#' \cr
#' \cr
#' Plugin Pseudo-Variances: \cr
#' \if{html}{\figure{pseudovar.png}{}}
#' \if{latex}{\figure{pseudovar.png}{}}
#' \cr
#' @param formula the model formula
#' @param data a data frame.
#' @param family one of "gaussian", "st" (Student-t with nu = 3), "binomial", or "poisson".
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
#' @references Li, Qing; Lin, Nan. The Bayesian elastic net. Bayesian Anal. 5 (2010), no. 1, 151--170. doi:10.1214/10-BA506. https://projecteuclid.org/euclid.ba/1340369796
#' 
#' @return A run.jags object
#' @export
#'
#' @examples
#' bayesEnet()
#'
bayesEnet  = function(formula, data, family = "gaussian", log_lik = FALSE, iter=10000, warmup=1000, adapt=2000, chains=4, thin=1, method = "parallel", cl = makeCluster(2), ...){

  X = model.matrix(formula, data)[,-1]
  y = model.frame(formula, data)[,1]
  
  if (family == "gaussian"){
  
    jags_elastic_net = "model{

              tau ~ dgamma(.01, .01)
              sigma <- sqrt(1/tau)
              lambdaL1 ~ dgamma(0.50, 0.2)
              lambdaL2 ~ dgamma(0.50, 0.2)

              Intercept ~ dnorm(0, 1e-10)

              for (p in 1:P){
                eta[p] ~ dgamma(.5, (8 * lambdaL2 * pow(sigma,2)) / pow(lambdaL1, 2)) T(1,)
                beta_prec[p] <- (lambdaL2/pow(sigma,2)) * (eta[p]/(eta[p]-1))
                beta[p] ~ dnorm(0, beta_prec[p])
              }

              for (i in 1:N){
                 y[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 log_lik[i] <- logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 ySim[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
              }
              Deviance <- -2 * sum(log_lik[1:N])
          }"

    P <- ncol(X)
    write_lines(jags_elastic_net, "jags_elastic_net.txt")
    jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X))
    monitor <- c("Intercept", "beta", "sigma", "lambdaL1", "lambdaL2", "Deviance",  "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
  inits <- lapply(1:chains, function(z) list("Intercept" = 0, 
                                             "beta" = lmSolve(formula, data)[-1], 
                                             "eta" = 1 + abs(jitter(rep(1, P), amount = .25)), 
                                             "lambdaL1" = 5, 
                                             "lambdaL2" = 15, 
                                             "tau" = 1, 
                                             "ySim" = sample(y, length(y)),
                                             .RNG.name= "lecuyer::RngStream", 
                                             .RNG.seed= sample(1:10000, 1)))

  out = run.jags(model = "jags_elastic_net.txt", modules = c("bugs on", "glm on", "dic off"), n.chains = chains, monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
  file.remove("jags_elastic_net.txt")
  if (!is.null(cl)) {
    parallel::stopCluster(cl = cl)
  }
  return(out)
  }
  
  
  if (family == "st"){
    
    jags_elastic_net = "model{

              tau ~ dgamma(.01, .01)
              sigma <- sqrt(1/tau)
              lambdaL1 ~ dgamma(0.50, 0.2)
              lambdaL2 ~ dgamma(0.50, 0.2)

              Intercept ~ dnorm(0, 1e-10)

              for (p in 1:P){
                eta[p] ~ dgamma(.5, (8 * lambdaL2 * pow(sigma,2)) / pow(lambdaL1, 2)) T(1,)
                beta_prec[p] <- (lambdaL2/pow(sigma,2)) * (eta[p]/(eta[p]-1))
                beta[p] ~ dnorm(0, beta_prec[p])
              }

              for (i in 1:N){
                 mu[i] <- Intercept + sum(beta[1:P] * X[i,1:P])
                 y[i] ~ dt(mu[i], tau, 3)
                 log_lik[i] <- logdensity.t(mu[i], tau, 3)
                 ySim[i] ~ dt(mu[i], tau, 3)
              }
              Deviance <- -2 * sum(log_lik[1:N])
          }"
    
    P <- ncol(X)
    write_lines(jags_elastic_net, "jags_elastic_net.txt")
    jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X))
    monitor <- c("Intercept", "beta", "sigma", "lambdaL1", "lambdaL2", "Deviance",  "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits <- lapply(1:chains, function(z) list("Intercept" = 0, 
                                               "beta" = lmSolve(formula, data)[-1], 
                                               "eta" = 1 + abs(jitter(rep(1, P), amount = .25)), 
                                               "lambdaL1" = 5, 
                                               "lambdaL2" = 15, 
                                               "tau" = 1, 
                                               "ySim" = sample(y, length(y)),
                                               .RNG.name= "lecuyer::RngStream", 
                                               .RNG.seed= sample(1:10000, 1)))
    
    out = run.jags(model = "jags_elastic_net.txt", modules = c("bugs on", "glm on", "dic off"), n.chains = chains, monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
    file.remove("jags_elastic_net.txt")
    if (!is.null(cl)) {
      parallel::stopCluster(cl = cl)
    }
    return(out)
  }
  
  
  if (family == "binomial"){
    
    jags_elastic_net = "model{

              lambdaL1 ~ dgamma(0.50, 0.2)
              lambdaL2 ~ dgamma(0.50, 0.2)

              Intercept ~ dnorm(0, 1e-10)

              for (p in 1:P){
                eta[p] ~ dgamma(.5, (8 * lambdaL2 * sigma2) / pow(lambdaL1, 2)) T(1,)
                beta_prec[p] <- (lambdaL2/sigma2) * (eta[p]/(eta[p]-1))
                beta[p] ~ dnorm(0, beta_prec[p])
              }

              for (i in 1:N){
                logit(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                y[i] ~ dbern(psi[i])
                log_lik[i] <- logdensity.bern(y[i], psi[i])
                ySim[i] ~ dbern(psi[i])
              }
    
              Deviance <- -2 * sum(log_lik[1:N])
          }"

    P <- ncol(X)
    
    write_lines(jags_elastic_net, "jags_elastic_net.txt")
    
    jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), sigma2 =pow(mean(y), -1) * pow(1 - mean(y), -1))
    
    monitor <- c("Intercept", "beta", "lambdaL1", "lambdaL2", "Deviance",  "ySim", "log_lik")
    
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    
    inits <- lapply(1:chains, function(z) list("Intercept" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "binomial", lambda = runif(1, .01, .15), alpha = runif(1, .01, .5), standardize = FALSE))[1,1]), 
                                               "beta" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "binomial", lambda = runif(1, .01, .15), alpha = runif(1, .01, .5), standardize = FALSE))[-1,1]), 
                                               "eta" = 1 + abs(jitter(rep(1, P), amount = .25)), 
                                               "lambdaL1" = 5, 
                                               "lambdaL2" = 15,
                                               "ySim" = sample(y, length(y)),
                                               .RNG.name= "lecuyer::RngStream", 
                                               .RNG.seed= sample(1:10000, 1)))
    
    out = run.jags(model = "jags_elastic_net.txt", modules = c("bugs on", "glm on", "dic off"), n.chains = chains, monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
    file.remove("jags_elastic_net.txt")
    if (!is.null(cl)) {
      parallel::stopCluster(cl = cl)
    }
    return(out)
  }
  
  if (family == "poisson"){
    
    jags_elastic_net = "model{

              lambdaL1 ~ dgamma(0.50, 0.2)
              lambdaL2 ~ dgamma(0.50, 0.2)

              Intercept ~ dnorm(0, 1e-10)

              for (p in 1:P){
                eta[p] ~ dgamma(.5, (8 * lambdaL2 * sigma2) / pow(lambdaL1, 2)) T(1,)
                beta_prec[p] <- (lambdaL2/sigma2) * (eta[p]/(eta[p]-1))
                beta[p] ~ dnorm(0, beta_prec[p])
              }

              for (i in 1:N){
                log(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                y[i] ~ dpois(psi[i])
                log_lik[i] <- logdensity.pois(y[i], psi[i])
                ySim[i] ~ dpois(psi[i])
              }
    
              Deviance <- -2 * sum(log_lik[1:N])
          }"
    
    P <- ncol(X)
    
    write_lines(jags_elastic_net, "jags_elastic_net.txt")
    
    jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), sigma2 = pow(mean(y), -1))
    
    monitor <- c("Intercept", "beta", "lambdaL1", "lambdaL2", "Deviance",  "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits <- lapply(1:chains, function(z) list("Intercept" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "poisson", lambda = runif(1, .01, .15), alpha = runif(1, .01, .5), standardize = FALSE))[1,1]), 
                                               "beta" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "poisson", lambda = runif(1, .01, .15), alpha = runif(1, .01, .5), standardize = FALSE))[-1,1]), 
                                               "eta" = 1 + abs(jitter(rep(1, P), amount = .25)), 
                                               "lambdaL1" = 5, 
                                               "lambdaL2" = 15, 
                                               "ySim" = sample(y, length(y)),
                                               .RNG.name= "lecuyer::RngStream", 
                                               .RNG.seed= sample(1:10000, 1)))
    
    out = run.jags(model = "jags_elastic_net.txt", modules = c("bugs on", "glm on", "dic off"), n.chains = chains, monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
    file.remove("jags_elastic_net.txt")
    if (!is.null(cl)) {
      parallel::stopCluster(cl = cl)
    }
    return(out)
  }
  
}
