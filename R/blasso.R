#' Bayesian Lasso 
#'
#' @description The Bayesian LASSO of Park & Casella (2008).  The Bayesian Lasso is equivalent to using independent double exponential 
#' (Laplace distribution) priors on the 
#' coefficients with a scale of sigma / lambda. However, doing this directly results in slow convergence and poor mixing.
#' The Laplace distribution can be expressed as a scale mixture of normals with an exponential distribution as the scale 
#' parameter. This is the method that Park & Casella (2008) utilize and the method that is utilized here. The hierarchical 
#' structure of the prior distribution is given below. \cr
#' \cr
#' Note that for the binomial and poisson likelihood functions 
#' the New Bayesian LASSO is used, which is a re-parameterization of the Bayesian LASSO utilizing a scale mixture of
#' uniform distributions to obtain the Laplacian priors (Mallick & Yi, 2014). I have found that this variant simply 
#' samples faster for the binomial and poisson models. Plug-in pseudovariances are used for these. 
#' \cr
#' Model Specification:
#' \cr
#' \if{html}{\figure{blasso.png}{}}
#' \if{latex}{\figure{blasso.png}{}}
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
#' @references Park, T., & Casella, G. (2008). The Bayesian Lasso. Journal of the American Statistical Association, 103(482), 681-686. Retrieved from http://www.jstor.org/stable/27640090 \cr
#' \cr
#' Mallick, H., & Yi, N. (2014). A New Bayesian Lasso. Statistics and its interface, 7(4), 571–582. doi:10.4310/SII.2014.v7.n4.a12 \cr
#' 
#' @return
#' a runjags object
#' @export

#' @examples
#' blasso()
blasso = function(formula, data, family = "gaussian",  log_lik = FALSE, iter=10000, warmup=1000, adapt=2000, chains=4, thin=1, method = "parallel", cl = makeCluster(2), ...){
    
    X = model.matrix(formula, data)[,-1]
    y = model.frame(formula, data)[,1]

    
if (family == "gaussian"){
      
  jags_blasso = "model{
  tau ~ dgamma(.01, .01) 
  sigma2 <- 1/tau
  lambda ~ dgamma(0.5 , 0.20)
  
  for (p in 1:P){
    eta[p] ~ dexp(lambda^2 / 2)
    omega[p] <- 1 / (sigma2 * eta[p])
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
      write_lines(jags_blasso, "jags_blasso.txt")
      jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X))
      monitor <- c("Intercept", "beta", "sigma", "lambda", "Deviance", "ySim", "log_lik")
      if (log_lik == FALSE){
        monitor = monitor[-(length(monitor))]
      }
      inits <- lapply(1:chains, function(z) list("Intercept" = lmSolve(formula, data)[1], 
                                                 "beta" = lmSolve(formula, data)[-1], 
                                                 "eta" = rep(1, P), 
                                                 "lambda" = 2, 
                                                 "tau" = 1, 
                                                 "ySim" = y, 
                                                 .RNG.name= "lecuyer::RngStream", 
                                                 .RNG.seed= sample(1:10000, 1)))
      
  out = run.jags(model = "jags_blasso.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
  file.remove("jags_blasso.txt")
  if (!is.null(cl)) {
      parallel::stopCluster(cl = cl)
  }
  return(out)
  
}  
  
  
  if (family == "binomial" || family == "logistic"){
    
    jags_bridge = "model{
    
  lambda ~ dgamma(0.5 , 0.20)
  
  for (i in 1:P){
    u[i] ~ dgamma( 2  , lambda )
    beta[i] ~ dunif(-1 * (sigma * u[i]), sigma * u[i])
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
    write_lines(jags_bridge, "jags_bridge.txt")
    jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), sigma = sqrt(pow(mean(y), -1) * pow(1 - mean(y), -1)))
    monitor <- c("Intercept", "beta", "lambda", "Deviance", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits <- lapply(1:chains, function(z) list("Intercept" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "binomial", lambda = 0.025, alpha = 0, standardize = FALSE))[1,1]), 
                                               "beta" = rep(0, P), 
                                               "u" = rgamma(P, 2, 1), 
                                               "lambda" = 1, 
                                               "ySim" = y, 
                                               .RNG.name= "lecuyer::RngStream", 
                                               .RNG.seed= sample(1:10000, 1)))
    
    out = run.jags(model = "jags_bridge.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
    file.remove("jags_bridge.txt")
    if (!is.null(cl)) {
      parallel::stopCluster(cl = cl)
    }
    return(out)
  }
  
  
  if (family == "poisson"){
    
    jags_bridge = "model{
    
  lambda ~ dgamma(0.5 , 0.20)
  
  for (i in 1:P){
    u[i] ~ dgamma( 2 , lambda )
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
  write_lines(jags_bridge, "jags_bridge.txt")
  jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), sigma = sqrt(pow(mean(y) , -1)))
  monitor <- c("Intercept", "beta", "lambda", "Deviance", "ySim", "log_lik")
  
  if (log_lik == FALSE){
    monitor = monitor[-(length(monitor))]
  }
  
  inits <- lapply(1:chains, function(z) list("Intercept" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "poisson", lambda = 0.025, alpha = 0, standardize = FALSE))[1,1]), 
                                             "beta" = rep(0, P), 
                                             "u" = rgamma(P, 2, 1), 
                                             "lambda" = 1, 
                                             "ySim" = y, 
                                             .RNG.name= "lecuyer::RngStream", 
                                             .RNG.seed= sample(1:10000, 1)))
  
  out = run.jags(model = "jags_bridge.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
  file.remove("jags_bridge.txt")
  if (!is.null(cl)) {
    parallel::stopCluster(cl = cl)
  }
  return(out)
  }
  
}
  
