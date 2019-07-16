#' Shape Adaptive Shrinkage Prior with Unpenalized Design Covariates
#'
#' @description The shape adaptive shrinkage prior of Sillanpää &  Mutshinda (2011). This is essentially bridge regression, but with a
#' shape parameter describing the Lp norm that is allowed to vary rather than stay fixed at a single value. The generalized gaussian 
#' prior is parameterized in the manner of Mallick, H. & Yi (2018) rather than the method originally described in Sillanpää &  Mutshinda (2011).
#' Analytically, this makes no difference, but computationally, it is much faster and more stable. This function has the further allowance for a set of covariates that are not penalized. 
#' For example, you may wish to include variables such as age and gender so that  the coefficients for the other variables are 
#' penalized while controlling for these. This is a common need in research.
#' \cr
#'
#' \cr
#' Model Specification: \cr
#' \cr
#' \if{html}{\figure{saspDC.png}{}}
#' \if{latex}{\figure{saspDC.png}{}}
#' \cr
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
#' Mallick, H. & Yi, N. (2018) Bayesian bridge regression, Journal of Applied Statistics, 45:6, 988-1008, DOI: 10.1080/02664763.2017.1324565 \cr
#' \cr
#' Sillanpää, S., &  Mutshinda, C., (2011) Bayesian shrinkage analysis of QTLs under shape-adaptive shrinkage priors, and accurate re-estimation of genetic effects. Heredity volume 107, pages 405–412. doi: 10.1038/hdy.2011.37 \cr
#' \cr
#' 
#' 
#' @return
#' a runjags object
#' @export
#'
#' @examples
#' saspDC()
#' 
saspDC = function(formula, design.formula, data, family = "gaussian", log_lik = FALSE, iter=10000, warmup=1000, adapt=2000, chains=4, thin=1, method = "parallel", cl = makeCluster(2), ...){
  
  X = model.matrix(formula, data)[,-1]
  y = model.frame(formula, data)[,1]
  FX <- as.matrix(model.matrix(design.formula, data)[, -1])
  
  if (family == "gaussian"){
    
    jags_bridge = "model{
  
  tau ~ dgamma(.01, .01) 
  sigma <- sqrt(1/tau)
  lambda ~ dgamma(0.50, 0.20)
  kappa ~ dgamma(3.0625, 2.1875)
  
  for (i in 1:P){
    u[i] ~ dgamma( (1/kappa) + 1  , lambda )
    beta[i] ~ dunif(-1 * pow(sigma * u[i], 1/kappa), pow(sigma * u[i], 1/kappa))
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
  
  Deviance <- -2 * sum(log_lik[1:N])
}"
    
    P <- ncol(X)
    write_lines(jags_bridge, "jags_bridge.txt")
    jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), FP = FP, FX = FX)
    monitor <- c("Intercept", "beta",  "design_beta", "sigma", "kappa", "lambda", "Deviance", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits <- lapply(1:chains, function(z) list("Intercept" = lmSolve(formula, data)[1], 
                                               "beta" = rep(0, P), 
                                               "design_beta" = coef(lm(design.formula ~ ., data))[-1],
                                               "u" = rgamma(P, (1 / .75) + 1, 1), 
                                               "lambda" = 1, 
                                               "kappa" = .75, 
                                               "tau" = 1, 
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
  
  
  if (family == "binomial" || family == "logistic"){
    
    jags_bridge = "model{


  lambda ~ dgamma(0.50, 0.20)
  kappa ~ dgamma(3.0625, 2.1875)
    
  for (f in 1:FP){
      design_beta[f] ~ dnorm(0, 1e-200)
  }

  for (i in 1:P){
    u[i] ~ dgamma( (1/kappa) + 1  , lambda )
    beta[i] ~ dunif(-1 * pow(sigma * u[i], 1/kappa), pow(sigma * u[i], 1/kappa))
  }
  
  Intercept ~ dnorm(0, 1e-10)
  
    for (i in 1:N){
      logit(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P]) + sum(design_beta[1:FP] * FX[i,1:FP])
      y[i] ~ dbern(psi[i])
      log_lik[i] <- logdensity.bern(y[i], psi[i])
      ySim[i] ~ dbern(psi[i])
    }

  Deviance <- -2 * sum(log_lik[1:N])
}"
    
    P <- ncol(X)
    write_lines(jags_bridge, "jags_bridge.txt")
    jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), FP = FP, FX = FX, sigma = sqrt(pow(mean(y), -1) * pow(1 - mean(y), -1)))
    monitor <- c("Intercept", "beta",  "design_beta", "kappa", "lambda", "Deviance", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits <- lapply(1:chains, function(z) list("Intercept" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "binomial", lambda = runif(1, .01, .15), alpha = 0, standardize = FALSE))[1,1]), 
                                               "beta" = rep(0, P), 
                                               "design_beta" = as.vector(coef(glm(design.formula ~ data, family = "binomial")))[-1],
                                               "u" = rgamma(P, (1 / .75) + 1, 1), 
                                               "lambda" = 1, 
                                               "kappa" = .75, 
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


  lambda ~ dgamma(0.50, 0.20)
  kappa ~ dgamma(3.0625, 2.1875)
  
  for (i in 1:P){
    u[i] ~ dgamma( (1/kappa) + 1  , lambda )
    beta[i] ~ dunif(-1 * pow(sigma * u[i], 1/kappa), pow(sigma * u[i], 1/kappa))
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
  
  P <- ncol(X)
  write_lines(jags_bridge, "jags_bridge.txt")
  jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), FP = FP, FX = FX, sigma = sqrt(pow(mean(y) , -1)))
  monitor <- c("Intercept", "beta",  "design_beta", "kappa", "lambda", "Deviance", "ySim", "log_lik")
  if (log_lik == FALSE){
    monitor = monitor[-(length(monitor))]
  }
  inits <- lapply(1:chains, function(z) list("Intercept" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "poisson", lambda = runif(1, .01, .15), alpha = 0, standardize = FALSE))[1,1]), 
                                             "beta" = rep(0, P), 
                                             "design_beta" = as.vector(coef(glm(design.formula ~ data, family = "poisson")))[-1],
                                             "u" = rgamma(P, (1 / .75) + 1, 1), 
                                             "lambda" = 1, 
                                             "kappa" = .75, 
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

