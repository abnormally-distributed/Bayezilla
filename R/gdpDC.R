#' Generalized double pareto shrinkage prior with unpenalized design covariates
#'
#' @description The generalized double pareto shrinkage prior of Armagan, Dunson, & Lee (2013). For the binomial
#' and poisson likelihoods the same pseudovariance functions used in all other models in this package are used. 
#' This variant grants the allowance for a set of covariates that are not penalized. 
#' For example, you may wish to include variables such as age and gender in all models so that 
#' the coefficients for the other variables are penalized while controlling for these. This
#' is a common need in research. \cr\cr
#' \cr
#' This model is parameterized similarly to the Bayesian LASSO of Park & Casella (2008), Normal-Exponential-Gamma prior of
#' Griffin, J. E. and Brown, P. J. (2011), and adaptive Bayesian LASSO of Leng, Tran and David Nott (2018).
#' The key feature is that this model explicitly utilizes generalized double pareto priors through a scale mixture
#' of normals, while the Bayesian LASSO utilizes double exponential priors through a scale mixture of normals.
#' The Bayesian adaptive LASSO also utilizes double exponential just as the BLASSO, but has coefficient
#' specific shrinkage parameters. The NEG-BLASSO utilizes (as the name suggests) normal-exponential-gamma
#' priors, which behave very similarly to the GDP. Both the NEG and GDP distributions
#' have a peak at zero, just as the double exponential distribution, but have very long, student-t-like tails. \cr
#' \cr
#' In this model, the coefficient specific shrinkage parameters are given gamma distributions that with shape
#' and rate parameters (alpha and zeta, respectively) each with independent gamma(4, 8) hyperpriors, 
#' which are very concentrated near 1. \cr
#' \cr
#' To quote directly from Armagan, Dunson, & Lee (2013): \cr
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
#' robustness properties.} \cr
#' \cr
#' \cr
#' The reason for not just fixing the values at 1 is that I have observed that this does not always result
#' in much, if any, shrinkage, and that using hyperpriors results in much better sampling (better chain mixing
#' and less autocorrelation). Furthermore, hyperpriors allow the data to speak as to which values are best. 
#' 
#' \cr
#' Model Specification:
#' \cr
#' \cr
#' \if{html}{\figure{gdpDC.png}{}}
#' \if{latex}{\figure{gdpDC.png}{}}
#' \cr
#' \cr
#' The marginal probability density function for the coefficients is of the form \cr
#' \if{html}{\figure{genparetoPDF.png}{}}
#' \if{latex}{\figure{genparetoPDF.png}{}}
#' \cr 
#' Which makes the implied prior on the coefficients \cr
#' \if{html}{\figure{gdpMarginal.png}{}}
#' \if{latex}{\figure{gdpMarginal.png}{}}
#' \cr 
#'
#' @param formula the model formula
#' @param design.formula formula for the design covariates.
#' @param data a data frame.
#' @param family one of "gaussian" (the default), "binomial", or "poisson"
#' @param log_lik Should the log likelihood be monitored? The default is FALSE.
#' @param iter How many post-warmup samples? Defaults to 10000.
#' @param warmup How many warmup samples? Defaults to 2500.
#' @param adapt How many adaptation steps? Defaults to 5000.
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
#' @examples
#' gdpDC()
#' 
gdpDC = function(formula, design.formula, data, family = "gaussian", log_lik = FALSE, iter=10000, warmup=1000, adapt=2000, chains=4, thin=1, method = "parallel", cl = makeCluster(2), ...){
  
  X = model.matrix(formula, data)[,-1]
  y = model.frame(formula, data)[,1]
  FX <- as.matrix(model.matrix(design.formula, data)[, -1])
  
  
  if (family == "gaussian"){
    
  jags_gdp = "model{
  
  tau ~ dgamma(.01, .01) 
  sigma2 <- 1/tau
  alpha ~ dgamma(4, 8)
  zeta ~ dgamma(4, 8)
  for (p in 1:P){
    lambda[p] ~ dgamma(alpha , zeta)
    eta[p] ~ dexp(lambda[p]^2 / 2)
    omega[p] <- 1 / (sigma2 * eta[p])
    beta[p] ~ dnorm(0, omega[p])
  }
  
  Intercept ~ dnorm(0, 1e-10)
  
  for (f in 1:FP){
    design_beta[f] ~ dnorm(0, 1e-200)
  }
  
  for (i in 1:N){
    y[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]) + sum(design_beta[1:FP] * FX[i,1:FP]), tau)
    log_lik[i] <- logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]) + sum(design_beta[1:FP] * FX[i,1:FP]), tau)
    ySim[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]) + sum(design_beta[1:FP] * FX[i,1:FP]), tau)
  }
  
  sigma <- sqrt(1/tau)
  Deviance <- -2 * sum(log_lik[1:N])
}"
  
  P <- ncol(X)
  FP <- ncol(FX)
  write_lines(jags_gdp, "jags_gdp.txt")
  jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), FP = FP, FX = FX)
  monitor <- c("Intercept", "beta", "design_beta", "sigma", "Deviance", "alpha", "zeta", "lambda", "eta", "ySim", "log_lik")
  if (log_lik == FALSE){
    monitor = monitor[-(length(monitor))]
  }
  inits <- lapply(1:chains, function(z) list("Intercept" = lmSolve(formula, data)[1], 
                                             "beta" = lmSolve(formula, data)[-1], 
                                             "design_beta" = lmSolve(design.formula, data)[-1], 
                                             "alpha" = 1, 
                                             "zeta" = 1, 
                                             "eta" = rep(1, P), 
                                             "lambda" = rep(1, P), 
                                             "tau" = 1, 
                                             "ySim" = sample(y, length(y)),
                                             .RNG.name= "lecuyer::RngStream", 
                                             .RNG.seed= sample(1:10000, 1)))
  
  out = run.jags(model = "jags_gdp.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, n.chains = chains, ...)
  file.remove("jags_gdp.txt")
  if (!is.null(cl)) {
    parallel::stopCluster(cl = cl)
  }
  return(out)
  }
  
  
  if (family == "binomial"){
    
    jags_gdp = "model{
    
  alpha ~ dgamma(4, 8)
  zeta ~ dgamma(4, 8)
  for (p in 1:P){
    lambda[p] ~ dgamma(alpha , zeta)
    eta[p] ~ dexp(lambda[p]^2 / 2)
    omega[p] <- 1 / (sigma2 * eta[p])
    beta[p] ~ dnorm(0, omega[p])
  }
  
  Intercept ~ dnorm(0, 1e-10)
  
  for (i in 1:N){
      logit(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P]) + sum(design_beta[1:FP] * FX[i,1:FP])
      y[i] ~ dbern(psi[i])
      log_lik[i] <- logdensity.bern(y[i], psi[i])
      ySim[i] ~ dbern(psi[i])
  }
  
  sigma <- sqrt(1/tau)
  Deviance <- -2 * sum(log_lik[1:N])
}"
    
    P <- ncol(X)
    write_lines(jags_gdp, "jags_gdp.txt")
    jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X),  sigma2 = pow(mean(y), -1) * pow(1 - mean(y), -1), FP = FP, FX = FX)
    monitor <- c("Intercept", "beta", "design_beta", "Deviance", "alpha", "zeta", "lambda", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits <- lapply(1:chains, function(z) list("Intercept" = coef(glm(design.formula, data, family = "binomial"))[1], 
                                               "design_beta" = coef(glm(design.formula, data, family = "binomial"))[-1], 
                                               "alpha" = 1, 
                                               "zeta" = 1, 
                                               "eta" = rep(1, P), 
                                               "lambda" = rep(1, P), 
                                               "ySim" = sample(y, length(y)),
                                               .RNG.name= "lecuyer::RngStream", 
                                               .RNG.seed= sample(1:10000, 1)))
    
    out = run.jags(model = "jags_gdp.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, n.chains = chains, summarise = FALSE, ...)
    file.remove("jags_gdp.txt")
    if (!is.null(cl)) {
      parallel::stopCluster(cl = cl)
    }
    return(out)
  }
  
  
  if (family == "poisson"){
    
    jags_gdp = "model{
    
  alpha ~ dgamma(4, 8)
  zeta ~ dgamma(4, 8)
  for (p in 1:P){
    lambda[p] ~ dgamma(alpha , zeta)
    eta[p] ~ dexp(lambda[p]^2 / 2)
    omega[p] <- 1 / (sigma2 * eta[p])
    beta[p] ~ dnorm(0, omega[p])
  }
  
  Intercept ~ dnorm(0, 1e-10)
  
  for (i in 1:N){
    log(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P]) + sum(design_beta[1:FP] * FX[i,1:FP])
    y[i] ~ dpois(psi[i])
    log_lik[i] <- logdensity.pois(y[i], psi[i])
    ySim[i] ~ dpois(psi[i])
  }
  
  sigma <- sqrt(1/tau)
  Deviance <- -2 * sum(log_lik[1:N])
}"
    
    P <- ncol(X)
    write_lines(jags_gdp, "jags_gdp.txt")
    jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), sigma2 = pow(mean(y), -1), FP = FP, FX = FX)
    monitor <- c("Intercept", "beta", "design_beta" , "Deviance", "alpha", "zeta", "lambda", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits <- lapply(1:chains, function(z) list("Intercept" = coef(glm(design.formula, data, family = "poisson"))[1], 
                                               "design_beta" = coef(glm(design.formula, data, family = "poisson"))[-1], 
                                               "beta" = rep(0, P), 
                                               "alpha" = 1, 
                                               "zeta" = 1, 
                                               "eta" = rep(1, P), 
                                               "ySim" = sample(y, length(y)),
                                               "lambda" = rep(1, P), 
                                               .RNG.name= "lecuyer::RngStream", 
                                               .RNG.seed= sample(1:10000, 1)))
    
    out = run.jags(model = "jags_gdp.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, n.chains = chains, summarise = FALSE, ...)
    file.remove("jags_gdp.txt")
    if (!is.null(cl)) {
      parallel::stopCluster(cl = cl)
    }
    return(out)
  }
  
}  
