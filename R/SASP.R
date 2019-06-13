#' Shape-adaptive shrinkage prior (SASP) 
#' 
#' 
#' @description This is the \Shape-adaptive shrinkage prior (SASP) 
#' model discussed by Mutshinda and Sillanpää (2011). 
#' \cr
#' 
#' Another model to consider is the  \code{\link[Bayezilla]{apcSpike}} (adaptive powered correlation prior) or 
#' Sillanpää's other variable selection algorithms, the extended Bayesian Lasso and IAt models
#' provided in this package, respectively at \code{\link[Bayezilla]{extLASSO}} and \code{\link[Bayezilla]{IAt}} \cr
#'
#' @references 
#  Mutshinda, C. M., and M. J. Sillanpää (2011) Bayesian shrinkage analysis of QTLs under shape-adaptive shrinkage priors, and accurate re-estimation of genetic effects. Heredity 107: 405-412. 
#' \cr
#' @param formula the model formula
#' @param data a data frame.
#' @param family one of "gaussian", "binomial", or "poisson".
#' @param df degrees of freedom on the prior thetas 
#' @param log_lik Should the log likelihood be monitored? The default is FALSE.
#' @param iter How many post-warmup samples? Defaults to 10000.
#' @param warmup How many warmup samples? Defaults to 5000.
#' @param adapt How many adaptation steps? Defaults to 15000.
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
#' SASP()
#'
SASP = function(formula, data, family = "gaussian", log_lik = FALSE, iter= 4000, warmup = 1500, adapt=5000, chains=4, thin = 1, method = "parallel", cl = makeCluster(3), ...){
  
  X = model.matrix(formula, data)[,-1]
  y = model.frame(formula, data)[,1]
  
  if (family == "gaussian"){
    jag_sasp = "model{
  
  Intercept ~ dnorm(0, 1)
  tau ~ dgamma(.1, .1)
  sigma <- sqrt(1/tau)
  lambda ~ dgamma(0.25, 0.05)
  kappa ~ dunif(0,2)
  inv.kappa <- 1/kappa 

  for (j in 1:P){
    eta[j]~ dnorm(1, 1) T(0,)
    delta[j] <- lambda*eta[j]
    pp[j]<- pow(delta[j], inv.kappa)
    u[j] ~ dgamma(inv.kappa, pp[j])
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
    monitor = c("Intercept", "beta", "sigma", "Deviance", "kappa", "lambda", "ySim" ,"log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    jagsdata = list(X = X, y = y, N = length(y), P = ncol(X))
    inits = lapply(1:chains, function(z) list("Intercept" = 0, "ySim" = y, "xi" = rep(1, P), "u" = rep(1, P), "eta" = rep(1, P), "lambda" = 10, "kappa"= 1, "tau" = 1 , .RNG.name="lecuyer::RngStream", .RNG.seed= sample(1:10000, 1)))  
    write_lines(jag_sasp, "jag_sasp.txt")
  }
  
  if (family == "binomial"){
    jag_sasp = "model{
    
  Intercept ~ dnorm(0, 1)
  lambda ~ dgamma(0.5, 0.1)
  kappa ~ dunif(0.75, 2.15)
  inv.kappa <- 1/kappa

  for (j in 1:P){
    eta[j]~ dnorm(1, 1) T(0,)
    delta[j] <- lambda*eta[j]
    pp[j]<- pow(delta[j], inv.kappa)
    u[j] ~ dgamma(inv.kappa, pp[j])
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
    monitor = c("Intercept", "beta", "Deviance", "kappa", "lambda", "ySim" ,"log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    jagsdata = list(X = X, y = y, N = length(y), P = ncol(X))
    inits = lapply(1:chains, function(z) list("Intercept" = 0, "ySim" = y, "xi" = rep(1, P), "u" = rep(1, P), "eta" = rep(1, P), "lambda" = .5, "kappa"= 1, "tau" = 1 , .RNG.name="lecuyer::RngStream", .RNG.seed= sample(1:10000, 1)))  
    write_lines(jag_sasp, "jag_sasp.txt")
  }
  
  if (family == "poisson"){
    jag_sasp = "model{
  
  Intercept ~ dnorm(0, 2)
  lambda ~ dgamma(0.5, 0.1)
  kappa ~ dunif(0.75, 2.5)
  inv.kappa <- 1/kappa

  for (j in 1:P){
    eta[j]~ dnorm(1, 1) T(0,)
    delta[j] <- lambda*eta[j]
    pp[j]<- pow(delta[j], inv.kappa)
    u[j] ~ dgamma(inv.kappa, pp[j])
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
    monitor = c("Intercept", "beta", "Deviance", "kappa", "lambda", "ySim" ,"log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    jagsdata = list(X = X, y = y, N = length(y), P = ncol(X))
    inits = lapply(1:chains, function(z) list("Intercept" = 0, "ySim" = y, "xi" = rep(1, P), "u" = rep(1, P), "eta" = rep(1, P), "lambda" = .5, "kappa"= 1, "tau" = 1 , .RNG.name="lecuyer::RngStream", .RNG.seed= sample(1:10000, 1)))  
    write_lines(jag_sasp, "jag_sasp.txt")
    
  }
  
  out = run.jags(model = "jag_sasp.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, n.chains = chains, method = method, cl = cl, summarise = FALSE, ...)
  file.remove("jag_sasp.txt")
  if (is.null(cl) == FALSE){
    parallel::stopCluster(cl = cl)
  }
  return(out)
}