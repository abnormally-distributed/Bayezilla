#' Shape-adaptive shrinkage prior (SASP) 
#' 
#' 
#' @description This is the Shape-adaptive shrinkage prior (SASP) 
#' model discussed by Mutshinda and Sillanpää (2011). This utilizes a power exponential (also known as generalized
#' normal or generalized gaussian) type I (kurtotic) distribution for the coefficient priors. The power exponential
#' distribution defines a continuum of distributional shapes, with the Laplacian, the Gaussian and continuous  
#' uniform distributions arising as particular cases when kappa=1, kappa=2 and kappa = ∞, respectively. For kappa < 1
#' a cauchy-like density is defined with extremely heavy tails and infinite density at zero.
#' 
#' Henceforth, powexp(mu, omega, kappa) denotes the probability density function of a power exponential distribution 
#' with a mean of mu = 0, scale parameter omega, and shape parameter kappa. \cr 
#' 
#' \cr 
#' The power exponential distribution is not implemented in JAGS, so it is defined through a change of variables, following
#' the example of Mutshinda and Sillanpää (2011) : \cr 
#' \cr
#' u_i ~ gamma(1 / kappa, 1 / omega^kappa) \cr
#' z_i <- u_i^kappa \cr
#' xi_i ~ bernoulli(.5) \cr
#' beta_i <- -1^xi_i * z_i \cr
#' 
#' \cr
#' The variable xi is introduced so that the factor -1^xi_i yields symmetry about zero. Without this the distribution
#' would be constrained to be positive. \cr
#' \cr
#' The remainder of the model is specified below: \cr
#' \cr
#' ##### top level parameters ##### \cr
#' kappa ~ shifted-gamma(3, 6, + .25) T(, 2) \cr
#' lambda ~ gamma(.5, .1) \cr
#' tau ~ gamma(.01, .01) \cr
#' Intercept ~ normal(0, 1) \cr
#' \cr
#' ##### Coefficient specific parameters ##### \cr
#' eta_i ~ half-normal(1, 1) # individual coefficient penalties \cr
#' omega_i ~ eta_i * lambda \cr
#' \cr
#' The prior for kappa is suggested to be uniform over a small interval in the original paper. However I find that this results in
#' really poor sampling efficiency and stability. I have found replacing the original uniform with a gamma(3, 6) distribution while truncating
#' at an upper limit of 2, and then shifting it by 0.25 improves the situation. Therefore, this is what is implemented here.\cr 
#' \cr 
#' !!!!! NOTE !!!!!!! \cr
#' I have an implementation of SASP in the Stan language in the companion package to this one,
#' BayezillaPlus. I find that it can work better than the JAGS implementation since in the Stan language
#' I am able to directly specify the log probability density function for the power exponential distribution. 
#' \cr
#' @references 
#  Mutshinda, C. M., and M. J. Sillanpää (2011) Bayesian shrinkage analysis of QTLs under shape-adaptive shrinkage priors, and accurate re-estimation of genetic effects. Heredity 107: 405-412. 
#' \cr
#' @param formula the model formula
#' @param data a data frame.
#' @param family one of "gaussian", "binomial", or "poisson".
#' @param df degrees of freedom on the prior thetas 
#' @param log_lik Should the log likelihood be monitored? The default is FALSE.
#' @param iter How many post-warmup samples? Defaults to 8000.
#' @param warmup How many warmup samples? Defaults to 2500.
#' @param adapt How many adaptation steps? Defaults to 2500.
#' @param chains How many chains? Defaults to 4.
#' @param thin Thinning interval. Defaults to 4, because this model tends to be difficult to sample from and high autocorrelation is common.
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
SASP = function(formula, data, family = "gaussian", log_lik = FALSE, iter= 8000, warmup = 2500, adapt=2500, chains=4, thin = 4, method = "parallel", cl = makeCluster(3), ...){
  
  X = model.matrix(formula, data)[,-1]
  y = model.frame(formula, data)[,1]
  
  if (family == "gaussian"){
    jag_sasp = "model{
  
  Intercept ~ dnorm(0, 1)
  tau ~ dgamma(.01, .01)
  sigma <- sqrt(1/tau)
  lambda ~ dgamma(1.25 , 0.25)
  kappa_raw ~ dgamma(3, 6) T(,2)
  kappa <- kappa_raw + .25
  inv.kappa <- 1/kappa 

  for (j in 1:P){
    eta[j]~ dnorm(1, 1) T(0,) # local penalty 
    omega[j] <- lambda*eta[j] # prior scale
    pp[j]<- pow(omega[j], inv.kappa) # prior precision
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
    inits = lapply(1:chains, function(z) list("Intercept" = 0, "ySim" = y, "xi" = sample(c(0,1), size = P, replace = TRUE), "u" = rexp(P, 1), "eta" = rgamma(P, 4, 8), "lambda" = 5, "kappa_raw"= 1, "tau" = 1 , .RNG.name="lecuyer::RngStream", .RNG.seed= sample(1:10000, 1)))  
    write_lines(jag_sasp, "jag_sasp.txt")
  }
  
  if (family == "binomial"){
    jag_sasp = "model{
    
  Intercept ~ dnorm(0, 1)
  lambda ~ dgamma(1.25 , 0.25)
  kappa_raw ~ dgamma(3, 6) T(,2)
  kappa <- kappa_raw + .25
  inv.kappa <- 1/kappa 

  for (j in 1:P){
    eta[j]~ dnorm(1, 1) T(0,) # local penalty 
    omega[j] <- lambda*eta[j] # prior scale
    pp[j]<- pow(omega[j], inv.kappa) # prior precision
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
    inits = lapply(1:chains, function(z) list("Intercept" = 0, "ySim" = y, "xi" = sample(c(0,1), size = P, replace = TRUE), "u" = rexp(P, 1), "eta" = rgamma(P, 4, 8), "lambda" = 5, "kappa_raw"= 1, .RNG.name="lecuyer::RngStream", .RNG.seed= sample(1:10000, 1)))  
    write_lines(jag_sasp, "jag_sasp.txt")
  }
  
  if (family == "poisson"){
    jag_sasp = "model{
  
  Intercept ~ dnorm(0, 1)
  lambda ~ dgamma(1.25 , 0.25)
  kappa_raw ~ dgamma(3, 6) T(,2)
  kappa <- kappa_raw + .25
  inv.kappa <- 1/kappa 

  for (j in 1:P){
    eta[j]~ dnorm(1, 1) T(0,) # local penalty 
    omega[j] <- lambda*eta[j] # prior scale
    pp[j]<- pow(omega[j], inv.kappa) # prior precision
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
    inits = lapply(1:chains, function(z) list("Intercept" = 0, "ySim" = y, "xi" = sample(c(0,1), size = P, replace = TRUE), "u" = rexp(P, 1), "eta" = rgamma(P, 4, 8), "lambda" = 5, "kappa_raw"= 1, .RNG.name="lecuyer::RngStream", .RNG.seed= sample(1:10000, 1)))  
    write_lines(jag_sasp, "jag_sasp.txt")
    
  }
  
  out = run.jags(model = "jag_sasp.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, n.chains = chains, method = method, cl = cl, summarise = FALSE, ...)
  file.remove("jag_sasp.txt")
  if (is.null(cl) == FALSE){
    parallel::stopCluster(cl = cl)
  }
  return(out)
}