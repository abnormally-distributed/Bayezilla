#' Stochastic Search Variable Selection (Indicator Variable and Adaptive Student-t model)
#'
#' @description 
#' This is the \strong{I}ndicator variables and \strong{A}daptive Student’s \strong{t}-distributions (IAt)
#' model discussed by Knürr, Läärä, and Sillanpää (2011). 
#' \cr
#' \cr
#' Model Specification:
#' \cr
#' \if{html}{\figure{IAt.png}{}}
#' \if{latex}{\figure{IAt.png}{}}
#'   
#' @references 
#  Knürr, T., E. Läärä, and M. J. Sillanpää (2011) Genetic analysis of complex traits via Bayesian variable selection: the utility of a mixture of uniform priors. Genetics Research 93: 303-318. doi:10.1017/S0016672311000164
#' 
#' @param formula the model formula
#' @param data a data frame.
#' @param family one of "gaussian", "binomial", or "poisson".
#' @param log_lik Should the log likelihood be monitored? The default is FALSE.
#' @param iter How many post-warmup samples? Defaults to 5000.
#' @param warmup How many warmup samples? Defaults to 5000.
#' @param adapt How many adaptation steps? Defaults to 10000.
#' @param chains How many chains? Defaults to 4.
#' @param thin Thinning interval. Defaults to 2.
#' @param method Defaults to "parallel". For an alternative parallel option, choose "rjparallel" or. Otherwise, "rjags" (single core run).
#' @param cl Use parallel::makeCluster(# clusters) to specify clusters for the parallel methods. Defaults to two cores.
#' @param ... Other arguments to run.jags.
#'
#' @return A run.jags object
#' @export
#'
#' @examples
#' IAt()
#'
IAt = function(formula, data, family = "gaussian", log_lik = FALSE, iter= 10000, warmup = 5000, adapt = 5000, chains=4, thin = 2, method = "parallel", cl = makeCluster(2), ...){
  
  X = model.matrix(formula, data)[,-1]
  y = model.frame(formula, data)[,1]
  
if (family == "gaussian"){
  jags_iat = "model{
  
  tau ~ dgamma(.01, .01)
  Intercept ~ dnorm(0, 1e-10)
  phi ~ dbeta(1, 1)
  df ~ dgamma(1, 0.33333333333333333333)
  
  for (p in 1:P){
    beta[p] <- delta[p]*theta[p]
    theta[p] ~ dnorm(0, omega[p])
    omega[p] ~ dgamma(df / 2, df / 2)
    delta[p] ~ dbern(phi)
  }
  
  for (i in 1:N){
    y[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
    ySim[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
    log_lik[i] <- logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]), tau)
  }
  sigma <- sqrt(1/tau)
  Deviance <- -2 * sum(log_lik[1:N])
  BIC <- (log(N) * sum(delta[1:P])) + Deviance
}"
  
  P = ncol(X)
  monitor = c("Intercept", "beta", "sigma", "df", "phi" , "delta", "Deviance", "BIC",  "ySim", "log_lik")
  if (log_lik == FALSE){
    monitor = monitor[-(length(monitor))]
  }
  jagsdata = list(X = X, y = y, N = length(y), P = ncol(X))
  inits = lapply(1:chains, function(z) list("Intercept" = 0, "ySim" = y, "omega" = rep(1, P), "df" = runif(1, 3, 30), "theta" = rep(0, P),  "phi" = .5, tau = 1 , .RNG.name="lecuyer::RngStream", .RNG.seed= sample(1:10000, 1)))  
  write_lines(jags_iat, "jags_iat.txt")
}

if (family == "binomial"){
  jags_iat = "model{
  
  Intercept ~ dnorm(0, 1e-10)
  phi ~ dbeta(1, 1)
  df ~ dgamma(1, 0.33333333333333333333)

  for (p in 1:P){
    beta[p] <- delta[p]*theta[p]
    theta[p] ~ dnorm(0, omega[p])
    omega[p] ~ dgamma(df / 2, df / 2)
    delta[p] ~ dbern(phi)
  }

  for (i in 1:N){
    logit(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
    y[i] ~ dbern(psi[i])
    log_lik[i] <- logdensity.bern(y[i], psi[i])
    ySim[i] ~ dbern(psi[i])
   }
   Deviance <- -2 * sum(log_lik[1:N])
   BIC <- (log(N) * sum(delta[1:P])) + Deviance
}"
  
  P = ncol(X)
  monitor = c("Intercept", "beta", "df", "phi", "delta", "Deviance", "BIC", "ySim" ,"log_lik")
  if (log_lik == FALSE){
    monitor = monitor[-(length(monitor))]
  }
  jagsdata = list(X = X, y = y, N = length(y), P = ncol(X))
  inits = lapply(1:chains, function(z) list("Intercept" = 0, "ySim" = y, "omega" = rep(1, P), "df" = runif(1, 3, 30), "theta" = rep(0, P),  "phi" = .5, .RNG.name="lecuyer::RngStream", .RNG.seed= sample(1:10000, 1)))  
  write_lines(jags_iat, "jags_iat.txt")
}

if (family == "poisson"){
  jags_iat = "model{
  
  Intercept ~ dnorm(0, 1e-10)
  phi ~ dbeta(1, 1)
  df ~ dgamma(1, 0.33333333333333333333)

  for (p in 1:P){
    beta[p] <- delta[p]*theta[p]
    theta[p] ~ dnorm(0, omega[p])
    omega[p] ~ dgamma(df / 2, df / 2)
    delta[p] ~ dbern(phi)
  }

  for (i in 1:N){
    log(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
    y[i] ~ dpois(psi[i])
    log_lik[i] <- logdensity.pois(y[i], psi[i])
    ySim[i] ~ dpois(psi[i])
   }
   Deviance <- -2 * sum(log_lik[1:N])
   BIC <- (log(N) * sum(delta[1:P])) + Deviance
}"
  
  P = ncol(X)
  monitor = c("Intercept", "beta", "df", "phi", "delta", "Deviance", "BIC", "ySim" ,"log_lik")
  if (log_lik == FALSE){
    monitor = monitor[-(length(monitor))]
  }
  jagsdata = list(X = X, y = y, N = length(y), P = ncol(X))
  inits = lapply(1:chains, function(z) list("Intercept" = 0, "ySim" = y, "omega" = rep(1, P), "df" = runif(1, 3, 30), "theta" = rep(0, P),  "phi" = .5, .RNG.name="lecuyer::RngStream", .RNG.seed= sample(1:10000, 1)))
  write_lines(jags_iat, "jags_iat.txt")

}

out = run.jags(model = "jags_iat.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, n.chains = chains, method = method, cl = cl, summarise = FALSE, ...)
file.remove("jags_iat.txt")
if (is.null(cl) == FALSE){
  parallel::stopCluster(cl = cl)
}
return(out)
}