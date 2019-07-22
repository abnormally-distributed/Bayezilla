#' Horseshoe
#'
#' @description This is the horseshoe model described by Carvalho et al. (2010). This tends to run very quickly
#' even for larger data sets or larger numbers of predictors and in my experience is faster and more stable (at least
#' on the tested data sets!) than the same model implemetned in Stan. \cr
#' \cr
#' \cr
#' Model Specification: \cr 
#' \cr
#' 
#' \if{html}{\figure{Horseshoe.png}{}}
#' \if{latex}{\figure{Horseshoe.png}{}}
#' \cr
#' \cr
#' Plugin Pseudo-Variances: \cr
#' \cr
#' \if{html}{\figure{pseudovar.png}{}}
#' \if{latex}{\figure{pseudovar.png}{}}
#' \cr
#'
#'
#' @references
#' Carvalho, C. M., Polson, N. G., and Scott, J. G. (2010). The horseshoe estimator for sparse signals. Biometrika, 97(2):465–480.
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
#' @param method Defaults to "rjparallel". For an alternative parallel option, choose "parallel" or. Otherwise, "rjags" (single core run).
#' @param cl Use parallel::makeCluster(# clusters) to specify clusters for the parallel methods. Defaults to two cores.
#' @param ... Other arguments to run.jags.
#'
#' @return
#' an rjags object
#' @export
#'
#' @examples
#' HS()
#'
HS = function(formula, data, family = "gaussian", log_lik = FALSE, iter = 4000, warmup=3000, adapt=3000, chains=4, thin=2, method = "rjparallel", cl = makeCluster(2), ...){
  
  X = model.matrix(formula, data)[,-1]
  y = model.frame(formula, data)[,1]
  
  if (family == "gaussian"){
  
  horseshoe = "model{
# tau is the precision, inverse of variance.
tau ~ dgamma(.01, .01) 
# lambda squared, the global penalty
global_lambda ~ dt(0, tau, 1) T(0, )
# Coefficients
Intercept ~ dnorm(0, 1e-10)
for (i in 1:P){
  local_lambda[i] ~ dt(0, 1, 1) T(0, )
  eta[i] <-  1 / (pow(global_lambda , 2) * pow(local_lambda[i], 2))
  beta[i] ~ dnorm(0, eta[i])
}
# Likelihood Function
for (i in 1:N){
  y[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
  ySim[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
  log_lik[i] <- logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]), tau)
}
Deviance <- -2 * sum(log_lik[1:N])
sigma <- sqrt(1/tau)
}"

write_lines(horseshoe, "horseshoe.txt")
jagsdata = list("y" = y, "X" = X, "N" = nrow(X), "P" = ncol(X))
inits = lapply(1:chains, function(z) list("beta" = lmSolve(formula, data)[-1], 
                                          "Intercept" = lmSolve(formula, data)[1], 
                                          "local_lambda" =  rep(1, ncol(X)),
                                          "global_lambda"= 1, 
                                          "tau" = 1,
                                          "ySim" = sample(y, length(y)),
                                          .RNG.name= "lecuyer::RngStream",
                                          .RNG.seed= sample(1:10000, 1)))
monitor = c("Intercept", "beta", "sigma", "Deviance" , "global_lambda", "local_lambda", "ySim", "log_lik")
if (log_lik == FALSE){
  monitor = monitor[-(length(monitor))]
}
out = run.jags(model = "horseshoe.txt", data = jagsdata, inits = inits, monitor = monitor, modules = c("bugs on", "glm on", "dic off"), n.chains = chains, 
               thin = thin, adapt = adapt, burnin = warmup, sample = iter, cl = cl, method = method)

file.remove("horseshoe.txt")
if (!is.null(cl)) {
  parallel::stopCluster(cl = cl)
}
return(out)

}

if (family == "binomial"){
  
  horseshoe = "model{
# lambda squared, the global penalty
global_lambda ~ dt(0, tau, 1) T(0, )
# Coefficients
Intercept ~ dnorm(0, 1e-10)
for (i in 1:P){
  local_lambda[i] ~ dt(0, 1, 1) T(0, )
  eta[i] <-  1 / (pow(global_lambda , 2) * pow(local_lambda[i], 2))
  beta[i] ~ dnorm(0, eta[i])
}
# Likelihood Function
for (i in 1:N){
      logit(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P]) 
      y[i] ~ dbern(psi[i])
      log_lik[i] <- logdensity.bern(y[i], psi[i])
      ySim[i] ~ dbern(psi[i])
}
Deviance <- -2 * sum(log_lik[1:N])
sigma <- sqrt(1/tau)
}"

write_lines(horseshoe, "horseshoe.txt")

jagsdata = list("y" = y, "X" = X, "N" = nrow(X), "P" = ncol(X), "tau" = 1 / (pow(mean(y), -1) * pow(1 - mean(y), -1)))

inits = lapply(1:chains, function(z) list("beta" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "binomial", lambda = runif(1, .01, .15), alpha = 1, standardize = FALSE))[-1,1]), 
                                          "Intercept" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "binomial", lambda = runif(1, .01, .15), alpha = 0, standardize = FALSE))[1,1]), 
                                          "local_lambda" =  rep(1, ncol(X)),
                                          "global_lambda"= 1, 
                                          "ySim" = sample(y, length(y)),
                                          .RNG.name= "lecuyer::RngStream",
                                          .RNG.seed= sample(1:10000, 1)))

monitor = c("Intercept", "beta", "Deviance" , "global_lambda", "local_lambda", "ySim", "log_lik")
if (log_lik == FALSE){
  monitor = monitor[-(length(monitor))]
}
out = run.jags(model = "horseshoe.txt", data = jagsdata, inits = inits, monitor = monitor, modules = c("bugs on", "glm on", "dic off"), n.chains = chains, 
               thin = thin, adapt = adapt, burnin = warmup, sample = iter, cl = cl, method = method)

file.remove("horseshoe.txt")
if (!is.null(cl)) {
  parallel::stopCluster(cl = cl)
}
return(out)

}

if (family == "poisson"){
  
  horseshoe = "model{
# lambda squared, the global penalty
global_lambda ~ dt(0, tau, 1) T(0, )
# Coefficients
Intercept ~ dnorm(0, 1e-10)
for (i in 1:P){
  local_lambda[i] ~ dt(0, 1, 1) T(0, )
  eta[i] <-  1 / (pow(global_lambda , 2) * pow(local_lambda[i], 2))
  beta[i] ~ dnorm(0, eta[i])
}
# Likelihood Function
for (i in 1:N){
    log(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P]) 
    y[i] ~ dpois(psi[i])
    log_lik[i] <- logdensity.pois(y[i], psi[i])
    ySim[i] ~ dpois(psi[i])
}
Deviance <- -2 * sum(log_lik[1:N])
sigma <- sqrt(1/tau)
}"

write_lines(horseshoe, "horseshoe.txt")

jagsdata = list("y" = y, "X" = X, "N" = nrow(X), "P" = ncol(X), "tau" = 1 / (pow(mean(y), -1)))

inits = lapply(1:chains, function(z) list("beta" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "poisson", lambda = runif(1, .01, .15), alpha = 1, standardize = FALSE))[-1,1]), 
                                          "Intercept" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "poisson", lambda = runif(1, .01, .15), alpha = 0, standardize = FALSE))[1,1]), 
                                          "local_lambda" =  rep(1, ncol(X)),
                                          "global_lambda"= 1, 
                                          "ySim" = sample(y, length(y)),
                                          .RNG.name= "lecuyer::RngStream",
                                          .RNG.seed= sample(1:10000, 1)))

monitor = c("Intercept", "beta", "Deviance" , "global_lambda", "local_lambda", "ySim", "log_lik")
if (log_lik == FALSE){
  monitor = monitor[-(length(monitor))]
}
out = run.jags(model = "horseshoe.txt", data = jagsdata, inits = inits, monitor = monitor, modules = c("bugs on", "glm on", "dic off"), n.chains = chains, 
               thin = thin, adapt = adapt, burnin = warmup, sample = iter, cl = cl, method = method)

file.remove("horseshoe.txt")
if (!is.null(cl)) {
  parallel::stopCluster(cl = cl)
}
return(out)

}

}
