#' Horseshoe with unpenalized design covariates
#'
#' @description This is the horseshoe model described by Carvalho et al. (2010), but with the allowance
#' for a set of covariates that are not penalized. For example, you may wish to include variables such
#' as age and gender in all models so that the coefficients for the other variables are penalized while
#' controlling for these. This is a common need in research. \cr 
#' \cr
#' This tends to run very quickly even for larger data sets or larger 
#' numbers of predictors and in my experience is faster and more stable (at least
#' on the tested data sets!) than the same model implemetned in Stan. \cr
#' \cr
#' Model Specification: \cr 
#' \cr
#' 
#' \if{html}{\figure{HorseshoeDC.png}{}}
#' \if{latex}{\figure{HorseshoeDC.png}{}}
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
#' @param design.formula formula for the design covariates.
#' @param data a data frame.
#' @param design.formula formula for the design covariates.
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
#' HSDC()
#'
HSDC = function(formula, design.formula, data, log_lik = FALSE, iter = 4000, warmup=3000, adapt=3000, chains=4, thin=2, method = "rjparallel", cl = makeCluster(2), ...){
  
  X = model.matrix(formula, data)[,-1]
  y = model.frame(formula, data)[,1]
  FX <- as.matrix(model.matrix(design.formula, data)[, -1])

if (family == "gaussian"){
    
horseshoe = 

"model{
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

# Design Variable Coefficients
for (f in 1:FP){
    design_beta[f] ~ dnorm(0, 1e-200)
}
  
# Likelihood Function
for (i in 1:N){
  mu[i] <- Intercept + sum(beta[1:P] * X[i,1:P]) + sum(design_beta[1:FP] * FX[i,1:FP]) 
  y[i] ~ dnorm(mu[i], tau)
  ySim[i] ~ dnorm(mu[i], tau)
  log_lik[i] <- logdensity.norm(y[i], mu[i], tau)
}

Deviance <- -2 * sum(log_lik[1:N])
sigma <- sqrt(1/tau)

}"

FP <- ncol(FX)
write_lines(horseshoe, "horseshoe.txt")
jagsdata = list("y" = y, "X" = X, "N" = nrow(X), "P" = ncol(X), "FP" = FP, "FX" = FX)
inits = lapply(1:chains, function(z) list("beta" = lmSolve(formula, data)[-1], 
                                          "design_beta" = as.vector(coef(lm(design.formula, data)))[-1],
                                          "Intercept" = as.vector(coef(lm(design.formula, data)))[1], 
                                          "local_lambda" = rep(1, ncol(X)),
                                          "global_lambda"= 1, 
                                          "tau" = 1,
                                          .RNG.name= "lecuyer::RngStream",
                                          .RNG.seed= sample(1:10000, 1)))

monitor = c("Intercept", "beta", "design_beta", "sigma", "Deviance", "global_lambda","local_lambda", "ySim", "log_lik")

if (log_lik == FALSE){
  monitor = monitor[-(length(monitor))]
}

out = run.jags(model = "horseshoe.txt", data = jagsdata, inits = inits, monitor = monitor, modules = c("bugs", "glm"), n.chains = chains, 
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
      logit(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P]) + sum(design_beta[1:FP] * FX[i,1:FP]) 
      y[i] ~ dbern(psi[i])
      log_lik[i] <- logdensity.bern(y[i], psi[i])
      ySim[i] ~ dbern(psi[i])
}
Deviance <- -2 * sum(log_lik[1:N])
}"

write_lines(horseshoe, "horseshoe.txt")

jagsdata = list("y" = y, "X" = X, "N" = nrow(X), "P" = ncol(X), "FP" = FP, "FX" = FX, "tau" = 1 / (pow(mean(y), -1) * pow(1 - mean(y), -1)))

inits = lapply(1:chains, function(z) list("beta" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "binomial", lambda = runif(1, .01, .15), alpha = 1, standardize = FALSE))[-1,1]), 
                                          "Intercept" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "binomial", lambda = runif(1, .01, .15), alpha = 0, standardize = FALSE))[1,1]), 
                                          "design_beta" = as.vector(coef(glm(design.formula, data, family = "binomial")))[-1],
                                          "local_lambda" =  rep(1, ncol(X)),
                                          "global_lambda"= 1, 
                                          .RNG.name= "lecuyer::RngStream",
                                          .RNG.seed= sample(1:10000, 1)))

monitor = c("Intercept", "beta", "Deviance", "global_lambda", "local_lambda", "ySim", "log_lik")

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
    log(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P]) + sum(design_beta[1:FP] * FX[i,1:FP]) 
    y[i] ~ dpois(psi[i])
    log_lik[i] <- logdensity.pois(y[i], psi[i])
    ySim[i] ~ dpois(psi[i])
}
Deviance <- -2 * sum(log_lik[1:N])
}"

write_lines(horseshoe, "horseshoe.txt")

jagsdata = list("y" = y, "X" = X, "N" = nrow(X), "P" = ncol(X), "tau" = 1 / (pow(mean(y), -1)), "FP" = FP, "FX" = FX)

inits = lapply(1:chains, function(z) list("beta" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "poisson", lambda = runif(1, .01, .15), alpha = 1, standardize = FALSE))[-1,1]), 
                                          "Intercept" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "poisson", lambda = runif(1, .01, .15), alpha = 0, standardize = FALSE))[1,1]), 
                                          "design_beta" = as.vector(coef(glm(design.formula, data, family = "poisson")))[-1],
                                          "local_lambda" =  rep(1, ncol(X)),
                                          "global_lambda"= 1, 
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

