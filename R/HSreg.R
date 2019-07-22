#' Regularized Horseshoe
#'
#' @description This is the horseshoe model described by Piironen & Vehtari (2017). This tends to run very quickly
#' even for larger data sets or larger numbers of predictors and in my experience is faster and more stable (at least
#' on the tested data sets!) than the same model implemetned in Stan. If the horseshoe+ is analagous to the 
#' adaptive Bayesian LASSO, then this could be compared to the Bayesian Elastic Net in that it imposes a combination
#' of different shrinkage penalties (the elastic net being a combination of L1 and L2, and the regularized horseshoe being
#' a combination of sub-L1 and student-t penalties). \cr
#' \cr
#' Model Specification: \cr 
#' \cr
#' 
#' \if{html}{\figure{regularizedHorseshoe.png}{}}
#' \if{latex}{\figure{regularizedHorseshoe.png}{}}
#' \cr
#' \cr
#' Plugin Pseudo-Variances: \cr
#' \cr
#' \if{html}{\figure{pseudovar.png}{}}
#' \if{latex}{\figure{pseudovar.png}{}}
#' \cr
#'
#' @references
#' Piironen, Juho; Vehtari, Aki. Sparsity information and regularization in the horseshoe and other shrinkage priors. Electron. J. Statist. 11 (2017), no. 2, 5018--5051. doi:10.1214/17-EJS1337SI. https://projecteuclid.org/euclid.ejs/1513306866 \cr
#' \cr
#' Carvalho, C. M., Polson, N. G., and Scott, J. G. (2010). The horseshoe estimator for sparse signals. Biometrika, 97(2):465–480. \cr
#'
#' @param formula the model formula
#' @param data a data frame.
#' @param family one of "gaussian", "binomial", or "poisson".
#' @param phi your prior guess on the inclusion probability. Defaults to .50. Best way to come up with a figure is a prior guess on how many coefficients are non-zero out of the total number of predictors.
#' @param slab_scale the standard deviation of the "slab". Defaults to 2.
#' @param slab_df the degrees of freedom fo the slab. Higher degrees of freedom give increased L2-like regularization. Defaults to 3.
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
#' HSreg()
#'
HSreg = function(formula, data, family = "gaussian", phi = .50, slab_scale = 2, slab_df = 3, log_lik = FALSE, 
                 iter = 10000, warmup = 4000, adapt = 5000, chains=4, thin=1, method = "parallel", cl = makeCluster(2), ...){
  
  X = model.matrix(formula, data)[,-1]
  y = model.frame(formula, data)[,1]
  
  if (family == "gaussian"){  
  
  horseshoe = "model{

# tau is the precision, inverse of variance.
tau ~ dgamma(.01, .01) 
sigma <- sqrt(1/tau)

# lambda squared, the global penalty
lambda0sqrd <- pow((phi / (1-phi)) * (sigma/sqrt(N)), 2)
lambda ~ dt(0, 1/lambda0sqrd, 1) T(0, )

# control parameter
c2_inv ~ dgamma(df / 2, (df*prior_variance)/2) 
c2 <- 1 / c2_inv

# Coefficients
Intercept ~ dnorm(0, 1e-10)
for (i in 1:P){
  eta[i] ~ dt(0,1,1) T(0, ) 
  eta_tilde[i] <- (pow(eta[i], 2)*c2) / (c2+(pow(lambda,2)*pow(eta[i], 2)))
  eta_inv[i] <- 1 / (eta_tilde[i] * lambda)
  beta[i] ~ dnorm(0, eta_inv[i])
}

# Likelihood Function
for (i in 1:N){
  y[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
  ySim[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
  log_lik[i] <- logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]), tau)
}
Deviance <- -2 * sum(log_lik[1:N])

}"

write_lines(horseshoe, "horseshoe.txt")
jagsdata = list("y" = y, "X" = X, "N" = nrow(X), "P" = ncol(X), "phi" = phi, "df" = slab_df, "prior_variance" = square(slab_scale))
inits = lapply(1:chains, function(z) list("beta" = rep(0, ncol(X)), 
                                          "Intercept" = 0, 
                                          "eta" =  rep(1, ncol(X)),
                                          "c2_inv" = 1/slab_df, 
                                          "lambda"= 10, 
                                          "tau" = 1,
                                          "ySim" = sample(y, length(y)),
                                          .RNG.name= "lecuyer::RngStream",
                                          .RNG.seed= sample(1:10000, 1)))
monitor = c("Intercept", "beta", "sigma", "Deviance", "lambda", "c2", "ySim", "log_lik")
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
lambda0sqrd <- pow((phi / (1-phi)) * (sigma/sqrt(N)), 2)
lambda ~ dt(0, 1/lambda0sqrd, 1) T(0, )

# control parameter
c2_inv ~ dgamma(df / 2, (df*prior_variance)/2) 
c2 <- 1 / c2_inv

# Coefficients
Intercept ~ dnorm(0, 1e-10)
for (i in 1:P){
  eta[i] ~ dt(0,1,1) T(0, ) 
  eta_tilde[i] <- (pow(eta[i], 2)*c2) / (c2+(pow(lambda,2)*pow(eta[i], 2)))
  eta_inv[i] <- 1 / (eta_tilde[i] * lambda)
  beta[i] ~ dnorm(0, eta_inv[i])
}

# Likelihood Function
for (i in 1:N){
      logit(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P]) 
      y[i] ~ dbern(psi[i])
      log_lik[i] <- logdensity.bern(y[i], psi[i])
      ySim[i] ~ dbern(psi[i])
}
Deviance <- -2 * sum(log_lik[1:N])

}"

write_lines(horseshoe, "horseshoe.txt")
jagsdata = list("y" = y, "X" = X, "N" = nrow(X), "P" = ncol(X), "phi" = phi, "df" = slab_df, "prior_variance" = square(slab_scale), sigma = sqrt(pow(mean(y), -1) * pow(1 - mean(y), -1)))
inits = lapply(1:chains, function(z) list("beta" = rep(0, ncol(X)), 
                                          "Intercept" = 0, 
                                          "eta" =  rep(1, ncol(X)),
                                          "c2_inv" = 1/slab_df, 
                                          "lambda"= 10, 
                                          "ySim" = sample(y, length(y)),
                                          .RNG.name= "lecuyer::RngStream",
                                          .RNG.seed= sample(1:10000, 1)))

monitor = c("Intercept", "beta", "Deviance" , "lambda", "c2", "ySim", "log_lik")
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
lambda0sqrd <- pow((phi / (1-phi)) * (sigma/sqrt(N)), 2)
lambda ~ dt(0, 1/lambda0sqrd, 1) T(0, )

# control parameter
c2_inv ~ dgamma(df / 2, (df*prior_variance)/2) 
c2 <- 1 / c2_inv

# Coefficients
Intercept ~ dnorm(0, 1e-10)
for (i in 1:P){
  eta[i] ~ dt(0,1,1) T(0, ) 
  eta_tilde[i] <- (pow(eta[i], 2)*c2) / (c2+(pow(lambda,2)*pow(eta[i], 2)))
  eta_inv[i] <- 1 / (eta_tilde[i] * lambda)
  beta[i] ~ dnorm(0, eta_inv[i])
}

# Likelihood Function
for (i in 1:N){
      logit(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P]) 
      y[i] ~ dbern(psi[i])
      log_lik[i] <- logdensity.bern(y[i], psi[i])
      ySim[i] ~ dbern(psi[i])
}
Deviance <- -2 * sum(log_lik[1:N])

}"

write_lines(horseshoe, "horseshoe.txt")
jagsdata = list("y" = y, "X" = X, "N" = nrow(X), "P" = ncol(X), "phi" = phi, "df" = slab_df, "prior_variance" = square(slab_scale), sigma = sqrt(pow(mean(y) , -1)))
inits = lapply(1:chains, function(z) list("beta" = rep(0, ncol(X)), 
                                          "Intercept" = 0, 
                                          "eta" =  rep(1, ncol(X)),
                                          "c2_inv" = 1/slab_df, 
                                          "lambda"= 10, 
                                          "ySim" = sample(y, length(y)),
                                          .RNG.name= "lecuyer::RngStream",
                                          .RNG.seed= sample(1:10000, 1)))

monitor = c("Intercept", "beta", "Deviance", "lambda", "c2", "ySim", "log_lik")
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