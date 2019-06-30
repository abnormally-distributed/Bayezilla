#' Horseshoe+
#'
#' @description This is the horseshoe+ model described by Bhadra et al (2017). This is essentially the horseshoe homologue
#' of the extended LASSO. The specification is identical to the horseshoe except an additional local shrinkage parameter is
#' added to the model. Note that in the model I use somewhat different names for the paramters. This is because a key parameter
#' in the model is often called "tau" in the literature (often geared towards Stan), but JAGS uses the precision parameterization
#' of the normal distribution. The precision is denoted as "tau", so this warranted some new terms. See the model specification below. \cr
#' \cr
#' Model Specification: \cr 
#' \cr
#' 
#' \if{html}{\figure{HorseshoePlus.png}{}}
#' \if{latex}{\figure{HorseshoePlus.png}{}}
#'
#'
#' @references
#' Bhadra, Anindya; Datta, Jyotishka; Polson, Nicholas G.; Willard, Brandon. The Horseshoe+ Estimator of Ultra-Sparse Signals. Bayesian Anal. 12 (2017), no. 4, 1105--1131. doi:10.1214/16-BA1028. https://projecteuclid.org/euclid.ba/1474572263 \cr
#' \cr
#' Carvalho, C. M., Polson, N. G., and Scott, J. G. (2010). The horseshoe estimator for sparse signals. Biometrika, 97(2):465–480. \cr
#'
#' @param formula the model formula
#' @param data a data frame.
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
#' HSplus()
#'
HSplus = function(formula, data, log_lik = FALSE, iter = 4000, warmup=3000, adapt=3000, chains=4, thin=2, method = "rjparallel", cl = makeCluster(2), ...){

  X = model.matrix(formula, data)[,-1]
  y = model.frame(formula, data)[,1]

  horseshoePlus = "model{
# tau is the precision, inverse of variance.
tau ~ dgamma(.01, .01) 
# lambda squared, the global penalty
global_lambda ~ dt(0, tau, 1) T(0, )
# Coefficients
Intercept ~ dnorm(0, 1e-10)
for (i in 1:P){
  local_lambda_A[i] ~ dt(0, 1, 1) T(0, )
  local_lambda_B[i] ~ dt(0, global_lambda * local_lambda_A[i], 1) T(0, )
  eta[i] <- 1 / pow(local_lambda_B[i], 2)
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

write_lines(horseshoePlus, "horseshoePlus.txt")
jagsdata = list("y" = y, "X" = X, "N" = nrow(X), "P" = ncol(X))
inits = lapply(1:chains, function(z) list("beta" = rep(0, ncol(X)),
                                          "Intercept" = 0,
                                          "local_lambda_A" =  rep(1, ncol(X)),
                                          "local_lambda_B" =  rep(1, ncol(X)),
                                          "global_lambda"= 1,
                                          "tau" = 1,
                                          .RNG.name= "lecuyer::RngStream",
                                          .RNG.seed= sample(1:10000, 1)))
monitor = c("Intercept", "beta", "sigma", "Deviance" , "global_lambda", "local_lambda_A", "local_lambda_B", "ySim", "log_lik")
if (log_lik == FALSE){
  monitor = monitor[-(length(monitor))]
}
out = run.jags(model = "horseshoePlus.txt", data = jagsdata, inits = inits, monitor = monitor, modules = c("bugs on", "glm on", "dic off"), n.chains = chains,
               thin = thin, adapt = adapt, burnin = warmup, sample = iter, cl = cl, method = method)

file.remove("horseshoePlus.txt")
if (!is.null(cl)) {
  parallel::stopCluster(cl = cl)
}
return(out)
}
