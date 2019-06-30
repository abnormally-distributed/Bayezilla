#' Group+Within Group Selection with the Group Horseshoe+
#'
#' @description This is the horseshoe+ model described by Bhadra et al (2017) adapted to the problem of group selection, and is
#' akin to the group Bayesian LASSO. The resulting model has three levels of shrinkage: a global shrinkage level, group level shrinkage,
#' and finally, coefficient level shrinkage. This is probably best applied when there are a rather large number of variables that can
#' be clustered into groups that themselves have fairly large numbers (at least 5-6 variables each). \cr
#' \cr
#' Note that this is a novel implementation of the horseshoe. 
#' While this does seem to work well in tests, it is not a canonical model so to speak.
#' Nevertheless results from tests seem to indicate the horseshoe family of models benefits
#' from group structure just as the LASSO family of models does. 
#' \cr
#' \cr
#' Model Specification: \cr 
#' \cr
#' 
#' \if{html}{\figure{groupHorseshoePlus.png}{}}
#' \if{latex}{\figure{groupHorseshoePlus.png}{}}
#'
#'
#' @references
#' Bhadra, Anindya; Datta, Jyotishka; Polson, Nicholas G.; Willard, Brandon. The Horseshoe+ Estimator of Ultra-Sparse Signals. Bayesian Anal. 12 (2017), no. 4, 1105--1131. doi:10.1214/16-BA1028. https://projecteuclid.org/euclid.ba/1474572263 \cr
#' \cr
#' Carvalho, C. M., Polson, N. G., and Scott, J. G. (2010). The horseshoe estimator for sparse signals. Biometrika, 97(2):465–480. \cr
#' \cr
#' Yuan, Ming; Lin, Yi (2006). Model Selection and Estimation in Regression with Grouped Variables. Journal of the Royal Statistical Society. Series B (statistical Methodology). Wiley. 68 (1): 49–67. doi:10.1111/j.1467-9868.2005.00532.x \cr
#' 
#'
#'
#' @param X the model matrix. Construct this manually with model.matrix()[,-1]
#' @param y the outcome variable
#' @param idx the group labels. Should be of length = to ncol(model.matrix()[,-1]) with the group assignments
#' for each covariate. Please ensure that you start numbering with 1, and not 0.
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
#' groupHSplus()
#'
groupHSplus = function(X, y, idx, log_lik = FALSE, iter = 4000, warmup=3000, adapt=3000, chains=4, thin=2, method = "rjparallel", cl = makeCluster(2), ...){
  

  horseshoePlus = "model{
# tau is the precision, inverse of variance.
tau ~ dgamma(.01, .01) 

# lambda squared, the global penalty
global_lambda ~ dt(0, tau, 1) T(0, )

# group lambda, the group level penalities
for (g in 1:nG){
  group_lambda[g] ~ dt(0, 1, 1) T(0, )
}

# Coefficients
Intercept ~ dnorm(0, 1e-10)
for (i in 1:P){
  local_lambda[i] ~ dt(0, global_lambda * group_lambda[idx[i]], 1) T(0, )
  eta[i] <- 1 / pow(local_lambda[i], 2)
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
jagsdata = list("y" = y, "X" = X, "N" = nrow(X), "P" = ncol(X), "idx" = idx, "nG" = max(idx))
inits = lapply(1:chains, function(z) list("beta" = rep(0, ncol(X)),
                                          "Intercept" = 0,
                                          "group_lambda" =  rep(1, max(idx)),
                                          "local_lambda" =  rep(1, ncol(X)),
                                          "global_lambda"= 1,
                                          "tau" = 1,
                                          .RNG.name= "lecuyer::RngStream",
                                          .RNG.seed= sample(1:10000, 1)))
monitor = c("Intercept", "beta", "sigma", "Deviance" , "global_lambda", "group_lambda", "local_lambda", "ySim", "log_lik")
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
