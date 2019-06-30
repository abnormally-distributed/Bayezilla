#' Group Horseshoe
#'
#' @description This is the horseshoe model described by Carvalho et al. (2010) adapted for group selection, akin to the
#' group Bayesian LASSO. Note that this is a novel implementation of the horseshoe. 
#' While this does seem to work well in tests, it is not a canonical model so to speak.
#' Nevertheless results from tests seem to indicate the horseshoe family of models benefits
#' from group structure just as the LASSO family of models does. 
#' \cr
#' \cr
#' Model Specification: \cr 
#' \cr
#' 
#' \if{html}{\figure{groupHorseshoe.png}{}}
#' \if{latex}{\figure{groupHorseshoe.png}{}}
#' 
#' @references
#' Carvalho, C. M., Polson, N. G., and Scott, J. G. (2010). The horseshoe estimator for sparse signals. Biometrika, 97(2):465–480. \cr
#' \cr
#' Yuan, Ming; Lin, Yi (2006). Model Selection and Estimation in Regression with Grouped Variables. Journal of the Royal Statistical Society. Series B (statistical Methodology). Wiley. 68 (1): 49–67. doi:10.1111/j.1467-9868.2005.00532.x \cr
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
#' groupHS()
#'
groupHS = function(X, y, idx, log_lik = FALSE, iter = 4000, warmup=3000, adapt=3000, chains=4, thin=2, method = "rjparallel", cl = makeCluster(2), ...){
  
  
horseshoe = "model{
  # tau is the precision, inverse of variance.
  tau ~ dgamma(.01, .01) 

  # lambda squared, the global penalty
  global_lambda ~ dt(0, tau, 1) T(0, )

  # Group Level shrinkage
  for (g in 1:nG){
    group_lambda[g] ~ dt(0, 1, 1) T(0, )
    eta[g] <-  1 / (pow(global_lambda , 2) * pow(group_lambda[g], 2))
  }
  
  # Coefficients
  Intercept ~ dnorm(0, 1e-10)
  for (p in 1:P){
    beta[p] ~ dnorm(0, eta[idx[p]])
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
jagsdata = list("y" = y, "X" = X, "N" = nrow(X), "P" = ncol(X), "idx" = idx, "nG" = max(idx))
inits = lapply(1:chains, function(z) list("beta" = rep(0, ncol(X)), 
                                          "Intercept" = 0, 
                                          "group_lambda" =  rep(1, max(idx)),
                                          "global_lambda"= 1, 
                                          "tau" = 1,
                                          .RNG.name= "lecuyer::RngStream",
                                          .RNG.seed= sample(1:10000, 1)))
monitor = c("Intercept", "beta", "sigma", "Deviance" , "global_lambda", "group_lambda", "ySim", "log_lik")
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
