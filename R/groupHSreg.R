#' Group Regularized Horseshoe
#'
#' @description This is the horseshoe model described by Piironen & Vehtari (2017) 
#' adapted for group selection, in the spirit of the group LASSO.
#' \cr
#' \cr
#' Model Specification: \cr 
#' \cr
#' \if{html}{\figure{groupregularizedHorseshoe.png}{}}
#' \if{latex}{\figure{groupregularizedHorseshoe.png}{}}
#'
#' \cr
#' @references
#' Piironen, Juho; Vehtari, Aki. Sparsity information and regularization in the horseshoe and other shrinkage priors. Electron. J. Statist. 11 (2017), no. 2, 5018--5051. doi:10.1214/17-EJS1337SI. https://projecteuclid.org/euclid.ejs/1513306866 \cr
#' \cr
#' Carvalho, C. M., Polson, N. G., and Scott, J. G. (2010). The horseshoe estimator for sparse signals. Biometrika, 97(2):465–480. \cr
#'
#' @param formula the model formula
#' @param data a data frame.
#' @param phi your prior guess on the inclusion probability. Defaults to .50. Best way to come up with a figure is a prior guess on how many groups are non-zero out of the total number of groups.
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
#' groupHSreg()
#'
groupHSreg = function(X, y, idx, phi = .50, slab_scale = 2, slab_df = 3, log_lik = FALSE, 
                 iter = 10000, warmup = 4000, adapt = 5000, chains=4, thin=1, method = "parallel", cl = makeCluster(2), ...){
  
horseshoe = "model{

# tau is the precision, inverse of variance.
tau ~ dgamma(.01, .01) 
sigma <- sqrt(1/tau)

# lambda squared, the global penalty
lambda0sqrd <- pow((phi / (1-phi)) * (sigma/sqrt(N)), 2)
lambda ~ dt(0, 1/lambda0sqrd, 1) T(0, )

# Group Penalties
for (g in 1:nG){
  # control parameter
  c2_inv[g] ~ dgamma(df / 2, (df*prior_variance)/2) 
  c2[g] <- 1 / c2_inv
  eta[g] ~ dt(0,1,1) T(0, ) 
  eta_tilde[g] <- (pow(eta[g], 2)*c2[g]) / (c2[g]+(pow(lambda,2)*pow(eta[g], 2)))
  eta_inv[g] <- 1 / (eta_tilde[g] * lambda)
}

# Coefficients
Intercept ~ dnorm(0, 1)
for (i in 1:P){
  beta[i] ~ dnorm(0, eta_inv[idx[i]])
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
jagsdata = list("y" = y, 
                "X" = X, 
                "N" = nrow(X), 
                "idx" = idx, 
                "nG" = max(idx),
                "P" = ncol(X), 
                "phi" = phi, 
                "df" = slab_df, 
                "prior_variance" = square(slab_scale))
inits = lapply(1:chains, function(z) list("beta" = rep(0, ncol(X)), 
                                          "Intercept" = 0, 
                                          "eta" =  rep(1, max(idx)),
                                          "c2_inv" = rep(1/slab_df, max(idx)),
                                          "lambda"= 10, 
                                          "tau" = 1,
                                          .RNG.name= "lecuyer::RngStream",
                                          .RNG.seed= sample(1:10000, 1)))
monitor = c("Intercept", "beta", "sigma", "Deviance" , "lambda", "c2", "eta_tilde", "eta", "ySim", "log_lik")
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
