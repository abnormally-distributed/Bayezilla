#' Regularized Horseshoe
#'
#' @description This is the horseshoe model described by Piironen & Vehtari (2017). This tends to run very quickly
#' even for larger data sets or larger numbers of predictors and in my experience is faster and more stable (at least
#' on the tested data sets!) than the same model implemetned in Stan. If the horseshoe+ is analagous to the 
#' extended Bayesian LASSO, then this could be compared to the Bayesian Elastic Net.
#' 
#'
#' @references
#' Piironen, Juho; Vehtari, Aki. Sparsity information and regularization in the horseshoe and other shrinkage priors. Electron. J. Statist. 11 (2017), no. 2, 5018--5051. doi:10.1214/17-EJS1337SI. https://projecteuclid.org/euclid.ejs/1513306866 \cr
#' 
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
#' regHS()
#'
regHS = function(formula, data, phi = .10, slab_scale = 2, slab_df = 6, log_lik = FALSE, 
                 iter = 10000, warmup = 4000, adapt = 5000, chains=4, thin=1, method = "parallel", cl = makeCluster(2), ...){
  
  X = model.matrix(formula, data)[,-1]
  y = model.frame(formula, data)[,1]
  
  horseshoe = "model{

# tau is the precision, inverse of variance.
tau ~ dscaled.gamma(1, 3) 
sigma <- sqrt(1/tau)

# lambda squared, the global penalty
lambda2_0 <- pow((phi / (1-phi)) * (sigma/sqrt(N)), 2)
lambda2 ~ dt(0, 1 / lambda2_0, 1) T(0, )

# control parameter
c2_inv ~ dgamma(df / 2, (df*prior_variance) / 2 ) 
c2 <- 1 / c2_inv

# Coefficients
Intercept ~ dnorm(0, 1)
for (i in 1:P){
  eta[i] ~ dt(0, 1 / lambda2, 1) T(0, ) # prior variance
  eta_tilde[i] <- (eta[i]*c2) / (c2 + (lambda2 * eta[i]))
  beta[i] ~ dnorm(0, 1 / (eta_tilde[i] * lambda2))
  delta[i] <- 1 - ( 1 / (1+eta_tilde[i]) ) 
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
                                          "lambda2"= 1, 
                                          "tau" = 1,
                                          .RNG.name= "lecuyer::RngStream",
                                          .RNG.seed= sample(1:10000, 1)))
monitor = c("Intercept", "beta", "sigma", "Deviance" , "lambda2", "c2", "delta", "eta_tilde", "eta", "ySim", "log_lik")
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
