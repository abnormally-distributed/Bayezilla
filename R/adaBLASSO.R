#' Adaptive Bayesian Lasso
#'
#' The Bayesian LASSO of Leng, Tran and David Nott (2018). Basically just the Bayesian Lasso of Park & Casella (2008) but with
#' individual lambdas on each parameter defined by a gamma(r, d) distribution, where r and d are hyperparameters. Here r and d
#' are given independent gamma(0.0001, 0.0001) priors which approximates a Jeffrey's prior. For alternatives that may peform better
#' see \code{\link[Bayezilla]{negLASSO}} or \code{\link[Bayezilla]{extLASSO}}. However this is provided for a rich choice of options,
#' as it is sometimes hard to tell a priori which LASSO variant will work the best.
#'
#' @param formula the model formula
#' @param data a data frame.
#' @param log_lik Should the log likelihood be monitored? The default is FALSE.
#' @param iter How many post-warmup samples? Defaults to 10000.
#' @param warmup How many warmup samples? Defaults to 1000.
#' @param adapt How many adaptation steps? Defaults to 2000.
#' @param chains How many chains? Defaults to 4.
#' @param thin Thinning interval. Defaults to 1.
#' @param method Defaults to "parallel". For an alternative parallel option, choose "rjparallel" or. Otherwise, "rjags" (single core run).
#' @param cl Use parallel::makeCluster(# clusters) to specify clusters for the parallel methods. Defaults to two cores.
#' @param ... Other arguments to run.jags.
#'
#' @references
#' Park, T., & Casella, G. (2008). The Bayesian Lasso. Journal of the American Statistical Association, 103(482), 681-686. Retrieved from http://www.jstor.org/stable/27640090 \cr
#' Leng, C., Tran, M.N., & Nott, D.J. (2014). Bayesian adaptive Lasso. arXiv:1009.2300 \cr
#' @return
#' a runjags object
#' @export
#'
#'
#' @examples
#' adaBLASSO()
#'
adaBLASSO = function(formula, data, log_lik = FALSE, iter=10000, warmup=1000, adapt=2000, chains=4, thin=1, method = "parallel", cl = makeCluster(2), ...){

  X = model.matrix(formula, data)[,-1]
  y = model.frame(formula, data)[,1]

  jags_blasso = "model{
  tau ~ dgamma(.001, .001)
  r ~ rgamma(.0001, .0001)
  d ~ rgamma(.0001, .0001)
  for (p in 1:P){
    lambda[p] ~ dgamma(r, d)
    eta[p] ~ dexp(lambda[p]^2 / 2)
    omega[p] <- 1 / ( (1 / tau) * eta[p])
    beta[p] ~ dnorm(0, omega[p])
  }
  Intercept ~ dnorm(0, 1)
  for (i in 1:N){
    y[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
    log_lik[i] <- logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]), tau)
    ySim[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
  }
  sigma <- sqrt(1/tau)
  Deviance <- -2 * sum(log_lik[1:N])
}"
  P <- ncol(X)
  write_lines(jags_blasso, "jags_blasso.txt")
  jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X))
  monitor <- c("Intercept", "beta", "sigma", "r" , "d", "lambda", "Deviance", "eta", "ySim", "log_lik")
  if (log_lik == FALSE){
    monitor = monitor[-(length(monitor))]
  }
  inits <- lapply(1:chains, function(z) list("Intercept" = 0, .RNG.name= "lecuyer::RngStream", .RNG.seed= sample(1:10000, 1), "beta" = rep(0, P), "eta" = rep(1, P), "lambda" = rep(1, P), "tau" = 1, "r" = .001, "d" = .001))

  out = run.jags(model = "jags_blasso.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
  file.remove("jags_blasso.txt")
  if (!is.null(cl)) {
    parallel::stopCluster(cl = cl)
  }
  return(out)
}
