#' Adaptive Bayesian Lasso
#'
#' @description The Bayesian LASSO of Leng, Tran and David Nott (2018). Basically just the Bayesian Lasso of Park & Casella (2008) but with
#' individual lambdas on each parameter defined by a gamma(sh, ra) distribution, where sh and ra are shape and rate hyperparameters. 
#' Here sh and ra are given independent gamma(4, 8) and gamma(0.0016, 0.1600) priors respectively. This places the expected
#' values for the shape and rate parameters at 0.50 and 0.01 respectively, which is consistent with the gamma(0.50, 0.01) prior on lambda
#' used for most other shrinkage models in this package. For alternatives 
#' see \code{\link[Bayezilla]{negLASSO}} (which is extremely similar) or \code{\link[Bayezilla]{extLASSO}}.
#' \cr
#' \cr 
#' Model Specification:
#' \cr
#' \if{html}{\figure{adaLASSO.png}{}}
#' \if{latex}{\figure{adaLASSO.png}{}}
#'
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
#' \cr
#' Leng, C., Tran, M.N., & Nott, D.J. (2014). Bayesian adaptive Lasso. arXiv:1009.2300 \cr
#' @return
#' a runjags object
#' @export
#'
#' @seealso 
#' \code{\link[Bayezilla]{negLASSO}} 
#' \code{\link[Bayezilla]{extLASSO}}
#' \code{\link[Bayezilla]{blasso}}
#' \code{\link[Bayezilla]{HS}}
#' \code{\link[Bayezilla]{HSplus}}
#' \code{\link[Bayezilla]{HSreg}}
#'
#' @examples
#' adaLASSO()
#'
adaLASSO = function(formula, data, log_lik = FALSE, iter=10000, warmup=1000, adapt=2000, chains=4, thin=1, method = "parallel", cl = makeCluster(2), ...){

  X = model.matrix(formula, data)[,-1]
  y = model.frame(formula, data)[,1]

  jags_blasso = "model{
  tau ~ dgamma(.01, .01)
  sh ~ dgamma(4, 8)
  ra ~ dgamma(0.0016, 0.1600)
    
  for (p in 1:P){
    lambda[p] ~ dgamma(sh , ra)
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
  monitor <- c("Intercept", "beta", "sigma", "sh" , "ra", "lambda", "Deviance", "ySim", "log_lik")
  if (log_lik == FALSE){
    monitor = monitor[-(length(monitor))]
  }
  inits <- lapply(1:chains, function(z) list("Intercept" = 0, .RNG.name= "lecuyer::RngStream", .RNG.seed= sample(1:10000, 1), "beta" = rep(0, P), "eta" = rep(1, P), "sh" = .5, "ra" = .05, "lambda" = sample(1:50, size = P, replace =TRUE), "tau" = 1))

  out = run.jags(model = "jags_blasso.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
  file.remove("jags_blasso.txt")
  if (!is.null(cl)) {
    parallel::stopCluster(cl = cl)
  }
  return(out)
}
