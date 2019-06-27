#' Bayesian Lasso with unpenalized design covariates
#'
#' The Bayesian LASSO of Park & Casella (2008), but with the allowance for a set of covariates that are not penalized. 
#' For example, you may wish to include variables such as age and gender in all models so that 
#' the coefficients for the other variables are penalized while controlling for these. This
#' is a common need in research. \cr
#' \cr
#' Note only the Gaussian likelihood is provided because the Bayesian LASSO requires conditioning
#' on the error variance, which GLM-families do not have. \cr
#'
#' @param formula the model formula.
#' @param design.formula formula for the design covariates.
#' @param data a data frame.
#' @param log_lik Should the log likelihood be monitored? The default is FALSE.
#' @param iter How many post-warmup samples? Defaults to 10000.
#' @param warmup How many warmup samples? Defaults to 1000.
#' @param adapt How many adaptation steps? Defaults to 2000.
#' @param chains How many chains? Defaults to 4.
#' @param thin Thinning interval. Defaults to 1.
#' @param method Defaults to "parallel". For an alternative parallel option, choose "rjparallel" or.
#' Otherwise, "rjags" (single core run).
#' @param cl Use parallel::makeCluster(# clusters) to specify clusters for the parallel methods. 
#' Defaults to two cores.
#' @param ... Other arguments to run.jags.
#'
#' @references Park, T., & Casella, G. (2008). The Bayesian Lasso. Journal of the American Statistical Association, 
#' 103(482), 681-686. Retrieved from http://www.jstor.org/stable/27640090 \cr
#'
#' @return
#' a runjags object
#' @export
#' 
#'
#' @examples
#' blassoDC()
blassoDC = function(formula, design.formula, data, log_lik = FALSE, iter=10000, warmup=1000, adapt=2000, chains=4, thin=1, method = "parallel", cl = makeCluster(2), ...){
  
  X = model.matrix(formula, data)[,-1]
  y = model.frame(formula, data)[,1]
  FX <- model.matrix(design.formula, data)[, -1]
  
  jags_blasso = "model{
  tau ~ dgamma(0.01, 0.01)
  sigma2 <- 1/tau
  lambda ~ dgamma(0.5 , 0.001)
  
  for (p in 1:P){
    eta[p] ~ dexp(lambda^2 / 2)
    omega[p] <- 1 / (sigma2 * eta[p])
    beta[p] ~ dnorm(0, omega[p])
  }
  
  for (f in 1:FP){
    design_beta[f] ~ dnorm(0, 1e-200)
  }
  
  Intercept ~ dnorm(0, 1)
  
  for (i in 1:N){
    y[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]) + sum(design_beta[1:FP] * FX[i,1:FP]) , tau) 
    log_lik[i] <- logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]) + sum(design_beta[1:FP] * FX[i,1:FP]), tau)
    ySim[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]) + sum(design_beta[1:FP] * FX[i,1:FP]) , tau)
  }
  sigma <- sqrt(1/tau)
  Deviance <- -2 * sum(log_lik[1:N])
}"
  
  P <- ncol(X)
  FP <- ncol(FX)
  write_lines(jags_blasso, "jags_blasso.txt")
  jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), FP = FP, FX = FX)
  monitor <- c("Intercept", "beta", "design_beta", "sigma", "lambda", "Deviance", "eta", "ySim", "log_lik")
  if (log_lik == FALSE){
    monitor = monitor[-(length(monitor))]
  }
  inits <- lapply(1:chains, function(z) list("Intercept" = 0, 
                                             "beta" = rep(0, P), 
                                             "design_beta" = rep(0, FP), 
                                             "eta" = rep(1, P), 
                                             "lambda" = 2, 
                                             "tau" = 1, 
                                             .RNG.name= "lecuyer::RngStream", 
                                             .RNG.seed= sample(1:10000, 1)))
  
  out = run.jags(model = "jags_blasso.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
  file.remove("jags_blasso.txt")
  if (!is.null(cl)) {
    parallel::stopCluster(cl = cl)
  }
  return(out)
}

