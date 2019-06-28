#' Group Bayesian Lasso 
#'
#' Group selection was introduced in the group LASSO by Yuan and Lin (2006) in
#' the context of the classical "frequentist" LASSO. The concept is adapted here to the Bayesian LASSO of Park & Casella (2008). Note only the Gaussian likelihood
#' is provided because the Bayesian LASSO requires conditioning on the error variance, which GLM-families
#' do not have. \cr
#'
#' \cr
#' Model Specification:
#' \cr
#' \if{html}{\figure{groupBlasso.png}{}}
#' \if{latex}{\figure{groupBlasso.png}{}}
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
#' @param method Defaults to "parallel". For an alternative parallel option, choose "rjparallel" or. Otherwise, "rjags" (single core run).
#' @param cl Use parallel::makeCluster(# clusters) to specify clusters for the parallel methods. Defaults to two cores.
#' @param ... Other arguments to run.jags.
#'
#' @references 
#' \cr
#' Yuan, Ming; Lin, Yi (2006). Model Selection and Estimation in Regression with Grouped Variables. Journal of the Royal Statistical Society. Series B (statistical Methodology). Wiley. 68 (1): 49–67. doi:10.1111/j.1467-9868.2005.00532.x \cr
#' \cr
#' Park, T., & Casella, G. (2008). The Bayesian Lasso. Journal of the American Statistical Association, 103(482), 681-686. Retrieved from http://www.jstor.org/stable/27640090 \cr
#'
#' @return
#' a runjags object
#' 
#' @seealso 
#' \code{\link[Bayezilla]{blassoDC}} 
#' \code{\link[Bayezilla]{adaLASSO}}
#' \code{\link[Bayezilla]{negLASSO}} 
#' \code{\link[Bayezilla]{blasso}}
#' \code{\link[Bayezilla]{HS}}
#' \code{\link[Bayezilla]{HSplus}}
#' \code{\link[Bayezilla]{HSreg}}
#'
#' @examples
#' groupBlasso()
#' 
#' @export
groupBlasso = function(X, y, idx, log_lik = FALSE, iter=10000, warmup=1000, adapt=2000, chains=4, thin=1, method = "parallel", cl = makeCluster(2), ...){
  

  jags_blasso = "model{
  tau ~ dgamma(.01, .01) 
  sigma2 <- 1/tau
  lambda ~ dexp(0.002)

  # Group Level shrinkage
  for (g in 1:nG){
    eta[g] ~ dgamma( (m[g] + 1) * 0.50 , pow(lambda, 2) * 0.50)
    omega[g] <- 1 / (sigma2 * eta[g])
  }
  
  for (p in 1:P){
    beta[p] ~ dnorm(0, omega[idx[p]])
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
  jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), idx = idx, nG = max(idx), m = as.vector(table(idx)))
  monitor <- c("Intercept", "beta", "sigma", "Deviance", "lambda", "eta", "ySim", "log_lik")
  if (log_lik == FALSE){
    monitor = monitor[-(length(monitor))]
  }
  inits <- lapply(1:chains, function(z) list("Intercept" = 0,
                                             "beta" = rep(0, P), 
                                             "eta" = rep(2, max(idx)), 
                                             "lambda" = 2, 
                                             "tau" = 1, 
                                             .RNG.name= "lecuyer::RngStream", 
                                             .RNG.seed= sample(1:10000, 1)
                                             ))
  
  out = run.jags(model = "jags_blasso.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
  file.remove("jags_blasso.txt")
  if (!is.null(cl)) {
    parallel::stopCluster(cl = cl)
  }
  return(out)
}

