#' Group Bayesian Lasso 
#'
#' Group selection was introduced in the group LASSO by Yuan and Lin (2006) in
#' the context of the classical "frequentist" LASSO. The concept is adapted here to the Bayesian LASSO 
#' following the example of Kyung et al. (2010)\cr
#' \cr
#' Note that for the binomial and poisson likelihood functions 
#' the New Bayesian LASSO is adapted for use here, which utilizes a scale mixture of
#' uniform distributions to obtain the Laplacian priors (Mallick & Yi, 2014). I have found that this parameterization
#' simply samples faster for the binomial and poisson models, but is logically equivalent to the normal-exponential
#' mixture parameterization. Plug-in pseudovariances are used for these. 
#' 
#'
#' \cr
#' Model Specification: \cr
#' \cr
#' \if{html}{\figure{groupBLASSO.png}{}}
#' \if{latex}{\figure{groupBLASSO.png}{}}
#' \cr
#' \cr
#' Plugin Pseudo-Variances: \cr
#' \cr
#' \if{html}{\figure{pseudovar.png}{}}
#' \if{latex}{\figure{pseudovar.png}{}}
#'
#' @param X the model matrix. Construct this manually with model.matrix()[,-1]
#' @param y the outcome variable
#' @param idx the group labels. Should be of length = to ncol(model.matrix()[,-1]) with the group assignments
#' for each covariate. Please ensure that you start numbering with 1, and not 0.
#' @param family one of "gaussian", "binomial", or "poisson".
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
#' \cr
#' Kyung, M., Gill, J., Ghosh, M., and Casella, G. (2010). Penalized regression, standard errors, and bayesian lassos. Bayesian Analysis, 5(2):369–411. \cr
#' \cr
#' Mallick, H., & Yi, N. (2014). A New Bayesian Lasso. Statistics and its interface, 7(4), 571–582. doi:10.4310/SII.2014.v7.n4.a12 \cr
#' 
#' @return
#' a runjags object
#' 
#'
#' @examples
#' groupBLASSO()
#' 
#' @export
groupBLASSO = function(X, y, idx, family = "gaussian", log_lik = FALSE, iter=10000, warmup=1000, adapt=2000, chains=4, thin=1, method = "parallel", cl = makeCluster(2), ...){
  
if (family == "gaussian"){
    
  jags_group_blasso = "model{
  
  tau ~ dgamma(0.01, 0.01) 
  sigma2 <- 1/tau
  lambda ~ dgamma(0.50 , 0.20)
  
  # Group Level shrinkage
  
  for (g in 1:nG){
    eta[g] ~ dgamma( (k[g] + 1) * 0.50 , pow(lambda, 2) * 0.50)
    omega[g] <- 1 / (sigma2 * eta[g])
  }
  
  for (p in 1:P){
    beta[p] ~ dnorm(0, omega[idx[p]])
  }
  
  Intercept ~ dnorm(0, 1e-10)
  
  for (i in 1:N){
    y[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
    log_lik[i] <- logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]), tau)
    ySim[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
  }
  
  sigma <- sqrt(1/tau)
  Deviance <- -2 * sum(log_lik[1:N])
  
}"
  
  P <- ncol(X)
  write_lines(jags_group_blasso, "jags_group_blasso.txt")
  jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), idx = idx, nG = max(idx), k = as.vector(table(idx)))
  monitor <- c("Intercept", "beta", "sigma", "Deviance", "lambda", "eta", "ySim", "log_lik")
  if (log_lik == FALSE){
    monitor = monitor[-(length(monitor))]
  }
  
  inits <- lapply(1:chains, function(z) list("Intercept" = lmSolve(y ~ ., data = data.frame(y = y, X))[1],
                                             "beta" = lmSolve(y ~ ., data = data.frame(y = y, X))[-1], 
                                             "eta" = rep(2, max(idx)), 
                                             "lambda" = 2, 
                                             "tau" = 1, 
                                             "ySim" = sample(y, length(y)),
                                             .RNG.name= "lecuyer::RngStream", 
                                             .RNG.seed= sample(1:10000, 1)
                                             ))
  
  out = run.jags(model = "jags_group_blasso.txt", n.chains = chains, modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
  file.remove("jags_group_blasso.txt")
  if (!is.null(cl)) {
    parallel::stopCluster(cl = cl)
  }
  return(out)
}


if (family == "binomial"){
  
  jags_group_blasso = "model{
  
  lambda ~ dgamma(0.50 , 0.20)
  
  for (i in 1:nG){
    u[g] ~ dgamma( k[g] + 1  , lambda)
  }
  
  for (p in 1:P){  
    beta[p] ~ dunif(-1 * (sigma * u[idx[p]]), (sigma * u[idx[p]]) )
  }
  
  Intercept ~ dnorm(0, 1e-10)
  
    for (i in 1:N){
      logit(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
      y[i] ~ dbern(psi[i])
      log_lik[i] <- logdensity.bern(y[i], psi[i])
      ySim[i] ~ dbern(psi[i])
    }
  }
  Deviance <- -2 * sum(log_lik[1:N])
}"
  
  P <- ncol(X)
  nG = max(idx)
  write_lines(jags_group_blasso, "jags_group_blasso.txt")
  jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), sigma = sqrt(pow(mean(y), -1) * pow(1 - mean(y), -1)), idx = idx, nG = max(idx), k = as.vector(table(idx)))
  monitor <- c("Intercept", "beta", "lambda", "Deviance", "ySim", "log_lik")
  if (log_lik == FALSE){
    monitor = monitor[-(length(monitor))]
  }
  inits <- lapply(1:chains, function(z) list("Intercept" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "binomial", lambda = 0.025, alpha = 0, standardize = FALSE))[1,1]), 
                                             "beta" = rep(0, P), 
                                             "u" = rgamma(nG, nG + 1, 1), 
                                             "lambda" = 1, 
                                             "ySim" = sample(y, length(y)),
                                             .RNG.name= "lecuyer::RngStream", 
                                             .RNG.seed= sample(1:10000, 1)))
  
  out = run.jags(model = "jags_group_blasso.txt", n.chains = chains, modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
  file.remove("jags_group_blasso.txt")
  if (!is.null(cl)) {
    parallel::stopCluster(cl = cl)
  }
  return(out)
}


if (family == "poisson"){
  
  jags_group_blasso = "model{
  
  lambda ~ dgamma(0.50 , 0.20)
  
  for (i in 1:nG){
    u[g] ~ dgamma( k[g] + 1  , lambda)
  }
  
  for (p in 1:P){  
    beta[p] ~ dunif(-1 * (sigma * u[idx[p]]), (sigma * u[idx[p]]) )
  }
  
  Intercept ~ dnorm(0, 1e-10)
  
  for (i in 1:N){
    log(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
    y[i] ~ dpois(psi[i])
    log_lik[i] <- logdensity.pois(y[i], psi[i])
    ySim[i] ~ dpois(psi[i])
}
              
  Deviance <- -2 * sum(log_lik[1:N])
}"

P <- ncol(X)
nG = max(idx)
write_lines(jags_group_blasso, "jags_group_blasso.txt")
jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), sigma = sqrt(pow(mean(y) , -1)), idx = idx, nG = max(idx), k = as.vector(table(idx)))
monitor <- c("Intercept", "beta", "lambda", "Deviance", "ySim", "log_lik")
if (log_lik == FALSE){
  monitor = monitor[-(length(monitor))]
}
inits <- lapply(1:chains, function(z) list("Intercept" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "poisson", lambda = 0.025, alpha = 0, standardize = FALSE))[1,1]), 
                                           "beta" = rep(0, P), 
                                           "u" = rgamma(nG, nG+1, 1), 
                                           "lambda" = 1, 
                                           "ySim" = sample(y, length(y)),
                                           .RNG.name= "lecuyer::RngStream", 
                                           .RNG.seed= sample(1:10000, 1)))

out = run.jags(model = "jags_group_blasso.txt", n.chains = chains, modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
file.remove("jags_group_blasso.txt")
if (!is.null(cl)) {
  parallel::stopCluster(cl = cl)
}
return(out)
}
}

