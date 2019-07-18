#' Bayesian Group Shape Adaptive Shrinkage Prior Regression
#'
#' @description The Group Bayesian Bridge model of Mallick & Yi (2018) adapted to the Shape Adaptive Shrinkage Prior (SASP) 
#' of Sillanpää &  Mutshinda (2011).
#' 
#' Bridge regression allows you to utilize different Lp norms for the shape  of the prior through the shape parameter kappa of 
#' the power exponential distribution (also known as generalized Gaussian). Norms of 1 and 2 give the Laplace and Gaussian 
#' distributions respectively (corresponding to the LASSO and Ridge Regression). Norms smaller than 1 are very difficult to 
#' estimate directly, but have very tall modes at zero and very long, cauchy like tails. Values greater than 2 become increasingly 
#' platykurtic, with the uniform distribution arising as it approaches infinity. \cr 
#' \cr
#' The benefit of the shape adaptive shrinkage prior is that one need not pick a specific norm. Hence, if there is uncertainty over
#' whether or not one wishes to choose the L1 norm (LASSO) or L2 norm (Ridge), this integrates over a reasonable range of values. The gamma
#' prior for the norm has an expected value of 1.4, which gives a reasonable compromise between the LASSO and Ridge. \cr
#' \cr
#' JAGS has no built in power exponential distribution, so the distribution is parameterized as a uniform-gamma mixture just as in Mallick & Yi (2018). 
#' The parameterization is given below. For generalized linear models plug-in pseudovariances are used. \cr
#' \cr
#' Model Specification:
#' \cr
#' \if{html}{\figure{groupSasp.png}{}}
#' \if{latex}{\figure{groupSasp.png}{}}
#'
#' \cr
#' Plugin Pseudo-Variances: \cr
#' \if{html}{\figure{pseudovar.png}{}}
#' \if{latex}{\figure{pseudovar.png}{}}
#'
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
#' Kyung, M., Gill, J., Ghosh, M., and Casella, G. (2010). Penalized regression, standard errors, and bayesian lassos. Bayesian Analysis, 5(2):369–411.\cr
#' \cr
#' Mallick, H. & Yi, N. (2018). Bayesian bridge regression, Journal of Applied Statistics, 45:6, 988-1008, DOI: 10.1080/02664763.2017.1324565 \cr
#' \cr
#' Mallick, H., & Yi, N. (2014). A New Bayesian Lasso. Statistics and its interface, 7(4), 571–582. doi:10.4310/SII.2014.v7.n4.a12 \cr
#' \cr
#' Sillanpää, S., &  Mutshinda, C., (2011). Bayesian shrinkage analysis of QTLs under shape-adaptive shrinkage priors, and accurate re-estimation of genetic effects. Heredity volume 107, pages 405–412. doi: 10.1038/hdy.2011.37 \cr
#' \cr
#' 
#' @return
#' a runjags object
#' @export
#'
#' @examples
#' groupSasp()
#' 
groupSasp = function(X, y, idx, family = "gaussian", log_lik = FALSE, iter=10000, warmup=1000, adapt=2000, chains=4, thin=1, method = "parallel", cl = makeCluster(2), ...){
  
  if (family == "gaussian"){
    
    jags_bridge = "model{
  
  kappa ~ dgamma(1.96, 1.40)
  tau ~ dgamma(.01, .01) 
  sigma <- sqrt(1/tau)
  
  for (g in 1:nG){
    lambda[g] ~ dgamma(0.50, 0.20)
    u[g] ~ dgamma( (k[g]/kappa) + 1  , lambda[g])
  }
  
  for (p in 1:P){  
    beta[p] ~ dunif(-1 * pow(sigma * u[idx[p]], 1/kappa), pow(sigma * u[idx[p]], 1/kappa))
  }
  
  Intercept ~ dnorm(0, 1e-10)
  
  for (i in 1:N){
    y[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
    log_lik[i] <- logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]), tau)
    ySim[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
  }
  
  Deviance <- -2 * sum(log_lik[1:N])
}"
    
    P <- ncol(X)
    nG = max(idx)
    write_lines(jags_bridge, "jags_bridge.txt")
    jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), idx = idx, nG = max(idx), k = as.vector(table(idx)))
    monitor <- c("Intercept", "beta", "sigma", "lambda", "kappa", "Deviance", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits <- lapply(1:chains, function(z) list("Intercept" = lmSolve(y ~ ., data = data.frame(y = y, X))[1], 
                                               "beta" = rep(0, P), 
                                               "u" = rgamma(nG, (1 / .75) + 1, 1), 
                                               "lambda" = rep(1, nG), 
                                               "tau" = 1, 
                                               "kappa" = .75, 
                                               "ySim" = y, 
                                               .RNG.name= "lecuyer::RngStream", 
                                               .RNG.seed= sample(1:10000, 1)))
    
    out = run.jags(model = "jags_bridge.txt", n.chains = chains, modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
    file.remove("jags_bridge.txt")
    if (!is.null(cl)) {
      parallel::stopCluster(cl = cl)
    }
    return(out)
  }
  
  if (family == "binomial" || family == "logistic"){
    
    jags_bridge = "model{


  kappa ~ dgamma(1.96, 1.40)

  for (g in 1:nG){
    lambda[g] ~ dgamma(0.50, 0.20)
    u[g] ~ dgamma( (k[g]/kappa) + 1  , lambda[g])
  }
  
  for (p in 1:P){  
    beta[p] ~ dunif(-1 * pow(sigma * u[idx[p]], 1/kappa), pow(sigma * u[idx[p]], 1/kappa))
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
    write_lines(jags_bridge, "jags_bridge.txt")
    jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), sigma = sqrt(pow(mean(y), -1) * pow(1 - mean(y), -1)), idx = idx, nG = max(idx), k = as.vector(table(idx)))
    monitor <- c("Intercept", "beta", "lambda", "kappa", "Deviance", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits <- lapply(1:chains, function(z) list("Intercept" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "binomial", lambda = 0.025, alpha = 0, standardize = FALSE))[1,1]), 
                                               "beta" = rep(0, P), 
                                               "u" = rgamma(nG, (1 / .75) + 1, 1), 
                                               "lambda" = rep(1, nG), 
                                               "kappa" = .75, 
                                               "ySim" = y, 
                                               .RNG.name= "lecuyer::RngStream", 
                                               .RNG.seed= sample(1:10000, 1)))
    
    out = run.jags(model = "jags_bridge.txt", n.chains = chains, modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
    file.remove("jags_bridge.txt")
    if (!is.null(cl)) {
      parallel::stopCluster(cl = cl)
    }
    return(out)
  }
  
  
  if (family == "poisson"){
    
    jags_bridge = "model{

  kappa ~ dgamma(1.96, 1.40)

  for (g in 1:nG){
    lambda[g] ~ dgamma(0.50, 0.20)
    u[g] ~ dgamma( (k[g]/kappa) + 1  , lambda[g])
  }
  
  for (p in 1:P){  
    beta[p] ~ dunif(-1 * pow(sigma * u[idx[p]], 1/kappa), pow(sigma * u[idx[p]], 1/kappa))
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
  write_lines(jags_bridge, "jags_bridge.txt")
  jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), sigma = sqrt(pow(mean(y) , -1)), idx = idx, nG = max(idx), k = as.vector(table(idx)))
  monitor <- c("Intercept", "beta", "lambda", "kappa", "Deviance", "ySim", "log_lik")
  if (log_lik == FALSE){
    monitor = monitor[-(length(monitor))]
  }
  inits <- lapply(1:chains, function(z) list("Intercept" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "poisson", lambda = 0.025, alpha = 0, standardize = FALSE))[1,1]), 
                                             "beta" = rep(0, P), 
                                             "u" = rgamma(nG, (1 / .75) + 1, 1), 
                                             "lambda" = rep(1, nG), 
                                             "ySim" = y, 
                                             "kappa" = .75, 
                                             .RNG.name= "lecuyer::RngStream", 
                                             .RNG.seed= sample(1:10000, 1)))
  
  out = run.jags(model = "jags_bridge.txt", n.chains = chains, modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
  file.remove("jags_bridge.txt")
  if (!is.null(cl)) {
    parallel::stopCluster(cl = cl)
  }
  return(out)
  }

}

