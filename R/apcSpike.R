#' Spike-Slab Adaptive Powered Correlation Prior for Gaussian likelihoods
#'
#' @param formula the model formula
#' @param data a data frame
#' @param lambda the power to use in the adaptive correlation prior. Default is -1, which gives the Zellner-Siow g prior. Setting
#' lambda to 0 results in a ridge-regression like prior. Setting lambda to a positive value adapts to collinearity by shrinking 
#' collinear variables towards a common value. Negative values of lambda pushes collinear variables further apart. I suggest 
#' fitting multiple values of lambda and selecting which is best via LOO-IC or WAIC. 
#' @param family one of "gaussian", "binomial", or "poisson".
#' @param log_lik Should the log likelihood be monitored? The default is FALSE.
#' @param iter How many post-warmup samples? Defaults to 10000.
#' @param warmup How many warmup samples? Defaults to 1000.
#' @param adapt How many adaptation steps? Defaults to 2000.
#' @param chains How many chains? Defaults to 4. 
#' @param thin Thinning interval. Defaults to 1.
#' @param method Defaults to "parallel". For an alternative parallel option, choose "rjparallel". Otherwise, "rjags" (single core run).
#' @param cl Use parallel::makeCluster(# clusters) to specify clusters for the parallel methods. Defaults to two cores.
#' @param ... Other arguments to run.jags.
#'
#' @return
#' A run.jags object.
#' @export
#'
#' @examples
#' apcLm()
apcSpike = function(formula, data, family = "gaussian", lambda = -1, log_lik = FALSE, iter=10000, warmup=1000, adapt=5000, chains=4, thin=1, method = "parallel", cl = makeCluster(2), ...)
{
  data <- as.data.frame(data)
  y <- as.numeric(model.frame(formula, data)[, 1])
  X <- model.matrix(formula, data)[, -1]
  ## Overkill to ensure that the correlation matrix is positive definite.
  cormat = cov2cor(pseudoinverse(pseudoinverse(cor(X, method = "spearman"))))
  L = eigen(cormat)$vectors
  D = eigen(cormat)$values
  Trace = function(mat){sum(diag(mat))}
  P = ncol(X)
  Dpower = rep(0, P)
  t = XtXinv(X, tol=1e-2)
  for(i in 1:P) {
    Dpower[i] <- (D[i]^lambda);
  }
  prior_cov = (L %*% diag(Dpower) %*% t(L)) / length(y)
  K = Trace(t) / Trace(prior_cov)
  prior_cov = K * (prior_cov)
  pd <- isPositiveDefinite(prior_cov)
  if (pd == FALSE){
    prior_cov <- solve(pseudoinverse(prior_cov))
  }
  
  if (family == "gaussian"){
    
    jags_apc = "model{
              
              phi ~ dbeta(1, 2)
              tau ~ dscaled.gamma(1, 3)
              g_inv ~ dgamma(.5, N * .5)
              g <- 1 / g_inv
              sigma <- sqrt(1/tau)
              
              for (j in 1:P){
                for (k in 1:P){
                  cov[j,k] = g * pow(sigma, 2) * prior_cov[j,k]
                }
              }
              
              omega <- inverse(cov) 
              theta[1:P] ~ dmnorm(rep(0,P), omega[1:P,1:P])
              for (i in 1:P){
                delta[i] ~ dbern(phi)
                beta[i] <- delta[i] * theta[i]
              }
              Intercept ~ dnorm(0, .01)
              for (i in 1:N){
                 y[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 log_lik[i] <- logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 ySim[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
              }
              Deviance <- -2 * sum(log_lik[1:N])
          }"
    
    P = ncol(X)
    write_lines(jags_apc, "jags_apc.txt")
    jagsdata = list(X = X, y = y,  N = length(y), P = ncol(X), prior_cov = prior_cov)
    monitor = c("Intercept", "beta", "sigma", "g", "Deviance", "phi", "delta", "ySim" ,"log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" =0, "phi" = .2 , "delta" = rep(0, P), "theta" = jitter(rep(0, P), amount = 1), "tau" = 1, "g_inv" = 1/length(y), "ySim" = y, .RNG.name= "lecuyer::RngStream", .RNG.seed = sample(1:10000, 1)))
  }
  
  if (family == "binomial" || family == "logistic"){
    
    jags_apc = "model{
    
              phi ~ dbeta(1, 2)
              g_inv ~ dgamma(.5, N * .5)
              g <- 1 / g_inv
              
              for (j in 1:P){
                for (k in 1:P){
                  cov[j,k] = g * 1.0 * prior_cov[j,k]
                }
              }
              
           
              omega <- inverse(cov) 
              theta[1:P] ~ dmnorm(rep(0,P), omega[1:P,1:P])
              for (i in 1:P){
                delta[i] ~ dbern(phi)
                beta[i] <- delta[i] * theta[i]
              }
              for (i in 1:N){
                 logit(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                 y[i] ~ dbern(psi[i])
                 log_lik[i] <- logdensity.bern(y[i], psi[i])
                 ySim[i] ~ dbern(psi[i])
              }
             Deviance <- -2 * sum(log_lik[1:N])
          }"
    
    P = ncol(X)
    write_lines(jags_apc, "jags_apc.txt")
    jagsdata = list(X = X, y = y, N = length(y), P = ncol(X), prior_cov = prior_cov)
    monitor = c("Intercept", "beta", "g", "Deviance", "phi", "delta","ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = 0, "phi" = .2 , "delta" = rep(0, P), "theta" = jitter(rep(0, P), amount = 1), "g_inv" = 1/length(y), "ySim" = y, .RNG.name= "lecuyer::RngStream", .RNG.seed= sample(1:10000, 1)))
  }
  
  if (family == "poisson"){
    
    jags_apc = "model{
    
              phi ~ dbeta(1, 2)

              g_inv ~ dgamma(.5, N * .5)
              g <- 1 / g_inv
              
              for (j in 1:P){
                for (k in 1:P){
                  cov[j,k] = g * 1.0 * prior_cov[j,k]
                }
              }
              
              omega <- inverse(cov) 
              theta[1:P] ~ dmnorm(rep(0,P), omega[1:P,1:P])
              for (i in 1:P){
                delta[i] ~ dbern(phi)
                beta[i] <- delta[i] * theta[i]
              }
              for (i in 1:N){
                 log(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                 y[i] ~ dpois(psi[i])
                 log_lik[i] <- logdensity.pois(y[i], psi[i])
                 ySim[i] ~ dpois(psi[i])
              }
              Deviance <- -2 * sum(log_lik[1:N])
          }"
    
    write_lines(jags_apc, "jags_apc.txt")
    P = ncol(X)
    jagsdata = list(X = X, y = y, N = length(y),  P = ncol(X), prior_cov = prior_cov)
    monitor = c("Intercept", "beta", "g", "Deviance", "phi", "delta", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = 0, "g_inv" = 1/length(y), "ySim" = y, .RNG.name= "lecuyer::RngStream", .RNG.seed= sample(1:10000, 1),"phi" = .2 , "delta" = rep(0, P), "theta" = jitter(rep(0, P), amount = 1)))
  }
  
  out = run.jags(model = "jags_apc.txt", modules = c("glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE,...)
  if (is.null(cl) == FALSE){
    parallel::stopCluster(cl = cl)
  }
  file.remove("jags_apc.txt")
  return(out)
}
