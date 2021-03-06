#' Group Normal-exponential-gamma Bayesian LASSO
#'
#' @description This implements the normal-exponential-gamma "hyperlasso" of Griffin & Brown (2011)
#' adapted here to the Kyung et al.'s (2010) group Bayesian LASSO. 
#' 
#' This model has normal priors
#' on each coefficient, whose precision is modeled by group specific exponential distributions. 
#' The exponential
#' distributions in turn have their respective rate parameters modeled through 
#' independent group-level gamma(k_g * 0.50, 1 / lambda^2) distributions, where k_g is the number of predictors
#' in group g. If there is no grouping then this reduces to the \code{\link[Bayezilla]{negLASSO}}.
#' Lambda is a single top-level hyperparameter here given a gamma(0.50 , 0.20) prior. \cr 
#' \cr
#' The model specification is given below: \cr
#' \cr
#' \cr
#' \cr 
#' Model Specification:
#' \cr
#' \if{html}{\figure{groupNegLASSO.png}{}}
#' \if{latex}{\figure{groupNegLASSO.png}{}}
#' \cr
#' \cr
#' @references 
#' \cr
#' Kyung, M., Gill, J., Ghosh, M., and Casella, G. (2010). Penalized regression, standard errors, and Bayesian lassos. Bayesian Analysis, 5(2):369–411. \cr
#' \cr
#' Griffin, J. E. and Brown, P. J. (2011), Bayesian Hyper‐LASSOs With Non-Convex Penalization. Australian & New Zealand Journal of Statistics, 53: 423-442. doi:10.1111/j.1467-842X.2011.00641.x \cr
#' \cr
#' Yuan, Ming; Lin, Yi (2006). Model Selection and Estimation in Regression with Grouped Variables. Journal of the Royal Statistical Society. Series B (statistical Methodology). Wiley. 68 (1): 49–67. doi:10.1111/j.1467-9868.2005.00532.x \cr
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
#' @param method Defaults to "parallel". For another parallel option, choose "rjparallel" or "rjags" for a single core run.
#' @param cl Use parallel::makeCluster(# clusters) to specify clusters for the parallel methods. Defaults to 3.
#' @param ... Other arguments to run.jags.
#'
#' @return A run.jags object
#'
#' @examples
#' groupNegLASSO()
#'
#'
#' @export
groupNegLASSO = function(X, y, idx, family = "gaussian", log_lik = FALSE, iter=10000, warmup=1000, adapt=2000, chains=4, thin=1, method = "parallel", cl = makeCluster(3), ...){
  
  if (family == "gaussian" || family == "normal") {
    
    jags_neg_LASSO = "model{

              tau ~ dgamma(.01, .01)

              lambda ~ dgamma(0.5 , 0.20)

              Intercept ~ dnorm(0, 1e-10)

              for (g in 1:nG){
               eta[g] ~ dgamma(k[g] * 0.50 , 1 / pow(lambda, 2))
               psi[g] ~ dexp(eta[g])
              }
              
              for (p in 1:P){
                beta[p] ~ dnorm(0, 1 / psi[idx[p]])
              }
              
              for (i in 1:N){
                 y[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 log_lik[i] <- logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 ySim[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
              }
              
              sigma <- sqrt(1/tau)
              Deviance <- -2 * sum(log_lik[1:N])
          }"
    
    P <- ncol(X)
    write_lines(jags_neg_LASSO, "jags_neg_LASSO.txt")
    jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), idx = idx, nG = max(idx), k = as.vector(table(idx)))
    monitor <- c("Intercept", "beta", "sigma", "lambda", "Deviance", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits <- lapply(1:chains, function(z) list("Intercept" = lmSolve(y ~ ., data = data.frame(y = y, X))[1], 
                                               "beta" = lmSolve(y ~ ., data = data.frame(y = y, X))[-1], 
                                               "eta" = abs(jitter(rep(1, max(idx)), amount = 1)), 
                                               "psi" = abs(jitter(rep(1, max(idx)), amount = 1)), 
                                               "lambda" = 2, 
                                               "tau" = 1, 
                                               "ySim" = sample(y, length(y)),
                                               .RNG.name= "lecuyer::RngStream", 
                                               .RNG.seed = sample(1:10000, 1)))
  }
  
  
  if (family == "binomial" || family == "logistic") {
    
    jags_neg_LASSO = "model{

              lambda ~ dgamma(0.5 , 0.20)

              Intercept ~ dnorm(0, 1e-10)

              for (g in 1:nG){
               eta[g] ~ dgamma( (k[g] + 1) * 0.50 , pow(lambda, 2) * 0.50)
               psi[g] ~ dexp(eta[g])
              }
              for (p in 1:P){
                beta[p] ~ dnorm(0, 1 / psi[idx[p]])
              }
              
              for (i in 1:N){
                 logit(phi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                 y[i] ~ dbern(phi[i])
                 log_lik[i] <- logdensity.bern(y[i], phi[i])
                 ySim[i] ~ dbern(phi[i])
              }
   
              Deviance <- -2 * sum(log_lik[1:N])
          }"
    
    P = ncol(X)
    write_lines(jags_neg_LASSO, "jags_neg_LASSO.txt")
    jagsdata = list(X = X, y = y, N = length(y), P = ncol(X), idx = idx, nG = max(idx), k = as.vector(table(idx)))
    monitor = c("Intercept", "beta", "lambda", "Deviance", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = 0, 
                                              "beta" = rep(0, P), 
                                              "eta" = abs(jitter(rep(1, max(idx)), amount = 1)), 
                                              "psi" = abs(jitter(rep(1, max(idx)), amount = 1)), 
                                              "lambda" = 2, 
                                              "ySim" = sample(y, length(y)),
                                              .RNG.name= "lecuyer::RngStream", 
                                              .RNG.seed = sample(1:10000, 1)))
  }
  
  if (family == "poisson") {
    
    jags_neg_LASSO = "model{

              lambda ~ dgamma(0.5 , 0.20)

              Intercept ~ dnorm(0, 1e-10)

              for (g in 1:nG){
               eta[g] ~ dgamma( (k[g] + 1) * 0.50 , pow(lambda, 2) * 0.50)
               psi[g] ~ dexp(eta[g])
              }
              
              for (p in 1:P){
                beta[p] ~ dnorm(0, 1 / psi[idx[p]])
              }
              
              for (i in 1:N){
                 log(phi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                 y[i] ~ dpois(phi[i])
                 log_lik[i] <- logdensity.pois(y[i], phi[i])
                 ySim[i] ~ dpois(phi[i])
              }
              
              Deviance <- -2 * sum(log_lik[1:N])
          }"
    
    P = ncol(X)
    write_lines(jags_neg_LASSO, "jags_neg_LASSO.txt")
    jagsdata = list(X = X, y = y, N = length(y), P = ncol(X), idx = idx, nG = max(idx), k = as.vector(table(idx)))
    monitor = c("Intercept", "beta", "lambda", "Deviance", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = 0, 
                                              "beta" = rep(0, P), 
                                              "eta" = abs(jitter(rep(1, max(idx)), amount = 1)), 
                                              "psi" = abs(jitter(rep(1, max(idx)), amount = 1)), 
                                              "lambda" = 2, 
                                              "ySim" = sample(y, length(y)),
                                              .RNG.name= "lecuyer::RngStream", 
                                              .RNG.seed = sample(1:10000, 1)))
  }
  out = run.jags(model = "jags_neg_LASSO.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
  file.remove("jags_neg_LASSO.txt")
  return(out)
}
