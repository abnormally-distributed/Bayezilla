#' Normal-exponential-gamma Bayesian LASSO
#'
#' @description This implements the normal-exponential-gamma "hyperlasso" of Griffin & Brown (2011). This model has independent normal priors
#' on each coefficient, whose precision is modeled by independent, predictor specific, exponential distributions. The exponential
#' distributions in turn have their respective rate parameters modeled through independent gamma(.5, 1 / lambda^2) distributions.
#' Lambda is a single top-level hyperparameter here given a gamma(0.50 , 0.20) prior. \cr 
#' \cr
#' The model specification is given below: \cr
#' \cr
#' \cr
#' \cr 
#' Model Specification:
#' \cr
#' \if{html}{\figure{negLASSO.png}{}}
#' \if{latex}{\figure{negLASSO.png}{}}
#' \cr
#' \cr
#' The normal-exponential-gamma (NEG) lasso
#' is very similar to the adaptive Bayesian Lasso (\code{\link[Bayezilla]{adaLASSO}}), which also makes use of a 
#' normal-exponential-gamma hierarchy, except that it is parameterized slightly differently. 
#' \cr
#' @references 
#' Griffin, J. E. and Brown, P. J. (2011), Bayesian Hyper‐LASSOs With Non-Convex Penalization. Australian & New Zealand Journal of Statistics, 53: 423-442. doi:10.1111/j.1467-842X.2011.00641.x
#'
#' @param formula the model formula
#' @param data a data frame.
#' @param family one of "gaussian", "binomial", or "poisson".
#' @param lambda.prior either "dmouch" (the default) or "gamma"
#' @param dof the degrees of freedom for the normal-exponential-gamma prior. Defaults to 0.50.
#' @param log_lik Should the log likelihood be monitored? The default is FALSE.
#' @param iter How many post-warmup samples? Defaults to 10000.
#' @param warmup How many warmup samples? Defaults to 1000.
#' @param adapt How many adaptation steps? Defaults to 2000.
#' @param chains How many chains? Defaults to 4.
#' @param thin Thinning interval. Defaults to 1.
#' @param method Defaults to "rjparallel". For another parallel option, choose "parallel" or "rjags" for a single core run.
#' @param cl Use parallel::makeCluster(# clusters) to specify clusters for the parallel methods. Defaults to 3.
#' @param ... Other arguments to run.jags.
#'
#' @return A run.jags object
#'
#' @examples
#' negLASSO()
#' 
#'
#' @export
negLASSO  = function(formula, data, family = "gaussian", lambda.prior = "dmouch" , dof = 0.5, log_lik = FALSE, iter=10000, warmup=1000, adapt=2000, chains=4, thin=1, method = "rjparallel", cl = makeCluster(3), ...){

  X = model.matrix(formula, data)[,-1]
  y = model.frame(formula, data)[,1]

  if (lambda.prior == "dmouch"){
    
    if (method == "parallel"){
      message("method switching to rjparallel to enable use of DuMouchley's prior")
      method <- "rjparallel"
    }
    
    if (family == "gaussian" || family == "normal") {
      
      jags_neg_LASSO = "model{

              tau ~ dgamma(.01, .01)

              lambda ~ dmouch(1)

              Intercept ~ dnorm(0, 1e-10)

              for (p in 1:P){
                eta[p] ~ dgamma(dof, 1 / pow(lambda,2))
                psi[p] ~ dexp(eta[p])
                beta[p] ~ dnorm(0, 1 / psi[p])
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
      jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), dof = dof)
      monitor <- c("Intercept", "beta", "sigma", "lambda", "Deviance", "ySim", "log_lik")
      if (log_lik == FALSE){
        monitor = monitor[-(length(monitor))]
      }
      inits <- lapply(1:chains, function(z) list("Intercept" = lmSolve(formula, data)[1], 
                                                 "beta" = lmSolve(formula, data)[-1], 
                                                 "eta" = abs(jitter(rep(1, P), amount = 1)), 
                                                 "psi" = abs(jitter(rep(1, P), amount = 1)), 
                                                 "lambda" = 2, 
                                                 "tau" = 1, 
                                                 "ySim" = sample(y, length(y)),
                                                 .RNG.name= "lecuyer::RngStream", 
                                                 .RNG.seed = sample(1:10000, 1)))
    }
    
    
    if (family == "binomial" || family == "logistic") {
      
      jags_neg_LASSO = "model{

              lambda ~ dmouch(1)

              Intercept ~ dnorm(0, 1e-10)

              for (p in 1:P){
                eta[p] ~ dgamma(dof, 1 / pow(lambda,2))
                psi[p] ~ dexp(eta[p])
                beta[p] ~ dnorm(0, 1 / psi[p])
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
      jagsdata = list(X = X, y = y, N = length(y), P = ncol(X), dof = dof)
      monitor = c("Intercept", "beta", "lambda", "Deviance", "ySim", "log_lik")
      if (log_lik == FALSE){
        monitor = monitor[-(length(monitor))]
      }
      inits = lapply(1:chains, function(z) list("Intercept" = 0, 
                                                "beta" = rep(0, P), 
                                                "eta" = abs(jitter(rep(1, P), amount = 1)), 
                                                "psi" = abs(jitter(rep(1, P), amount = 1)), 
                                                "lambda" = 2, 
                                                "ySim" = sample(y, length(y)),
                                                .RNG.name= "lecuyer::RngStream", 
                                                .RNG.seed = sample(1:10000, 1)))
    }
    
    if (family == "poisson") {
      
      jags_neg_LASSO = "model{

              lambda ~ dmouch(1)

              Intercept ~ dnorm(0, 1e-10)

              for (p in 1:P){
                eta[p] ~ dgamma(dof, 1 / pow(lambda,2))
                psi[p] ~ dexp(eta[p])
                beta[p] ~ dnorm(0, 1 / psi[p])
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
      jagsdata = list(X = X, y = y, N = length(y), P = ncol(X), dof = dof)
      monitor = c("Intercept", "beta", "lambda", "Deviance", "ySim", "log_lik")
      if (log_lik == FALSE){
        monitor = monitor[-(length(monitor))]
      }
      inits = lapply(1:chains, function(z) list("Intercept" = 0, 
                                                "beta" = rep(0, P), 
                                                "eta" = abs(jitter(rep(1, P), amount = 1)), 
                                                "psi" = abs(jitter(rep(1, P), amount = 1)), 
                                                "lambda" = 2, 
                                                "ySim" = sample(y, length(y)),
                                                .RNG.name= "lecuyer::RngStream", 
                                                .RNG.seed = sample(1:10000, 1)))
    
    }
    out = run.jags(model = "jags_neg_LASSO.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
    file.remove("jags_neg_LASSO.txt")
    return(out)
  }
if (lambda.prior == "gamma"){
  
  if (family == "gaussian" || family == "normal") {

  jags_neg_LASSO = "model{

              tau ~ dgamma(.01, .01)

              lambda ~ dgamma(0.5 , 0.20)

              Intercept ~ dnorm(0, 1e-10)

              for (p in 1:P){
                eta[p] ~ dgamma(dof, 1 / pow(lambda,2))
                psi[p] ~ dexp(eta[p])
                beta[p] ~ dnorm(0, 1 / psi[p])
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
  jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), dof = dof)
  monitor <- c("Intercept", "beta", "sigma", "lambda", "Deviance", "ySim", "log_lik")
  if (log_lik == FALSE){
    monitor = monitor[-(length(monitor))]
  }
  inits <- lapply(1:chains, function(z) list("Intercept" = lmSolve(formula, data)[1], 
                                             "beta" = lmSolve(formula, data)[-1], 
                                             "eta" = abs(jitter(rep(1, P), amount = 1)), 
                                             "psi" = abs(jitter(rep(1, P), amount = 1)), 
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

              for (p in 1:P){
                eta[p] ~ dgamma(dof, 1 / pow(lambda,2))
                psi[p] ~ dexp(eta[p])
                beta[p] ~ dnorm(0, 1 / psi[p])
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
    jagsdata = list(X = X, y = y, N = length(y), P = ncol(X), dof = dof)
    monitor = c("Intercept", "beta", "lambda", "Deviance", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = 0, 
                                              "beta" = rep(0, P), 
                                              "eta" = abs(jitter(rep(1, P), amount = 1)), 
                                              "psi" = abs(jitter(rep(1, P), amount = 1)), 
                                              "lambda" = 2, 
                                              "ySim" = sample(y, length(y)),
                                              .RNG.name= "lecuyer::RngStream", 
                                              .RNG.seed = sample(1:10000, 1)))
  }

  if (family == "poisson") {

    jags_neg_LASSO = "model{

              lambda ~ dgamma(0.5 , 0.20)

              Intercept ~ dnorm(0, 1e-10)

              for (p in 1:P){
                eta[p] ~ dgamma(dof, 1 / pow(lambda,2))
                psi[p] ~ dexp(eta[p])
                beta[p] ~ dnorm(0, 1 / psi[p])
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
    jagsdata = list(X = X, y = y, N = length(y), P = ncol(X), dof = dof)
    monitor = c("Intercept", "beta", "lambda", "Deviance", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = 0, 
                                              "beta" = rep(0, P), 
                                              "eta" = abs(jitter(rep(1, P), amount = 1)), 
                                              "psi" = abs(jitter(rep(1, P), amount = 1)), 
                                              "lambda" = 2, 
                                              "ySim" = sample(y, length(y)),
                                              .RNG.name= "lecuyer::RngStream", 
                                              .RNG.seed = sample(1:10000, 1)))
  }
  out = run.jags(model = "jags_neg_LASSO.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
  file.remove("jags_neg_LASSO.txt")
  return(out)
  }
}
