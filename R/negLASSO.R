#' Negative-exponential-gamma Bayesian LASSO
#'
#' This implements the normal-exponential-gamma "hyperlasso" of Griffin & Brown (2011). This model has independent normal priors
#' on each coefficient, whose precision is modeled by independent, predictor specific, exponential distributions. The exponential
#' distributions in turn have their respective rate parameters modeled through independent gamma(.5, 1 / lambda^2) distributions.
#' Lambda is a single top-level hyperparameter here given a gamma(0.5 , 0.05) prior. \cr 
#' \cr
#' The model specification is given below: \cr
#' \cr
#' #### Top Level Parameters #### \cr 
#' \cr
#' tau ~ gamma(.01, .01) # only for Gaussian likelihood \cr
#' lambda ~ gamma(0.5 , 0.05) \cr
#' Intercept ~ normal(0, 1) \cr
#' 
#' #### Predictor Level Parameters #### \cr
#' eta_i ~ gamma(.5, 1 / lambda^2) \cr
#' psi_i ~ exponential(eta_i) \cr
#' beta_i ~ normal(0, 1 / psi_i) # psi is inverted here because JAGS parameterizes by the precision \cr
#' \cr
#' \cr
#' The normal-exponential-gamma (NEG) lasso
#' is very similar to the adaptive Bayesian Lasso (\code{\link[Bayezilla]{adaBLASSO}}), which also makes use of a 
#' normal-exponential-gamma hierarchy, except that it is parameterized slightly differently. 
#' \cr
#' @references 
#' Griffin, J. E. and Brown, P. J. (2011), Bayesian Hyper‐LASSOs With Non-Convex Penalization. Australian & New Zealand Journal of Statistics, 53: 423-442. doi:10.1111/j.1467-842X.2011.00641.x
#'
#' @param formula the model formula
#' @param data a data frame.
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
#' @export
#'
#' @examples
#' negLASSO()
#'
negLASSO  = function(formula, data, family = "gaussian", log_lik = FALSE, iter=10000, warmup=1000, adapt=2000, chains=4, thin=1, method = "parallel", cl = makeCluster(3), ...){

  X = model.matrix(formula, data)[,-1]
  y = model.frame(formula, data)[,1]

  if (family == "gaussian" || family == "normal") {

  jags_neg_LASSO = "model{

              tau ~ dgamma(.01, .01)

              lambda ~ dgamma(0.5 , 0.05)

              Intercept ~ dnorm(0, 1)

              for (p in 1:P){
                eta[p] ~ dgamma(.5, 1 / pow(lambda,2))
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
  jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X))
  monitor <- c("Intercept", "beta", "sigma", "lambda", "eta", "Deviance", "psi", "ySim", "log_lik")
  if (log_lik == FALSE){
    monitor = monitor[-(length(monitor))]
  }
  inits <- lapply(1:chains, function(z) list("Intercept" = 0, "beta" = rep(0, P), "eta" = abs(jitter(rep(1, P), amount = 1)), "psi" = abs(jitter(rep(1, P), amount = 1)), "lambda" = 10, "tau" = 1, .RNG.name= "lecuyer::RngStream", .RNG.seed = sample(1:10000, 1)))
}


  if (family == "binomial" || family == "logistic") {

    jags_neg_LASSO = "model{

              lambda ~ dgamma(.5, .05)

              Intercept ~ dnorm(0, 1)

              for (p in 1:P){
                eta[p] ~ dgamma(.5, 1 / lambda^2)
                psi[p] ~ dexp(eta[p])
                beta[p] ~ dnorm(0, 1 / psi[p])
              }
              for (i in 1:N){
                 logit(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                 y[i] ~ dbern(psi[i])
                 log_lik[i] <- logdensity.bern(y[i], psi[i])
                 ySim[i] ~ dbern(psi[i])
              }
              sigma <- sqrt(1/tau)
              Deviance <- -2 * sum(log_lik[1:N])
          }"

    P = ncol(X)
    write_lines(jags_neg_LASSO, "jags_neg_LASSO.txt")
    jagsdata = list(X = X, y = y, N = length(y), P = ncol(X))
    monitor = c("Intercept", "beta", "lambda", "eta",  "Deviance", "psi", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = 0, "beta" = rep(0, P), "eta" = abs(jitter(rep(1, P), amount = 1)), "psi" = abs(jitter(rep(1, P), amount = 1)), "lambda" = 10, .RNG.name= "lecuyer::RngStream", .RNG.seed = sample(1:10000, 1)))
  }

  if (family == "poisson") {

    jags_neg_LASSO = "model{

              lambda ~ dgamma(.5, .05)

              Intercept ~ dnorm(0, 1)

              for (p in 1:P){
                eta[p] ~ dgamma(.5, 1 / lambda^2)
                psi[p] ~ dexp(eta[p])
                beta[p] ~ dnorm(0, 1 / psi[p])
              }
              for (i in 1:N){
                 log(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                 y[i] ~ dpois(psi[i])
                 log_lik[i] <- logdensity.pois(y[i], psi[i])
                 ySim[i] ~ dpois(psi[i])
              }
              sigma <- sqrt(1/tau)
              Deviance <- -2 * sum(log_lik[1:N])
          }"

    P = ncol(X)
    write_lines(jags_neg_LASSO, "jags_neg_LASSO.txt")
    jagsdata = list(X = X, y = y, N = length(y), P = ncol(X))
    monitor = c("Intercept", "beta", "lambda", "eta",  "Deviance", "psi", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = 0, "beta" = rep(0, P), "eta" = abs(jitter(rep(1, P), amount = 1)), "psi" = abs(jitter(rep(1, P), amount = 1)), "lambda" = 10, .RNG.name= "lecuyer::RngStream", .RNG.seed = sample(1:10000, 1)))
  }
  out = run.jags(model = "jags_neg_LASSO.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
  file.remove("jags_neg_LASSO.txt")
  return(out)
}
