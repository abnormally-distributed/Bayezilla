#' Bayesian Basic GLMs
#'
#' This model utilizes normal-gamma mixture priors. A student-t prior can be parameterized as a norma-gamma mixture 
#' by utilizing a gamma(df/2, df/2) distribution where df is the desired degrees of freedom. 
#' This model utilizes a single degree of freedom, which gives marginal cauchy distributions. The cauchy distribution has no defined 
#' first or second moments (mean and variance), hence does not require subjective specification of the prior means. The cauchy
#' distribution's extremely long tails allow coefficients with strong evidence of being large to not be shrunk too strongly, while the 
#' large probability mass at the mode of zero captures small noisy coefficients and regularizes them. The scale is modeled as a single hyperparameter "omega" via the gamma(.5, .5)
#' prior and this pooled-scale induces shrinkage when neccessary. 
#' This setup yields an ideal proper reference / objective prior that is data-driven and adaptive.
#' 
#' 
#' Standard gaussian, binomial, and poisson likelihood functions are available. 
#' 
#' Note that if you do not scale and center your numeric predictors, this will likely not perform well or
#' give reasonable results. The mixing hyperparameter omega assumes all covariates are on the same scale.
#' For an alternative proper reference prior see \code{\link[Bayezilla]{apcGlm}}
#' 
#'
#' @param formula the model formula
#' @param data a data frame
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
#' @return A run.jags object
#' @export
#'
#' @examples
#' glmBayes()
#'
#' @seealso \code{\link[Bayezilla]{apcGlm}}
#' 
glmBayes  = function(formula, data, family = "gaussian", log_lik = FALSE, iter=10000, warmup=1000, adapt=2000, chains=4, thin=1, method = "parallel", cl = makeCluster(2), ...){

  X = as.matrix(model.matrix(formula, data)[,-1])
  y = model.frame(formula, data)[,1]

  
  if (family == "gaussian"){

    jags_glm = "model{
              tau ~ dgamma(.01, .01)

              ## omega is the hyper prior on the beta precision specified to yield
              ## independent marginal student-t distributions
              omega ~ dgamma(.5,  .5)

              for (p in 1:P){
                beta[p] ~ dnorm(0, omega)
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

    P = ncol(X)
    write_lines(jags_glm, "jags_glm.txt")
    jagsdata = list(X = X, y = y,  N = length(y), P = ncol(X))
    monitor = c("Intercept", "beta", "sigma", "omega", "Deviance", "ySim" ,"log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = 0, "beta" = jitter(rep(0, P), amount = 1), "tau" = 1, "omega" = 1, "ySim" = y, .RNG.name= "lecuyer::RngStream", .RNG.seed = sample(1:10000, 1)))
  }

  if (family == "binomial" || family == "logistic"){

    jags_glm = "model{

              ## omega is the hyper prior on the beta precision specified to yield
              ## independent marginal cauchy distributions
              omega ~ dgamma(.5,  .5)

              for (p in 1:P){
                beta[p] ~ dnorm(0, omega)
              }

              Intercept ~ dnorm(0, 1)

              for (i in 1:N){
                 logit(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                 y[i] ~ dbern(psi[i])
                 log_lik[i] <- logdensity.bern(y[i], psi[i])
                 ySim[i] ~ dbern(psi[i])
              }
             Deviance <- -2 * sum(log_lik[1:N])
          }"

    P = ncol(X)
    write_lines(jags_glm, "jags_glm.txt")
    jagsdata = list(X = X, y = y, N = length(y),  P = ncol(X))
    monitor = c("Intercept", "beta", "omega", "Deviance", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = 0, "beta" = jitter(rep(0, P), amount = 1), "omega" = 1, "ySim" = y, .RNG.name= "lecuyer::RngStream", .RNG.seed= sample(1:10000, 1)))
  }

  if (family == "poisson"){

    jags_glm = "model{

              ## omega is the hyper prior on the beta precision specified to yield
              ## independent marginal cauchy distributions
              omega ~ dgamma(.5,  .5)

              for (p in 1:P){
                  beta[p] ~ dnorm(0, omega)
              }

              Intercept ~ dnorm(0, 1)

              for (i in 1:N){
                 log(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                 y[i] ~ dpois(psi[i])
                 log_lik[i] <- logdensity.pois(y[i], psi[i])
                 ySim[i] ~ dpois(psi[i])
              }

              Deviance <- -2 * sum(log_lik[1:N])
          }"

    write_lines(jags_glm, "jags_glm.txt")
    P = ncol(X)
    jagsdata = list(X = X, y = y, N = length(y),  P = ncol(X))
    monitor = c("Intercept", "beta", "omega", "Deviance", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = 0, "omega" = 1, "ySim" = y, .RNG.name= "lecuyer::RngStream", .RNG.seed= sample(1:10000, 1),"beta" = jitter(rep(0, P), amount = 1)))
  }

  out = run.jags(model = "jags_glm.txt", modules = c("glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE,...)
  if (is.null(cl) == FALSE){
    parallel::stopCluster(cl = cl)
  }
  file.remove("jags_glm.txt")
  return(out)
}
