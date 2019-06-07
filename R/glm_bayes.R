#' Bayesian Basic GLMs
#'
#' This model utilizes normal-gamma mixture priors. A student-t prior can be parameterized as a norma-gamma mixture 
#' by utilizing a gamma(nu/2, nu/2) distribution where nu is the desired degrees of freedom. By default,
#' this model utilizes a single degree of freedom. One degree of freedom yields gamma(.5, .5), which gives marginal cauchy distributions. Hence, 
#' this model results in  marginal independent cauchy priors on each coefficient. The cauchy distribution has no defined 
#' first or second moments (mean and variance) and hence is an ideal proper reference prior. The cauchy distribution's 
#' extremely long tails allow coefficients with strong evidence of being large to not be shrunk too strongly, while the 
#' large probability mass at the mode of zero captures small noisy coefficients and regularizes them. However,
#' at times the cauchy distribution can be difficult to sample from, so I have enabled the user to select from either 1, 3, 6,
#' or 12 degrees of freedom. Larger degrees of freedom result in stronger regularization. Essentially, the degrees
#' of freedom controls the prior values of the omega hyperparameter (precision on the coefficient priors). Nevertheless,
#' this is only a *prior*. The model will learn from the data, and this adaptive shrinkage property results in an ideal prior
#' that, other than subjective choice of selecting the degrees of freedom, is completely data driven. However, one could always
#' usem mode comparison to select a degrees of freedom.
#' 
#' Standard gaussian, binomial, and poisson likelihood functions are available. 
#' 
#' Note that if you do not scale and center your numeric predictors, this will likely not perform well or
#' give reasonable results. The mixing hyperparameter omega assumes all covariates are on the same scale.
#' 
#'
#' @param formula the model formula
#' @param data a data frame
#' @param family one of "gaussian", "binomial", or "poisson".
#' @param nu the degrees of freedom for the marginal student-t priors on the coefficients 
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
glmBayes  = function(formula, data, family = "gaussian", log_lik = FALSE, iter=10000, warmup=1000, adapt=2000, chains=4, thin=1, method = "parallel", cl = makeCluster(2), ...){

  X = as.matrix(model.matrix(formula, data)[,-1])
  y = model.frame(formula, data)[,1]
  
  if (family == "gaussian"){

    jags_glm = "model{
              tau ~ dgamma(.001, .001)

              ## omega is the hyper prior on the beta precision specified to yield
              ## independent marginal cauchy distributions for non-informativeness
              omega ~ dgamma(nu * .5, nu *  .5)

              for (p in 1:P){
                beta[p] ~ dnorm(0, omega)
              }

              Intercept ~ dnorm(0, .01)

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
    inits = lapply(1:chains, function(z) list("Intercept" =0, "beta" = jitter(rep(0, P), amount = 1), "tau" = 1, "omega" = 1, "ySim" = y, .RNG.name= "lecuyer::RngStream", .RNG.seed = sample(1:10000, 1)))
  }

  if (family == "binomial" || family == "logistic"){

    jags_glm = "model{

              ## omega is the hyper prior on the beta precision specified to yield
              ## independent marginal cauchy distributions for non-informativeness
              omega ~ dgamma(nu * .5, nu *  .5)

              for (p in 1:P){
                beta[p] ~ dnorm(0, omega)
              }

              Intercept ~ dnorm(0, .01)

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
    jagsdata = list(X = X, y = y, N = length(y), P = ncol(X))
    monitor = c("Intercept", "beta", "omega", "Deviance", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = 0, "beta" = jitter(rep(0, P), amount = 1), "omega" = 1, "ySim" = y, .RNG.name= "lecuyer::RngStream", .RNG.seed= sample(1:10000, 1)))
  }

  if (family == "poisson"){

    jags_glm = "model{

              ## omega is the hyper prior on the beta precision specified to yield
              ## independent marginal cauchy distributions for non-informativeness
              omega ~ dgamma(nu * .5, nu *  .5)

              for (p in 1:P){
                  beta[p] ~ dnorm(0, omega)
              }

              Intercept ~ dnorm(0, .01)

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
