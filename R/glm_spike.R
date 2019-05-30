#' Bernoulli-Normal Mixture Selection for GLMs
#'
#' @description This is the most basic type of Bayesian variable selection. This models the
#' regression coefficients as coming from either a null distribution represented
#' by a probability mass of 100% at zero, or as coming from a normal distribution. A student-t family can be obtained
#' as a norma-gamma mixture by utilizing a gamma(nu/2, nu/2) distribution where nu is the desired degrees of freedom.
#' One degree of freedom yields gamma(.5, .5), which is the cauchy distribution. Hence, this model results in 
#' marginal independent cauchy priors on each coefficient. The cauchy distribution has no defined first or second moments
#' (mean and variance) and hence is an ideal proper reference prior. The cauchy distribution's extremely
#' long tails allow coefficients with strong evidence of being large to not be shrunk too strongly, while the large
#' probability mass at the mode of zero captures small noisy coefficients and regularizes them. This adaptive shrinkage
#' property results in an ideal prior.
#' 
#' Note that you should center and scale your predictor variables before using this function.
#' If you do not scale and center your numeric predictors, this will likely not perform well or
#' give reasonable results. The mixing hyperparameter omega assumes all covariates are on the same scale.
#'
#' This model works best on smaller to medium sized data sets with a small number of variables (less than 20). 
#' If you experience difficulty with running times or obtaining a good effective sample size consider using the extended LASSO. 
#' Another tip that may improve performance is using a beta(1, 1) or beta(1.5, 1.5) prior on phi. If all coefficients are set to
#' zero and there is truly good reason to believe this is not correct (ie, you aren't just hunting for "statistical
#' significance", are you?) a prior such as beta(4, 2) may help. If there is no sparsity, and you have genuine reason
#' to believe there should be try a beta(2,8) prior (if not, ask yourself why are you using variable selection?
#' If it is to deal with collinearity or other foibles for least squares, consider trying the glm_bayes function instead).
#'
#' @param formula the model formula
#' @param data a data frame
#' @param family one of "gaussian", "binomial", or "poisson".
#' @param phi_prior The beta distribution parameters on the inclusion probabilities. Default is Jeffrey's prior, c(0.5, 0.5).
#' @param log_lik Should the log likelihood be monitored? The default is FALSE.
#' @param iter How many post-warmup samples? Defaults to 10000.
#' @param warmup How many warmup samples? Defaults to 1000.
#' @param adapt How many adaptation steps? Defaults to 2000.
#' @param chains How many chains? Defaults to 4.
#' @param thin Thinning interval. Defaults to 3.
#' @param method Defaults to "parallel". For an alternative parallel option, choose "rjparallel" or. Otherwise, "rjags" (single core run).
#' @param cl Use parallel::makeCluster(# clusters) to specify clusters for the parallel methods. Defaults to two cores.
#' @param ... Other arguments to run.jags.
#'
#' @return A run.jags object
#' @export
#'
#' @examples
#' glmSpike()
#'
Spike  = function(formula, data, family = "gaussian", phi_prior = c(.5, .5), log_lik = FALSE, iter=10000, warmup=1000, adapt=2000, chains=4, thin=3, method = "parallel", cl = makeCluster(2), ...){

  X = model.matrix(formula, data)[,-1]
  y = model.frame(formula, data)[,1]

  RNGlist = c("base::Wichmann-Hill", "base::Marsaglia-Multicarry", "base::Super-Duper", "base::Mersenne-Twister")
  if (chains > 4){
    chains = 4
  }
  
  if (family == "gaussian"){

    jags_glm_spike = "model{
              tau ~ dgamma(.001, .001)
              phi ~ dbeta(a, b)
              omega ~ dgamma(.5, .5)
              for (p in 1:P){
                theta[p] ~ dnorm(0, omega)
                delta[p] ~ dbern(phi)
                beta[p] <- delta[p] * theta[p]
              }

              Intercept ~ dnorm(0, .01)

              for (i in 1:N){
                 y[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 log_lik[i] <- logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 ySim[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
              }
              sigma <- sqrt(1/tau)
              deviance <- -2 * sum(log_lik[1:N])
          }"

    P = ncol(X)
    write_lines(jags_glm_spike, "jags_glm_spike.txt")
    jagsdata = list(X = X, y = y, N = length(y), P = ncol(X), a = phi_prior[1], b = phi_prior[2])
    monitor = c("Intercept", "beta", "sigma", "phi", "omega", "delta", "deviance", "ySim" , "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = 0, .RNG.name= RNGlist[z], .RNG.seed= sample(1:10000, 1), "omega" = .0001,  "ySim" = y, "delta"=rep(1, P), "phi" = .20 , "theta" = jitter(rep(0, P), amount = .25), "tau" = 1))
    out = run.jags(model = "jags_glm_spike.txt", modules = "glm", monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, ...)
    return(out)
  }

  if (family == "binomial" || family == "logistic"){

    jags_glm_spike = "model{
              phi ~ dbeta(a, b)
              omega ~ dgamma(.5, .5)
              for (p in 1:P){
                theta[p] ~ dnorm(0, omega)
                delta[p] ~ dbern(phi)
                beta[p] <- delta[p] * theta[p]
              }
              Intercept ~ dnorm(0, .01)
              for (i in 1:N){
                 logit(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                 y[i] ~ dbern(psi[i])
                 log_lik[i] <- logdensity.bern(y[i], psi[i])
                 ySim[i] ~ dbern(psi[i])
              }
             deviance <- -2 * sum(log_lik[1:N])
          }"

    P = ncol(X)
    write_lines(jags_glm_spike, "jags_glm_spike.txt")
    y = as.numeric(as.factor(y)) - 1
    jagsdata = list(X = X, y = y, N = length(y), P = ncol(X), a = phi_prior[1], b = phi_prior[2])
    monitor = c("Intercept", "beta", "phi", "omega", "delta", "deviance", "ySim" , "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = 0, .RNG.name= RNGlist[z], .RNG.seed= sample(1:10000, 1), "omega" = .0001,  "ySim" = y, "delta" = rep(1, P), "phi" = .20 , "theta" = jitter(rep(0, P), amount = .25)))
    out = run.jags(model = "jags_glm_spike.txt", modules = "glm", monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, ...)
    return(out)
  }

  if (family == "poisson"){

    jags_glm_spike = "model{
              phi ~ dbeta(a, b)
              omega ~ dgamma(.5, .5)
              for (p in 1:P){
                theta[p] ~ dnorm(0, omega)
                delta[p] ~ dbern(phi)
                beta[p] <- delta[p] * theta[p]
              }
              Intercept ~ dnorm(0, .01)
              for (i in 1:N){
                 log(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                 y[i] ~ dpois(psi[i])
                 log_lik[i] <- logdensity.pois(y[i], psi[i])
                 ySim[i] ~ dpois(psi[i])
              }
              deviance <- -2 * sum(log_lik[1:N])
          }"

    write_lines(jags_glm_spike, "jags_glm_spike.txt")
    P = ncol(X)
    jagsdata = list(X = X, y = y, N = length(y),  P = ncol(X), a = phi_prior[1], b = phi_prior[2])
    monitor = c("Intercept", "beta", "phi", "omega", "delta", "deviance", "ySim" , "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = 0, .RNG.name= RNGlist[z], .RNG.seed= sample(1:10000, 1), "omega" = .0001,  "ySim" = y, "delta"=rep(1, P), "phi" = .20 , "theta" = jitter(rep(0, P), amount = .25)))
    out = run.jags(model = "jags_glm_spike.txt", modules = "glm", monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, ...)
    file.remove("jags_glm_spike.txt")
    if (!is.null(cl)) {
      parallel::stopCluster(cl = cl)
    }
    return(out)
  }
}
