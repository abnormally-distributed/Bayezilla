#' Bayesian GLMs
#'
#' @description This model utilizes cauchy priors for the coefficients. The cauchy distribution
#' is a Student-t distribution with 1 degree of freedom. It has undefined moments, meaning that it has
#' undefined mean and undefined variance (although it does have a scale). Here the scale is set to a precision
#' of 0.01, meaning a scale of 10. However, the extremely long tails of the cauchy 
#' distribution make it extremely uninformative and ideal for use as a weakly informative prior to obtain
#' unbiased estimates so long as your data are standardized to mean zero and standard deviation of 1. However,
#' it should still be weakly informative for unstandardized data as well, but it's a good idea to standardize
#' when using MCMC methods for computational efficiency. 
#' \cr
#' For an alternative prior that also performs very well see \code{\link[Bayezilla]{apcGlm}}, which when
#' lambda = -1 gives the Zellner-Siow Cauchy g-prior.
#' \cr
#' The model structure is given below: \cr
#' \cr
#' \cr
#' \if{html}{\figure{glm.png}{}}
#' \if{latex}{\figure{glm.png}{}}
#' \cr
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
#' 
glmBayes  = function(formula, data, family = "gaussian", log_lik = FALSE, iter=10000, warmup=1000, adapt=2000, chains=4, thin=1, method = "parallel", cl = makeCluster(2), ...){

  X = as.matrix(model.matrix(formula, data)[,-1])
  y = model.frame(formula, data)[,1]

  if (family == "gaussian"){

    jags_glm = "model{
              tau ~ dgamma(.01, .01) 

              for (p in 1:P){
                beta[p] ~ dt(0, .01, 1)
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

    P = ncol(X)
    write_lines(jags_glm, "jags_glm.txt")
    jagsdata = list(X = X, y = y,  N = length(y), P = ncol(X))
    monitor = c("Intercept", "beta", "sigma", "Deviance", "ySim" ,"log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = lmSolve(formula, data)[1], 
                                              "beta" = lmSolve(formula, data)[-1], 
                                              "tau" = 1, 
                                              "ySim" = y, 
                                              .RNG.name= "lecuyer::RngStream", 
                                              .RNG.seed = sample(1:10000, 1)))
  }

  if (family == "binomial" || family == "logistic"){

    jags_glm = "model{

              for (p in 1:P){
                beta[p] ~ dt(0, .01, 1)
              }

              Intercept ~ dnorm(0, 1e-10)

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
    monitor = c("Intercept", "beta", "Deviance", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = coef(glm(formula, data, family = "binomial"))[1], 
                                              "beta" = coef(glm(formula, data, family = "binomial"))[-1],
                                              "ySim" = y, 
                                              .RNG.name= "lecuyer::RngStream", 
                                              .RNG.seed= sample(1:10000, 1)))
  }

  if (family == "poisson"){

    jags_glm = "model{

              for (p in 1:P){
                  beta[p] ~ dt(0, .01, 1)
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

    write_lines(jags_glm, "jags_glm.txt")
    P = ncol(X)
    jagsdata = list(X = X, y = y, N = length(y),  P = ncol(X))
    monitor = c("Intercept", "beta", "Deviance", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = coef(glm(formula, data, family = "poissson"))[1], 
                                              "ySim" = y,
                                              "beta" = coef(glm(formula, data, family = "poissson"))[-1], 
                                              .RNG.name= "lecuyer::RngStream", 
                                              .RNG.seed= sample(1:10000, 1)))
  }

  out = run.jags(model = "jags_glm.txt", modules = c("glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE,...)
  if (is.null(cl) == FALSE){
    parallel::stopCluster(cl = cl)
  }
  file.remove("jags_glm.txt")
  return(out)
}
