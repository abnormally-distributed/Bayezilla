#' Bayesian UN-regularized GLMs
#'
#' @description This model utilizes extremely vague uniform priors on all parameters. 
#' These are essentially, for practical purposes, improper flat priors. While JAGS 
#' does not support an improper prior over the whole real line, the coefficients 
#' have normal(0, 1e-200) priors. A precision of 1e-200 is equivalent to a standard 
#' deviation of 1e100, or a googol.  No data that would ever be observed would be 
#' anywhere near this scale, so this can be used 
#' with scaled or unscaled data with equal validity, though I still recommend 
#' standardizing because MCMC samplers work better with data all on the same scale. 
#' The results from this should be equivalent to using lm() or glm()
#' in R but with the benefit of being able to make probability statements from the
#' marginal posterior ditributions. \cr
#' \cr
#' 
#' Standard gaussian, binomial, and poisson likelihood functions are available. \cr
#' \cr
#' Model Specification:
#' \cr
#' \if{html}{\figure{flat.png}{}}
#' \if{latex}{\figure{flat.png}{}}
#' 
#' @param formula the model formula
#' @param data a data frame
#' @param family one of "gaussian", "binomial", or "poisson".
#' @param log_lik Should the log likelihood be monitored? The default is FALSE.
#' @param iter How many post-warmup samples? Defaults to 10000.
#' @param warmup How many warmup samples? Defaults to 1000.
#' @param adapt How many adaptation steps? Defaults to 1000.
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
#' glmFlat()
#'
#' @seealso 
#' \code{\link[Bayezilla]{apcGlm}}
#' \code{\link[Bayezilla]{glmBayes}}
#' 
glmFlat  = function(formula, data, family = "gaussian", log_lik = FALSE, iter=10000, warmup=1000, adapt=1000, chains=4, thin=1, method = "parallel", cl = makeCluster(2), ...){
  
  X = as.matrix(model.matrix(formula, data)[,-1])
  y = model.frame(formula, data)[,1]
  
  
  if (family == "gaussian"){
    
    jags_glm = "model{
              tau ~  dunif(0, 1e200)
              
              for (p in 1:P){
                beta[p] ~ dnorm(0, 1e-200)
              }
              
              Intercept ~ dnorm(0, 1e-200)

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
    inits = lapply(1:chains, function(z) list("Intercept" = as.vector(coef(lm(formula, data)))[1], 
                                              "beta" = as.vector(coef(lm(formula, data)))[-1], 
                                              "tau" = 1, 
                                              "ySim" = y, 
                                              .RNG.name= "lecuyer::RngStream", 
                                              .RNG.seed = sample(1:20000, 1)))
  }
  
  if (family == "binomial" || family == "logistic"){
    
    jags_glm = "model{
    
              for (p in 1:P){
                beta[p] ~ dnorm(0, 1e-200)
              }

              Intercept ~ dnorm(0, 1e-200)

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
    inits = lapply(1:chains, function(z) list("Intercept" = as.vector(coef(glm(formula, data, family = "binomial")))[1], 
                                              "beta" = as.vector(coef(glm(formula, data, family = "binomial")))[-1], 
                                              "ySim" = y, 
                                              .RNG.name= "lecuyer::RngStream", 
                                              .RNG.seed= sample(1:20000, 1)))
  }
  
  if (family == "poisson"){
    
    jags_glm = "model{

              for (p in 1:P){
                  beta[p] ~ dnorm(0, 1e-200)
              }

              Intercept ~ dnorm(0, 1e-200)

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
    inits = lapply(1:chains, function(z) list("Intercept" = as.vector(coef(glm(formula, data, family = "poisson")))[1], 
                                              "ySim" = y,
                                              "beta" = as.vector(coef(glm(formula, data, family = "poisson")))[-1],
                                              .RNG.name= "lecuyer::RngStream", 
                                              .RNG.seed= sample(1:20000, 1)))
  }
  
  out = run.jags(model = "jags_glm.txt", modules = c("glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE,...)
  if (is.null(cl) == FALSE){
    parallel::stopCluster(cl = cl)
  }
  file.remove("jags_glm.txt")
  return(out)
}
