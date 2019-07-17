#' Bayesian GLMs
#'
#' @description This model utilizes Student-t priors for the coefficients. It is parameterized
#' as a scale mixture of Gaussians, which works by setting a gamma prior on the precision of the 
#' normal with shape and rate parameters equal to the desired degrees of freedom divided by two (the
#' shape parameter can also be multiplied by a variance to give the resulting student-t distribution
#' the appropriately scaled precision) :  
#' \if{html}{\figure{dof.png}{}}
#' \if{latex}{\figure{dof.png}{}}
#' \cr \cr
#' By default, the Cauchy distribution is obtained by setting the degrees of freedom to 1. The default 
#' prior standard deviation is also 1. While JAGS has a Student-t distribution built in, 
#' the samplers in JAGS tend to sample from Student-t 
#' distributions very inefficiently at times. By parameterizing the Student-t 
#' distribution this way JAGS can take advantage of the normal-gamma conjugacy with 
#' Gibbs sampling and sample very quickly and accurately. \cr
#' \cr
#' The default prior settings assume your data are standardized to mean zero and standard deviation of 1. 
#' \cr \cr
#' For an alternative prior that also performs very well see \code{\link[Bayezilla]{apcGlm}}, which when
#' lambda = -1 gives the Zellner-Siow Cauchy g-prior. \cr
#' \cr
#' The full model structure is given below: \cr
#' \cr
#' \cr
#' \if{html}{\figure{glm.png}{}}
#' \if{latex}{\figure{glm.png}{}}
#' \cr
#' \cr
#' The relationship between the degrees of freedom and the prior precision 
#' (inverse of squared prior standard deviation)
#' is represented in the figure below. Essentially, if the precision 
#' is sufficiently small (or equivalently, if prior s.d. is large) the amount of shrinkage
#' is essentially none regardless of the degrees of freedom. Increasing the degrees of freedom
#' also leads to greater regularization. However, the pattern of regularization differs for different
#' values, which is not represented in the figure. Supposing the precision is held constant at 1, 
#' degrees of freedom  < 8 will tend to shrink smaller 
#' coefficients to a greater degree, while larger coefficients will be affected less due to the long tails. However, 
#' as the degrees of freedom increase the Student-t distribution will take on an increasingly gaussian shape and the
#' tails will be pulled in, giving more uniform shrinkage (at this point it is effectively a ridge regression). If the
#' precision is increased, the contrast between small and large coefficients will tend to be even greater for small
#' degrees of freedom, while higher degrees of freedom approach a highly regularized ridge solution. \cr  \cr
#' Note that the figure below is for conceptual illustrative purposes, and does not correspond to an
#' exact mathematical function. \cr
#' \cr
#' \if{html}{\figure{shrinkagelowres.png}{}}
#' \if{latex}{\figure{shrinkagelowres.png}{}}
#' 
#' 
#' @param formula the model formula
#' @param data a data frame
#' @param family one of "gaussian", "binomial", or "poisson".
#' @param df degrees of freedom for prior.
#' @param s The desired prior scale. Defaults to 1. Is automatically squared within the model so
#' select a number here on the standard deviation scale.
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
glmBayes  = function(formula, data, family = "gaussian", s = 1, df = 1, log_lik = FALSE, iter=10000, warmup=1000, adapt=2000, chains=4, thin=1, method = "parallel", cl = makeCluster(2), ...){
  
  X = as.matrix(model.matrix(formula, data)[,-1])
  y = model.frame(formula, data)[,1]
  
  if (family == "gaussian"){
    
    jags_glm = "model{
              tau ~ dgamma(0.01, 0.01) 
              eta ~ dgamma(df * 0.50, pow(s, 2) * (df * 0.50))
              
              for (p in 1:P){
                beta[p] ~ dnorm(0, eta)
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
    jagsdata = list(X = X, y = y,  N = length(y), P = ncol(X), s = s, df = df)
    monitor = c("Intercept", "beta", "sigma", "Deviance", "ySim" ,"log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = lmSolve(formula, data)[1], 
                                              "beta" = lmSolve(formula, data)[-1], 
                                              "tau" = 1, 
                                              "ySim" = y, 
                                              "eta" = 2, 
                                              .RNG.name= "lecuyer::RngStream", 
                                              .RNG.seed = sample(1:10000, 1)))
  }
  
  if (family == "binomial"){
    
    jags_glm = "model{
    
              eta ~ dgamma(df * 0.50, pow(s, 2) * (df * 0.50))
              
              for (p in 1:P){
                beta[p] ~ dnorm(0, eta)
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
    jagsdata = list(X = X, y = y, N = length(y),  P = ncol(X), s = s, df = df)
    monitor = c("Intercept", "beta", "Deviance", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = coef(glm(formula, data, family = "binomial"))[1], 
                                              "beta" = coef(glm(formula, data, family = "binomial"))[-1],
                                              "ySim" = y, 
                                              "eta" = 2, 
                                              .RNG.name= "lecuyer::RngStream", 
                                              .RNG.seed= sample(1:10000, 1)))
  }
  
  
  if (family == "poisson"){
    
    jags_glm = "model{

              eta ~ dgamma(df * 0.50, pow(s, 2) * (df * 0.50))
              
              for (p in 1:P){
                beta[p] ~ dnorm(0, eta)
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
    jagsdata = list(X = X, y = y, N = length(y),  P = ncol(X), s = s, df = df)
    monitor = c("Intercept", "beta", "Deviance", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = coef(glm(formula, data, family = "poissson"))[1], 
                                              "ySim" = y,
                                              "eta" = 2, 
                                              "beta" = coef(glm(formula, data, family = "poissson"))[-1], 
                                              .RNG.name= "lecuyer::RngStream", 
                                              .RNG.seed= sample(1:10000, 1)))
  }
  
  out = run.jags(model = "jags_glm.txt", modules = c("glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, n.chains = chains, summarise = FALSE,...)
  if (is.null(cl) == FALSE){
    parallel::stopCluster(cl = cl)
  }
  file.remove("jags_glm.txt")
  return(out)
}
