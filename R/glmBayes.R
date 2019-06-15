#' Bayesian Regularized GLMs
#'
#' This model utilizes normal-gamma mixture priors. A student-t prior can be parameterized as a norma-gamma mixture 
#' by utilizing a gamma(df/2, df/2) distribution on the precision, where df is the desired degrees of freedom. Smaller
#' degrees of freedom 
#' result in long tails which can capture larger coefficients with sufficient evidence while also providing regularization
#' for coefficients which lack strong evidence of being large. This prevents unrealistically large coefficient estimates
#' due to collinearity. Collinearity results in increased posterior variance due to a flatter likelihood function, and a side
#' effect of this is inflated coefficient estimates. \cr
#' \cr
#' This model utilizes a prior on the degrees of freedom, meaning the prior distribution is a mixture of student-t distributions
#' with different degrees of freedom. This setup yields a prior that is data-driven and adaptive. In other words, only as much
#' regularization is applied as is needed. The behavior of this estimator is nearly identical to ridge regression. \cr
#' \cr
#' Standard gaussian, binomial, and poisson likelihood functions are available. \cr
#' \cr
#' For an alternative prior that also performs very well see \code{\link[Bayezilla]{apcGlm}}
#' \cr
#' Note that if you do not scale and center your numeric predictors, this will likely not perform well or
#' give reasonable results. The mixing hyperparameters assume all covariates are on the same scale.
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
              
              df ~ dgamma(.01, .01)
              omega ~ dgamma(.5 * df,  .5 * df)

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
    monitor = c("Intercept", "beta", "sigma", "omega", "df", "Deviance", "ySim" ,"log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = 0, "beta" = jitter(rep(0, P), amount = 1), "df" = 3, "tau" = 1, "omega" = 1, "ySim" = y, .RNG.name= "lecuyer::RngStream", .RNG.seed = sample(1:10000, 1)))
  }
  
  if (family == "student_t"){
    
    jags_glm = "model{
              nu ~ dexp(1/27) T(3, )
              tau_raw ~ dgamma(.01, .01)
              tau <- pow((1/tau_raw) * ((nu-2)/nu), -1)
              
              df ~ dgamma(.01, .01)
              omega ~ dgamma(.5 * df,  .5 * df)

              for (p in 1:P){
                beta[p] ~ dnorm(0, omega)
              }

              Intercept ~ dnorm(0, 1)

              for (i in 1:N){
                 y[i] ~ dt(Intercept + sum(beta[1:P] * X[i,1:P]), tau, nu)
                 log_lik[i] <- logdensity.t(y[i], Intercept + sum(beta[1:P] * X[i,1:P]), tau, nu)
                 ySim[i] ~ dt(Intercept + sum(beta[1:P] * X[i,1:P]), tau, nu)
              }
              sigma <- sqrt(1/tau)
              Deviance <- -2 * sum(log_lik[1:N])
          }"
    
    P = ncol(X)
    write_lines(jags_glm, "jags_glm.txt")
    jagsdata = list(X = X, y = y,  N = length(y), P = ncol(X))
    monitor = c("Intercept", "beta", "sigma", "nu", "omega", "df", "Deviance", "ySim" ,"log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = 0, "beta" = jitter(rep(0, P), amount = 1), "df" = 3, "nu" = 3, "tau_raw" = 1, "omega" = 1, "ySim" = y, .RNG.name= "lecuyer::RngStream", .RNG.seed = sample(1:10000, 1)))
  }

  if (family == "binomial" || family == "logistic"){

    jags_glm = "model{

              df ~ dgamma(.01, .01)
              omega ~ dgamma(.5 * df,  .5 * df)

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
    monitor = c("Intercept", "beta", "omega", "df", "Deviance", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = 0, "beta" = jitter(rep(0, P), amount = 1),"df" = 3, "omega" = 1, "ySim" = y, .RNG.name= "lecuyer::RngStream", .RNG.seed= sample(1:10000, 1)))
  }

  if (family == "poisson"){

    jags_glm = "model{

              df ~ dgamma(.01, .01)
              omega ~ dgamma(.5 * df,  .5 * df)

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
    monitor = c("Intercept", "beta", "omega", "df", "Deviance", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = 0, "omega" = 1, "df" = 3, "ySim" = y, .RNG.name= "lecuyer::RngStream", .RNG.seed= sample(1:10000, 1),"beta" = jitter(rep(0, P), amount = 1)))
  }

  out = run.jags(model = "jags_glm.txt", modules = c("glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE,...)
  if (is.null(cl) == FALSE){
    parallel::stopCluster(cl = cl)
  }
  file.remove("jags_glm.txt")
  return(out)
}
