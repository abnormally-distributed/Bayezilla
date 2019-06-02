#' Regression Models for Location, Scale, and Shape
#'
#' This is a very limited set models for estimating generalized linear models for location, scale, and shape.
#' Essentially what one does is model the scale parameter (see the examples for usage and syntax).
#' This is very useful if you have a heteroscedastic relationship between the linear predictor 
#' and the response variable. One can also model the scale of the response as a function of
#' a grouping variable, such as that shown in the examples where the Species factor variable
#' is used to model the standard deviation. Due to the limited set of distributions available
#' in JAGS, the support for the glmlss / gamlss family of models is extremely limited. However,
#' for more flexibility check out the bamlss package, which allows estimation of a wide class
#' of likelihood functions including skewed normal, power exponential (generalized normal), and
#' many more. This package only offers "gaussian" (the default), "student_t", and "laplace".
#' 
#' NOTE: Make sure to center and scale the predictor variables. All priors on coefficients
#' for all model parameters are set to a default precision of .0625, implying a standard 
#' deviation of 4 on each coefficient's prior. This is weakly informative when the predictors
#' are unit-scale, but cease to be weakly informative when the predictors are not unit scale.
#'
#' @param formula the model formula
#' @param sigma.formula the formula for the scale parameter. Defaults to ~ 1 (intercept only)
#' @param nu.formula the formula for the normality parameter, if using family = "student_t". Defaults to ~ 1 (intercept only)
#' @param data a data frame
#' @param family one of "gaussian", "laplace", or "student_t".
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
#' lssBayes(Petal.Width ~ . - Species, sigma.formula = ~ Species, data = scale(iris))
#' lssBayes(Petal.Width ~ . - Sepal.Length - Species, sigma.formula = ~ Species, nu.formula = ~ Sepal.Length, data = scale(iris), family = "student_t")
#' 
lssBayes  = function(formula, sigma.formula = ~ 1, nu.formula = ~ 1, data, family = "gaussian", log_lik = FALSE, iter=10000, warmup=1000, adapt=2000, chains=4, thin=3, method = "parallel", cl = makeCluster(2), ...){
  
  RNGlist = c("base::Wichmann-Hill", "base::Marsaglia-Multicarry", "base::Super-Duper", "base::Mersenne-Twister")
  if (chains > 4){
    chains = 4
  }  
  
if (family == "gaussian"){
  
  jags_gamlss = "model{
              
              Intercept ~ dnorm(0, .0625)
              
              for (p in 1:P){
                beta[p] ~ dnorm(0, .0625)
              }
              
              for (s in 1:S){
                sigma_beta[s] ~ dnorm(0, .0625)
              }
              
              for (i in 1:N){
                 log(sigma[i]) <- sum(sigma_beta[1:S] * SX[i, 1:S])
                 tau[i] <- 1 / pow(sigma[i], 2)
                 y[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau[i])
                 log_lik[i] <- logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]), tau[i])
                 ySim[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau[i])
              }
              deviance <- -2 * sum(log_lik[1:N])
          }"
  
  X = as.matrix(model.matrix(formula, data)[,-1])
  y = model.frame(formula, data)[,1]
  SX = as.matrix(model.matrix(sigma.formula, data))
  P = ncol(X)
  S = ncol(SX)
  write_lines(jags_gamlss, "jags_gamlss.txt")
  jagsdata = list(X = X, SX = SX, y = y,  N = length(y), P = P, S = S)
  monitor = c("Intercept", "beta", "sigma_beta", "deviance", "ySim" ,"log_lik")
  if (log_lik == FALSE){
    monitor = monitor[-(length(monitor))]
  }
  inits = lapply(1:chains, function(z) list("Intercept" =0, "beta" = jitter(rep(0, P), amount = 1), "sigma_beta" = jitter(rep(0, S), amount = 1), "ySim" = y, .RNG.name=RNGlist[z], .RNG.seed = sample(1:10000, 1)))
}
  
  if (family == "laplace"){
    
    jags_gamlss = "model{
              
              Intercept ~ ddexp(0, .0625)
              
              for (p in 1:P){
                beta[p] ~ dnorm(0, .0625)
              }
              
              for (s in 1:S){
                scale_beta[s] ~ ddexp(0, .0625)
              }
              
              for (i in 1:N){
                 log(scale[i]) <- sum(scale_beta[1:S] * SX[i, 1:S])
                 tau[i] <- 1 / pow(scale[i], 2)
                 y[i] ~ ddexp(Intercept + sum(beta[1:P] * X[i,1:P]), tau[i])
                 log_lik[i] <- logdensity.dexp(y[i], Intercept + sum(beta[1:P] * X[i,1:P]), tau[i])
                 ySim[i] ~ ddexp(Intercept + sum(beta[1:P] * X[i,1:P]), tau[i])
              }
              deviance <- -2 * sum(log_lik[1:N])
          }"
    
    X = as.matrix(model.matrix(formula, data)[,-1])
    y = model.frame(formula, data)[,1]
    SX = as.matrix(model.matrix(sigma.formula, data))
    P = ncol(X)
    S = ncol(SX)
    write_lines(jags_gamlss, "jags_gamlss.txt")
    jagsdata = list(X = X, SX = SX, y = y,  N = length(y), P = P, S = S)
    monitor = c("Intercept", "beta", "scale_beta", "deviance", "ySim" ,"log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" =0, "beta" = jitter(rep(0, P), amount = 1), "scale_beta" = jitter(rep(0, S), amount = 1), "ySim" = y, .RNG.name=RNGlist[z], .RNG.seed = sample(1:10000, 1)))
  }
  
  
  if (family == "student_t"){
    
    jags_gamlss = "model{
              
              Intercept ~ dnorm(0, .0625)
              
              for (p in 1:P){
                beta[p] ~ dnorm(0, .0625)
              }
              
              for (s in 1:S){
                nu_beta[s] ~ dnorm(0, .0625)
              }
              
              for (v in 1:V){
                scale_beta[v] ~ dnorm(0, .0625)
              }
              
              for (i in 1:N){
                 log(scale[i]) <- sum(scale_beta[1:S] * SX[i, 1:S])
                 tau[i] <- 1 / pow(scale[i], 2)
                 log(nu[i]) <- sum(nu_beta[1:V] * VX[i, 1:V])
                 y[i] ~ ddexp(Intercept + sum(beta[1:P] * X[i,1:P]), tau[i])
                 log_lik[i] <- logdensity.dexp(y[i], Intercept + sum(beta[1:P] * X[i,1:P]), tau[i])
                 ySim[i] ~ ddexp(Intercept + sum(beta[1:P] * X[i,1:P]), tau[i])
              }
              deviance <- -2 * sum(log_lik[1:N])
          }"
    
    X = as.matrix(model.matrix(formula, data)[,-1])
    y = model.frame(formula, data)[,1]
    SX = as.matrix(model.matrix(sigma.formula, data))
    VX = as.matrix(model.matrix(nu.formula, data))
    P = ncol(X)
    S = ncol(SX)
    V = nocl(VX)
    write_lines(jags_gamlss, "jags_gamlss.txt")
    jagsdata = list(X = X, SX = SX, VX = VX, y = y,  N = length(y), P = P, S = S, V=V)
    monitor = c("Intercept", "beta", "scale_beta", "nu_beta", "deviance", "ySim" ,"log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" =0, "beta" = jitter(rep(0, P), amount = 1), "scale_beta" = jitter(rep(0, S), amount = 1), "nu_beta" = jitter(rep(0, V), amount = 1), "ySim" = y, .RNG.name=RNGlist[z], .RNG.seed = sample(1:10000, 1)))
  }
  
  out = run.jags(model = "jags_gamlss.txt", modules = "glm", monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE,...)
  if (is.null(cl) == FALSE){
    parallel::stopCluster(cl = cl)
  }
  file.remove("jags_gamlss.txt")
  return(out)
}
