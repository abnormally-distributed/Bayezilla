#' Bayesian Elastic Net with additional unpenalized design covariates
#'
#' @description The Bayesian elastic net described by Li and Lin (2010). This function has the further allowance for a set of covariates that are not penalized. 
#' For example, you may wish to include variables such as age and gender so that  the coefficients for the other variables are 
#' penalized while controlling for these. This is a common need in research.  For the binomial and poisson likelihoods
#' plug-in pseudo-variances are used. 
#' 
#' 
#' \cr
#' The model structure is given below: \cr
#' \cr
#' \cr
#' \if{html}{\figure{elasticNetDC.png}{}}
#' \if{latex}{\figure{elasticNetDC.png}{}}
#' \cr
#' \cr
#' \cr
#' Plugin Pseudo-Variances: \cr
#' \if{html}{\figure{pseudovar.png}{}}
#' \if{latex}{\figure{pseudovar.png}{}}
#' \cr
#' 
#' @param formula the model formula
#' @param design.formula formula for the design covariates.
#' @param data a data frame.
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
#' @references Li, Qing; Lin, Nan. The Bayesian elastic net. Bayesian Anal. 5 (2010), no. 1, 151--170. doi:10.1214/10-BA506. https://projecteuclid.org/euclid.ba/1340369796
#' 
#' @return A run.jags object
#' @export
#'
#' @examples
#' bayesEnetDC()
#'
bayesEnetDC  = function(formula, design.formula, data,  family = "gaussian", log_lik = FALSE, iter=10000, warmup=1000, adapt=2000, chains=4, thin=1, method = "parallel", cl = makeCluster(2), ...){
  
  X = model.matrix(formula, data)[,-1]
  y = model.frame(formula, data)[,1]
  FX <- as.matrix(model.matrix(design.formula, data)[, -1])

  if (family == "gaussian"){
    
    jags_elastic_net = "model{
              tau ~ dgamma(.01, .01)
              sigma <- sqrt(1/tau)
              lambdaL1 ~ dgamma(0.25 , 0.2)
              lambdaL2 ~ dgamma(0.25 , 0.2)

              Intercept ~ dnorm(0, 1e-10)

              for (p in 1:P){
                eta[p] ~ dgamma(.5, (8 * lambdaL2 * pow(sigma,2)) / pow(lambdaL1, 2)) T(1,)
                beta_prec[p] <- (lambdaL2/pow(sigma,2)) * (eta[p]/(eta[p]-1))
                beta[p] ~ dnorm(0, beta_prec[p])
              }
              
              for (f in 1:FP){
                 design_beta[f] ~ dnorm(0, 1e-200)
              }

              for (i in 1:N){
                 y[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P])  + sum(design_beta[1:FP] * FX[i,1:FP]), tau)
                 log_lik[i] <- logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]) + sum(design_beta[1:FP] * FX[i,1:FP]), tau)
                 ySim[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]) + sum(design_beta[1:FP] * FX[i,1:FP]), tau)
              }

              Deviance <- -2 * sum(log_lik[1:N])
          }"
  
  P <- ncol(X)
  FP <- ncol(FX)
  write_lines(jags_elastic_net, "jags_elastic_net.txt")
  jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), FP = FP, FX = FX)
  monitor <- c("Intercept", "beta", "design_beta", "sigma", "lambdaL1", "lambdaL2", "Deviance",  "ySim", "log_lik")
  if (log_lik == FALSE){
    monitor = monitor[-(length(monitor))]
  }
  inits <- lapply(1:chains, function(z) list("Intercept" = 0, 
                                             .RNG.name= "lecuyer::RngStream", 
                                             .RNG.seed= sample(1:10000, 1), 
                                             "beta" = rep(0, P), 
                                             "design_beta" = rep(0, FP),
                                             "eta" = 1 + abs(jitter(rep(1, P), amount = .25)), 
                                             "lambdaL1" = 2, 
                                             "lambdaL2" = 5, 
                                             "tau" = 1))
  
  out = run.jags(model = "jags_elastic_net.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
  file.remove("jags_elastic_net.txt")
  if (!is.null(cl)) {
    parallel::stopCluster(cl = cl)
  }
  return(out)
}


if (family == "binomial" || family == "logistic"){
  
  jags_elastic_net = "model{

              lambdaL1 ~ dgamma(0.25 , 0.2)
              lambdaL2 ~ dgamma(0.25 , 0.2)

              Intercept ~ dnorm(0, 1e-10)

              for (p in 1:P){
                eta[p] ~ dgamma(.5, (8 * lambdaL2 * sigma2) / pow(lambdaL1, 2)) T(1,)
                beta_prec[p] <- (lambdaL2/sigma2) * (eta[p]/(eta[p]-1))
                beta[p] ~ dnorm(0, beta_prec[p])
              }

              for (i in 1:N){
                logit(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])  + sum(design_beta[1:FP] * FX[i,1:FP])
                y[i] ~ dbern(psi[i])
                log_lik[i] <- logdensity.bern(y[i], psi[i])
                ySim[i] ~ dbern(psi[i])
              }
    
              Deviance <- -2 * sum(log_lik[1:N])
          }"
  
  P <- ncol(X)
  FP <- ncol(FX)
  write_lines(jags_elastic_net, "jags_elastic_net.txt")
  
  jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), sigma2 =pow(mean(y), -1) * pow(1 - mean(y), -1), FP = FP, FX = FX)
  
  monitor <- c("Intercept", "beta", "lambdaL1", "lambdaL2", "Deviance",  "ySim", "log_lik")
  
  if (log_lik == FALSE){
    monitor = monitor[-(length(monitor))]
  }
  
  inits <- lapply(1:chains, function(z) list("Intercept" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "binomial", lambda = 0.025, alpha = 0, standardize = FALSE))[1,1]), 
                                             "beta" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "binomial", lambda = 0.025, alpha = 0, standardize = FALSE))[-1,1]), 
                                             "eta" = 1 + abs(jitter(rep(1, P), amount = .25)), 
                                             "design_beta" = as.vector(coef(glmnet::glmnet(x = FX, y = y, family = "binomial", lambda = 0, alpha = 0, standardize = FALSE))[-1,1]),
                                             "lambdaL1" = 1, 
                                             "lambdaL2" = 1, 
                                             .RNG.name= "lecuyer::RngStream", 
                                             .RNG.seed= sample(1:10000, 1)))
  
  out = run.jags(model = "jags_elastic_net.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
  file.remove("jags_elastic_net.txt")
  if (!is.null(cl)) {
    parallel::stopCluster(cl = cl)
  }
  return(out)
}

if (family == "poisson"){
  
  jags_elastic_net = "model{

              lambdaL1 ~ dgamma(0.25 , 0.2)
              lambdaL2 ~ dgamma(0.25 , 0.2)

              Intercept ~ dnorm(0, 1e-10)

              for (p in 1:P){
                eta[p] ~ dgamma(.5, (8 * lambdaL2 * sigma2) / pow(lambdaL1, 2)) T(1,)
                beta_prec[p] <- (lambdaL2/sigma2) * (eta[p]/(eta[p]-1))
                beta[p] ~ dnorm(0, beta_prec[p])
              }

              for (i in 1:N){
                log(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P]) + sum(design_beta[1:FP] * FX[i,1:FP])
                y[i] ~ dpois(psi[i])
                log_lik[i] <- logdensity.pois(y[i], psi[i])
                ySim[i] ~ dpois(psi[i])
              }
    
              Deviance <- -2 * sum(log_lik[1:N])
          }"
  
  P <- ncol(X)
  FP <- ncol(FX)
  write_lines(jags_elastic_net, "jags_elastic_net.txt")
  
  jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), sigma2 = pow(mean(y), -1), FP = FP, FX = FX)
  
  monitor <- c("Intercept", "beta", "lambdaL1", "lambdaL2", "Deviance",  "ySim", "log_lik")
  if (log_lik == FALSE){
    monitor = monitor[-(length(monitor))]
  }
  inits <- lapply(1:chains, function(z) list("Intercept" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "poisson", lambda = 0.025, alpha = 0, standardize = FALSE))[1,1]), 
                                             "beta" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "poisson", lambda = 0.025, alpha = 0, standardize = FALSE))[-1,1]), 
                                             "design_beta" = as.vector(coef(glmnet::glmnet(x = FX, y = y, family = "poisson", lambda = 0, alpha = 0, standardize = FALSE))[-1,1]),
                                             "eta" = 1 + abs(jitter(rep(1, P), amount = .25)), 
                                             "lambdaL1" = 1, 
                                             "lambdaL2" = 1, 
                                             .RNG.name= "lecuyer::RngStream", 
                                             .RNG.seed= sample(1:10000, 1)))
  
  out = run.jags(model = "jags_elastic_net.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
  file.remove("jags_elastic_net.txt")
  if (!is.null(cl)) {
    parallel::stopCluster(cl = cl)
  }
  return(out)
}

}
