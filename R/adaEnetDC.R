#' Bayesian Adaptive Elastic Net with additional unpenalized design covariates
#'
#' @description This is an adaptation of the frequentist adaptive elastic net of Zou & Zhang (2009) to the Bayesian paradigm through a modification of the Bayesian elastic
#' net (Li & Lin, 2010). This function has the further allowance for a set of covariates that are not penalized. 
#' For example, you may wish to include variables such as age and gender so that  the coefficients for the other variables are 
#' penalized while controlling for these. This is a common need in research.
#' 
#' For the binomial and poisson likelihood functions the uniform-gamma scale mixture for the
#' variant of the Bayesian LASSO is adapted for use here. 
#' 
#' \cr
#' The model structure is given below: \cr
#' \cr
#' \cr
#' \if{html}{\figure{adaElasticNetDC.png}{}}
#' \if{latex}{\figure{adaElasticNetDC.png}{}}
#' \cr
#'
#' \cr
#' Plugin Pseudo-Variances: \cr
#' \if{html}{\figure{pseudovar.png}{}}
#' \if{latex}{\figure{pseudovar.png}{}}
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
#' @references 
#' Li, Qing; Lin, Nan. The Bayesian elastic net. Bayesian Anal. 5 (2010), no. 1, 151--170. doi:10.1214/10-BA506. https://projecteuclid.org/euclid.ba/1340369796 \cr
#' \cr
#' Zou, H.; Zhang, H. (2009) On the adaptive elastic-net with a diverging number of parameters, Ann. Statist. 37 , no. 4, 1733–1751, DOI 10.1214/08-AOS625. MR2533470 (2010j:62210) \cr
#' \cr
#' @return A run.jags object
#' @export
#'
#' @examples
#' adaEnetDC()
#'
adaEnetDC  = function(formula, design.formula, data, family = "gaussian", log_lik = FALSE, iter=10000, warmup=1000, adapt=2000, chains=4, thin=1, method = "parallel", cl = makeCluster(2), ...){
  
  X = as.matrix(model.matrix(formula, data)[,-1])
  y = model.frame(formula, data)[,1]
  FX <- as.matrix(model.matrix(design.formula, data)[, -1])
  
  if (family == "gaussian"){
  
  jags_adaptive_elastic_net = "model{

              tau ~ dgamma(.01, .01)
              sigma <- sqrt(1/tau)
              lambdaL2 ~ dgamma(0.50, 0.20)
              Intercept ~ dnorm(0, 1e-10)

              for (p in 1:P){
                lambdaL1[p] ~ dgamma(0.50, 0.20)
                eta[p] ~ dgamma(.5, (8 * lambdaL2 * pow(sigma,2)) / pow(lambdaL1[p], 2)) T(1,)
                beta_prec[p] <- (lambdaL2/pow(sigma,2)) * (eta[p]/(eta[p]-1))
                beta[p] ~ dnorm(0, beta_prec[p])
              }
              
              for (f in 1:FP){
                 design_beta[f] ~ dnorm(0, 1e-200)
              }

              for (i in 1:N){
                 y[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]) + sum(design_beta[1:FP] * FX[i,1:FP]), tau)
                 log_lik[i] <- logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]) + sum(design_beta[1:FP] * FX[i,1:FP]), tau)
                 ySim[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]) + sum(design_beta[1:FP] * FX[i,1:FP]), tau)
              }

              Deviance <- -2 * sum(log_lik[1:N])
          }"
  
  P <- ncol(X)
  FP <- ncol(FX)
  write_lines(jags_adaptive_elastic_net, "jags_adaptive_elastic_net.txt")
  jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), FP = FP, FX = FX)
  monitor <- c("Intercept", "beta", "design_beta", "sigma", "lambdaL1", "lambdaL2", "Deviance",  "ySim", "log_lik")
  if (log_lik == FALSE){
    monitor = monitor[-(length(monitor))]
  }
  inits <- lapply(1:chains, function(z) list("Intercept" = 0, 
                                             "beta" = rep(0, P), 
                                             "design_beta" = rep(0, FP),
                                             "eta" = 1 + abs(jitter(rep(1, P), amount = .25)), 
                                             "lambda1" = rep(1, P), 
                                             "lambda2" = 20, 
                                             "tau" = 1,
                                             .RNG.name= "lecuyer::RngStream", 
                                             .RNG.seed= sample(1:10000, 1)))
  
  out = run.jags(model = "jags_adaptive_elastic_net.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, 
                 data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, 
                 n.chains = chains, method = method, cl = cl, summarise = FALSE, ...)
  
  file.remove("jags_adaptive_elastic_net.txt")
  if (!is.null(cl)) {
    parallel::stopCluster(cl = cl)
  }
  return(out)
}


if (family == "binomial" || family == "logistic"){
  
  jags_adaptive_elastic_net = "model{

              
              lambdaL2 ~ dgamma(0.50, 0.20)

              Intercept ~ dnorm(0, 1e-10)

              for (p in 1:P){
                lambdaL1[p] ~ dgamma(0.50, 0.20)
                eta[p] ~ dgamma(.5, (8 * lambdaL2 * sigma2) / pow(lambdaL1[p], 2)) T(1,)
                beta_prec[p] <- (lambdaL2/sigma2) * (eta[p]/(eta[p]-1))
                beta[p] ~ dnorm(0, beta_prec[p])
              }
              
              for (f in 1:FP){
                 design_beta[f] ~ dnorm(0, 1e-200)
              }

              for (i in 1:N){
                logit(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P]) + sum(design_beta[1:FP] * FX[i,1:FP])
                y[i] ~ dbern(psi[i])
                log_lik[i] <- logdensity.bern(y[i], psi[i])
                ySim[i] ~ dbern(psi[i])
              }
    
              Deviance <- -2 * sum(log_lik[1:N])
          }"
  
  P <- ncol(X)
  
  write_lines(jags_adaptive_elastic_net, "jags_adaptive_elastic_net.txt")
  
  jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), sigma2 =pow(mean(y), -1) * pow(1 - mean(y), -1), FP = FP, FX = FX)
  
  monitor <- c("Intercept", "beta", "lambdaL1", "lambdaL2", "Deviance",  "ySim", "log_lik")
  
  if (log_lik == FALSE){
    monitor = monitor[-(length(monitor))]
  }
  
  inits <- lapply(1:chains, function(z) list("Intercept" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "binomial", lambda = 0.025, alpha = 0, standardize = FALSE))[1,1]), 
                                             "beta" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "binomial", lambda = 0.025, alpha = 0, standardize = FALSE))[-1,1]), 
                                             "eta" = 1 + abs(jitter(rep(1, P), amount = .25)), 
                                             "lambdaL1" = rep(1, P), 
                                             "lambdaL2" = 2, 
                                             "design_beta" = rep(0, FP),
                                             .RNG.name= "lecuyer::RngStream", 
                                             .RNG.seed= sample(1:10000, 1)))
  
  out = run.jags(model = "jags_adaptive_elastic_net.txt", modules = c("bugs on", "glm on", "dic off"),
                 monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, 
                 thin = thin, adapt = adapt, method = method,  n.chains = chains, cl = cl, summarise = FALSE, ...)
  
  file.remove("jags_adaptive_elastic_net.txt")
  if (!is.null(cl)) {
    parallel::stopCluster(cl = cl)
  }
  return(out)
}

if (family == "poisson"){
  
  jags_adaptive_elastic_net = "model{
    
              lambdaL2 ~ dgamma(0.50, 0.20)

              Intercept ~ dnorm(0, 1e-10)

              for (p in 1:P){
                lambdaL1[p] ~ dgamma(0.50, 0.20)
                eta[p] ~ dgamma(.5, (8 * lambdaL2 * sigma2) / pow(lambdaL1[p], 2)) T(1,)
                beta_prec[p] <- (lambdaL2/sigma2) * (eta[p]/(eta[p]-1))
                beta[p] ~ dnorm(0, beta_prec[p])
              }
              
              for (f in 1:FP){
                 design_beta[f] ~ dnorm(0, 1e-200)
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
  
  write_lines(jags_adaptive_elastic_net, "jags_adaptive_elastic_net.txt")
  
  jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), sigma2 = pow(mean(y), -1), FP = FP, FX = FX)
  
  monitor <- c("Intercept", "beta", "lambdaL1", "lambdaL2", "Deviance",  "ySim", "log_lik")
  
  if (log_lik == FALSE){
    monitor = monitor[-(length(monitor))]
  }
  
  inits <- lapply(1:chains, function(z) list("Intercept" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "poisson", lambda = 0.025, alpha = 0, standardize = FALSE))[1,1]), 
                                             "beta" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "poisson", lambda = 0.025, alpha = 0, standardize = FALSE))[-1,1]), 
                                             "eta" = 1 + abs(jitter(rep(1, P), amount = .25)), 
                                             "lambdaL1" = rep(1, P), 
                                             "lambdaL2" = 2, 
                                             "design_beta" = rep(0, FP),
                                             .RNG.name= "lecuyer::RngStream", 
                                             .RNG.seed= sample(1:10000, 1)))
  
  out = run.jags(model = "jags_adaptive_elastic_net.txt", modules = c("bugs on", "glm on", "dic off"), 
                 monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, 
                 thin = thin, adapt = adapt, method = method, cl = cl,  n.chains = chains, summarise = FALSE, ...)
  
  file.remove("jags_adaptive_elastic_net.txt")
  if (!is.null(cl)) {
    parallel::stopCluster(cl = cl)
  }
  return(out)
}
}

