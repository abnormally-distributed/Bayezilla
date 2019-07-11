#' Group+Within Group Selection with Bayesian Adaptive Elastic Net
#'
#' @description This is an adaptation of the frequentist adaptive elastic net of Ghosh (2007, 2011) and Zou & Zhang (2009) to the Bayesian paradigm through a modification of the Bayesian elastic
#' net (Li & Lin, 2010). It is further adapted such that coefficients can be assigned groups in the 
#' spirit of the Group LASSO (Yuan & Lin, 2006) and Group Bayesian LASSO (Kyung et al., 2010). Each group receives an 
#' independent L1 norm penalty, which is combined with the top level L2 penalty on a coefficient specific basis via
#' the truncated gamma priors. 
#'
#' \cr
#' The model structure is given below: \cr
#' \cr
#' \cr
#' \if{html}{\figure{groupadaptiveelasticNet.png}{}}
#' \if{latex}{\figure{groupadaptiveelasticNet.png}{}}
#' \cr
#'
#' @param X the model matrix. Construct this manually with model.matrix()[,-1]
#' @param y the outcome variable
#' @param idx the group labels. Should be of length = to ncol(model.matrix()[,-1]) with the group assignments
#' for each covariate. Please ensure that you start numbering with 1, and not 0.
#' @param family one of "gaussian" (default), "binomial", or "poisson"
#' @param log_lik Should the log likelihood be monitored? The default is FALSE.
#' @param iter How many post-warmup samples? Defaults to 10000.
#' @param warmup How many warmup samples? Defaults to 5000.
#' @param adapt How many adaptation steps? Defaults to 5000.
#' @param chains How many chains? Defaults to 4.
#' @param thin Thinning interval. Defaults to 3.
#' @param method Defaults to "parallel". For an alternative parallel option, choose "rjparallel" or. Otherwise, "rjags" (single core run).
#' @param cl Use parallel::makeCluster(# clusters) to specify clusters for the parallel methods. Defaults to two cores.
#' @param ... Other arguments to run.jags.
#'
#' @references 
#' Ghosh, S. (2007) Adaptive Elastic Net: A Doubly Regularized method for variable selection to Achieve Oracle Properties. Tech. Rep. pr07-01, available at http://www.math.iupui.edu/research/preprints.php, IUPUI  \cr
#' \cr
#' Ghosh, S. (2011) On the grouped selection and model complexity of the adaptive elastic net. Statistics and Computing 21, no. 3, 451. https://doi.org/10.1007/s11222-010-9181-4 \cr
#' \cr
#' Kyung, M., Gill, J., Ghosh, M., and Casella, G. (2010). Penalized regression, standard errors, and bayesian lassos. Bayesian Analysis, 5(2):369–411. \cr
#' \cr
#' Li, Qing; Lin, Nan. (2010) The Bayesian elastic net. Bayesian Anal. 5, no. 1, 151--170. doi:10.1214/10-BA506. https://projecteuclid.org/euclid.ba/1340369796 \cr
#' \cr
#' Yuan, Ming; Lin, Yi (2006) Model Selection and Estimation in Regression with Grouped Variables. Journal of the Royal Statistical Society. Series B (statistical Methodology). Wiley. 68 (1): 49–67. doi:10.1111/j.1467-9868.2005.00532.x \cr
#' \cr
#' Zou, H.; Zhang, H. (2009) On the adaptive elastic-net with a diverging number of parameters, Ann. Statist. 37 , no. 4, 1733–1751, DOI 10.1214/08-AOS625. MR2533470 (2010j:62210) \cr
#' \cr
#' @return A run.jags object
#' @export
#'
#' @examples
#' groupAdaEnet()
#'
groupAdaEnet  = function(X, y, idx, family = "gaussian", log_lik = FALSE, iter=10000, warmup=5000, adapt=5000, chains=4, thin=3, method = "parallel", cl = makeCluster(2), ...){

  if (family == "gaussian"){
  
  jags_adaptive_elastic_net = "model{

              tau ~ dgamma(.01, .01)
              sigma <- sqrt(1/tau)
              
              Intercept ~ dnorm(0, 1e-10)
              lambdaL2 ~ dgamma(1.5, 0.25)
              
              for (g in 1:nG){
                  lambdaL1[g] ~ dgamma(1.5, 0.25)
              }
              
              for (p in 1:P){
          
                eta[p] ~ dgamma(k[idx[p]] *.5, (8 * lambdaL2 * pow(sigma,2)) / pow(lambdaL1[idx[p]], 2)) T(1,)
                beta_prec[p] <- (lambdaL2/pow(sigma,2)) * (eta[p]/(eta[p]-1))
                beta[p] ~ dnorm(0, beta_prec[p])
              }

              for (i in 1:N){
                 y[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 log_lik[i] <- logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 ySim[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
              }

              Deviance <- -2 * sum(log_lik[1:N])
          }"
  
  P <- ncol(X)
  write_lines(jags_adaptive_elastic_net, "jags_adaptive_elastic_net.txt")
  jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), idx = idx, nG = max(idx), k = as.vector(table(idx)))
  monitor <- c("Intercept", "beta", "sigma", "lambdaL1", "lambdaL2", "Deviance",  "ySim", "log_lik")
  
  if (log_lik == FALSE){
    monitor = monitor[-(length(monitor))]
  }
  
  inits <- lapply(1:chains, function(z) list("Intercept" = lmSolve(y ~ ., data = data.frame(y = y, X))[1], 
                                             "beta" = lmSolve(y ~ ., data = data.frame(y = y, X))[-1], 
                                             "eta" = 1 + abs(jitter(rep(1, P), amount = .25)), 
                                             "lambdaL1" = rep(2, max(idx)), 
                                             "lambdaL2" = 3,  
                                             "tau" = 1,
                                             .RNG.name = "lecuyer::RngStream", 
                                             .RNG.seed = sample(1:1000, 1) 
                                             ))
  
  out = run.jags(model = "jags_adaptive_elastic_net.txt", 
                 modules = c("bugs on", "glm on", "dic off"), 
                 monitor = monitor, 
                 data = jagsdata, 
                 inits = inits, 
                 burnin = warmup, 
                 sample = iter, 
                 thin = thin, 
                 adapt = adapt, 
                 method = method, 
                 n.chains = chains, 
                 cl = cl, 
                 summarise = FALSE, 
                 ...)
  file.remove("jags_adaptive_elastic_net.txt")
  if (!is.null(cl)) {
    parallel::stopCluster(cl = cl)
  }
  return(out)
  
  }

  if (family == "binomial"){
    
    jags_elastic_net = "model{
              
              Intercept ~ dnorm(0, 1e-10)
              lambdaL2 ~ dgamma(0.25 , 0.20)
              
              for (g in 1:nG){
                  lambdaL1[g] ~ dgamma(0.25 , 0.20)
              }
              
              for (p in 1:P){
                eta[p] ~ dgamma(k[idx[p]] *.5, (8 * lambdaL2 * sigma2) / pow(lambdaL1[idx[p]], 2)) T(1,)
                beta_prec[p] <- (lambdaL2/sigma2) * (eta[p]/(eta[p]-1))
                beta[p] ~ dnorm(0, beta_prec[p])
              }
              
              for (i in 1:N){
                logit(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                y[i] ~ dbern(psi[i])
                log_lik[i] <- logdensity.bern(y[i], psi[i])
                ySim[i] ~ dbern(psi[i])
              }

              Deviance <- -2 * sum(log_lik[1:N])
          }"
    
    P <- ncol(X)
    write_lines(jags_elastic_net, "jags_elastic_net.txt")
    jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), idx = idx, nG = max(idx), k = as.vector(table(idx)), sigma2 = pow(mean(y), -1) * pow(1 - mean(y), -1))
    monitor <- c("Intercept", "beta", "lambdaL1", "lambdaL2", "Deviance",  "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits <- lapply(1:chains, function(z) list("Intercept" = 0, 
                                               "beta" = rep(0, P), 
                                               "eta" = 1 + abs(jitter(rep(1, max(idx)), amount = .25)), 
                                               "lambdaL1" = rep(2, max(idx)), 
                                               "lambdaL2" = 5, 
                                               .RNG.name= "lecuyer::RngStream", 
                                               .RNG.seed= sample(1:10000, 1)))
    
    out = run.jags(model = "jags_elastic_net.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, n.chains = chains, method = method, cl = cl, summarise = FALSE, ...)
    file.remove("jags_elastic_net.txt")
    if (!is.null(cl)) {
      parallel::stopCluster(cl = cl)
    }
    return(out)
  }
  
  
  if (family == "poisson"){
    
    jags_elastic_net = "model{

              Intercept ~ dnorm(0, 1e-10)
              lambdaL2 ~ dgamma(0.25 , 0.20)
              
              for (g in 1:nG){
                  lambdaL1[g] ~ dgamma(0.25 , 0.20)
              }
              
              for (p in 1:P){
                eta[p] ~ dgamma(k[idx[p]] *.5, (8 * lambdaL2 * sigma2) / pow(lambdaL1[idx[p]], 2)) T(1,)
                beta_prec[p] <- (lambdaL2/sigma2) * (eta[p]/(eta[p]-1))
                beta[p] ~ dnorm(0, beta_prec[p])
              }
              
              for (i in 1:N){
                log(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                y[i] ~ dpois(psi[i])
                log_lik[i] <- logdensity.pois(y[i], psi[i])
                ySim[i] ~ dpois(psi[i])
              }

              Deviance <- -2 * sum(log_lik[1:N])
          }"
    
    P <- ncol(X)
    write_lines(jags_elastic_net, "jags_elastic_net.txt")
    jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), idx = idx, nG = max(idx), k = as.vector(table(idx)), sigma2 = pow(mean(y), -1))
    monitor <- c("Intercept", "beta", "lambdaL1", "lambdaL2", "Deviance",  "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits <- lapply(1:chains, function(z) list("Intercept" = 0, 
                                               "beta" = rep(0, P), 
                                               "eta" = 1 + abs(jitter(rep(1, max(idx)), amount = .25)), 
                                               "lambdaL1" = rep(2, max(idx)), 
                                               "lambdaL2" = 5, 
                                               .RNG.name= "lecuyer::RngStream", 
                                               .RNG.seed= sample(1:10000, 1)))
    
    out = run.jags(model = "jags_elastic_net.txt", modules = c("bugs on", "glm on", "dic off"), 
                   monitor = monitor, data = jagsdata, 
                   inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt,
                   n.chains = chains, method = method, cl = cl, summarise = FALSE, ...)
    file.remove("jags_elastic_net.txt")
    if (!is.null(cl)) {
      parallel::stopCluster(cl = cl)
    }
    return(out)
  }
  
}
