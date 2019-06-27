#' Group+Within Group Selection with Bayesian Adaptive Elastic Net
#'
#' @description This is an adaptation of the frequentist adaptive elastic net of Zou & Zhang (2009) to the Bayesian paradigm through a modification of the Bayesian elastic
#' net (Li & Lin, 2010). It is further adapted such that coefficients can be assigned groups. Each group receives an independent
#' L2 norm penalty, followed by coefficient specific L1 norm penalities. 
#' 
#' Note only the Gaussian likelihood is 
#' provided because the Bayesian elastic net requires conditioning on the error variance, which GLM-families
#' do not have.
#' 
#' \cr
#' The model structure is given below: \cr
#' \cr
#' \cr
#' \if{html}{\figure{adaElasticNet.png}{}}
#' \if{latex}{\figure{adaElasticNet.png}{}}
#' \cr
#'
#' @param X the model matrix. Construct this manually with model.matrix()[,-1]
#' @param y the outcome variable
#' @param idx the group labels. Should be of length = to ncol(model.matrix()[,-1]) with the group assignments
#' for each covariate. Please ensure that you start numbering with 1, and not 0.
#' @param log_lik Should the log likelihood be monitored? The default is FALSE.
#' @param iter How many post-warmup samples? Defaults to 10000.
#' @param warmup How many warmup samples? Defaults to 2000.
#' @param adapt How many adaptation steps? Defaults to 5000.
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
#' groupAdaEnet()
#'
groupAdaEnet  = function(X, y, idx, log_lik = FALSE, iter=10000, warmup=2000, adapt=5000, chains=4, thin=1, method = "parallel", cl = makeCluster(2), ...){

  jags_adaptive_elastic_net = "model{

              tau ~ dgamma(.01, .01)
              sigma <- sqrt(1/tau)
              
              Intercept ~ dnorm(0, 1)

              for (g in 1:nG){
                  lambdaL2[g] ~ dgamma(0.5 , 0.001)
              }
              
              for (p in 1:P){
                lambdaL1[p] ~ dgamma(0.5 , 0.001)
                eta[p] ~ dgamma(.5, (8 * lambdaL2[idx[p]] * pow(sigma,2)) / pow(lambdaL1[p], 2)) T(1,)
                beta_prec[p] <- (lambdaL2[idx[p]]/pow(sigma,2)) * (eta[p]/(eta[p]-1))
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
  jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), idx = idx, nG = max(idx))
  monitor <- c("Intercept", "beta", "sigma", "lambdaL1", "lambdaL2", "Deviance", "eta", "ySim", "log_lik")
  if (log_lik == FALSE){
    monitor = monitor[-(length(monitor))]
  }
  inits <- lapply(1:chains, function(z) list("Intercept" = 0, 
                                             "beta" = rep(0, P), 
                                             "eta" = 1 + abs(jitter(rep(1, P), amount = .25)), 
                                             "lambdaL1" = rep(10, P), 
                                             "lambdaL2" = rep(20, max(idx),  
                                             "tau" = 1,
                                             .RNG.name= "lecuyer::RngStream", 
                                             .RNG.seed= sample(1:10000, 1) 
                                             )))
  
  out = run.jags(model = "jags_adaptive_elastic_net.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
  file.remove("jags_adaptive_elastic_net.txt")
  if (!is.null(cl)) {
    parallel::stopCluster(cl = cl)
  }
  return(out)
}
