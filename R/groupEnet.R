#' Group Elastic Net for Gaussian Likelihood
#'
#' @description The Bayesian elastic net described by Li and Lin (2010) modified for use as a group selection model,
#' akin to the Group Bayesian LASSO described by Kyung et al. (2010). Group selection is a method
#' described first by Yuan & Lin (2006) for applying shrinkage penalties to coefficients that have some natural grouping,
#' such as refelcting dummy variables of a single factor, or coefficients corresponding to related variables (for example, predictors
#' derived from a single brain region). Note only the Gaussian likelihood is 
#' provided because the Bayesian elastic net requires conditioning on the error variance, which GLM-families
#' do not have.
#' 
#' \cr
#' The model structure is given below: \cr
#' \cr
#' \cr
#' \if{html}{\figure{groupelasticNet.png}{}}
#' \if{latex}{\figure{groupelasticNet.png}{}}
#' \cr
#'
#' @param X the model matrix. Construct this manually with model.matrix()[,-1]
#' @param y the outcome variable
#' @param idx the group labels. Should be of length = to ncol(model.matrix()[,-1]) with the group assignments
#' for each covariate. Please ensure that you start numbering with 1, and not 0.
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
#' Yuan, Ming; Lin, Yi (2006). Model Selection and Estimation in Regression with Grouped Variables. Journal of the Royal Statistical Society. Series B (statistical Methodology). Wiley. 68 (1): 49–67. doi:10.1111/j.1467-9868.2005.00532.x \cr
#' \cr
#' Kyung, M., Gill, J., Ghosh, M., and Casella, G. (2010). Penalized regression, standard errors, and bayesian lassos. Bayesian Analysis, 5(2):369–411. \cr
#' \cr
#' Li, Qing; Lin, Nan. The Bayesian elastic net. Bayesian Anal. 5 (2010), no. 1, 151--170. doi:10.1214/10-BA506. https://projecteuclid.org/euclid.ba/1340369796
#' \cr
#' @return A run.jags object
#' @export
#'
#' @examples
#' groupEnet()
#'
groupEnet  = function(X, y, idx, log_lik = FALSE, iter=10000, warmup=1000, adapt=2000, chains=4, thin=1, method = "parallel", cl = makeCluster(2), ...){
  
  jags_elastic_net = "model{

              tau ~ dgamma(.01, .01)
              sigma <- sqrt(1/tau)
              lambdaL1 ~ dgamma(.25, .01)
              lambdaL2 ~ dgamma(.25, .01)

              Intercept ~ dnorm(0, 1)

              for (g in 1:nG){
                eta[g] ~ dgamma(k[g] * .5, (8 * lambdaL2 * pow(sigma,2)) / pow(lambdaL1, 2)) T(1,)
                beta_prec[g] <- (lambdaL2/pow(sigma,2)) * (eta[g]/(eta[g]-1))
              }
              
              for (p in 1:P){
                beta[p] ~ dnorm(0, beta_prec[idx[p]])
              }
              
              for (i in 1:N){
                 y[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 log_lik[i] <- logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 ySim[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
              }

              Deviance <- -2 * sum(log_lik[1:N])
          }"
  
  P <- ncol(X)
  write_lines(jags_elastic_net, "jags_elastic_net.txt")
  jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), idx = idx, nG = max(idx), k = as.vector(table(idx)))
  monitor <- c("Intercept", "beta", "sigma", "lambdaL1", "lambdaL2", "Deviance", "eta", "ySim", "log_lik")
  if (log_lik == FALSE){
    monitor = monitor[-(length(monitor))]
  }
  inits <- lapply(1:chains, function(z) list("Intercept" = 0, 
                                             "beta" = rep(0, P), 
                                             "eta" = 1 + abs(jitter(rep(1, max(idx)), amount = .25)), 
                                             "lambdaL1" = 20, 
                                             "lambdaL2" = 5, 
                                             "tau" = 1,
                                             .RNG.name= "lecuyer::RngStream", 
                                             .RNG.seed= sample(1:10000, 1)))
  
  out = run.jags(model = "jags_elastic_net.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
  file.remove("jags_elastic_net.txt")
  if (!is.null(cl)) {
    parallel::stopCluster(cl = cl)
  }
  return(out)
}
