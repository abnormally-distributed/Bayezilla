#' Bayesian Elastic Net for Gaussian Likelihood
#'
#' @description The Bayesian elastic net described by Li and Lin (2010). Note only the Gaussian likelihood is 
#' provided because the Bayesian elastic net requires conditioning on the error variance, which GLM-families
#' do not have.
#' 
#' \cr
#' The model structure is given below: \cr
#' \cr
#' \cr
#' \if{html}{\figure{elasticNet.png}{}}
#' \if{latex}{\figure{elasticNet.png}{}}
#' \cr
#'
#' @param formula the model formula
#' @param data a data frame.
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
#' bayesEnet()
#'
bayesEnet  = function(formula, data, log_lik = FALSE, iter=10000, warmup=1000, adapt=2000, chains=4, thin=1, method = "parallel", cl = makeCluster(2), ...){

  X = model.matrix(formula, data)[,-1]
  y = model.frame(formula, data)[,1]

    jags_elastic_net = "model{

              tau ~ dgamma(.01, .01)
              sigma <- sqrt(1/tau)
              lambda1 ~ dgamma(0.5 , 0.001)
              lambda2 ~ dgamma(0.5 , 0.001)

              Intercept ~ dnorm(0, 1)

              for (p in 1:P){
                eta[p] ~ dgamma(.5, (8 * lambda2 * pow(sigma,2)) / pow(lambda1, 2)) T(1,)
                beta_prec[p] <- (lambda2/pow(sigma,2)) * (eta[p]/(eta[p]-1))
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
    write_lines(jags_elastic_net, "jags_elastic_net.txt")
    jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X))
    monitor <- c("Intercept", "beta", "sigma", "lambda1", "lambda2", "Deviance", "eta", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
  inits <- lapply(1:chains, function(z) list("Intercept" = 0, .RNG.name= "lecuyer::RngStream", .RNG.seed= sample(1:10000, 1), "beta" = rep(0, P), "eta" = 1 + abs(jitter(rep(1, P), amount = .25)), "lambda1" = 20, "lambda2" = 5, "tau" = 1))

  out = run.jags(model = "jags_elastic_net.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
  file.remove("jags_elastic_net.txt")
  if (!is.null(cl)) {
    parallel::stopCluster(cl = cl)
  }
  return(out)
}
