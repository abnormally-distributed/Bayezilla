#' Bayesian Elastic Net for Gaussian Likelihood
#'
#' @param formula the model formula
#' @param data a data frame.
#' @param log_lik Should the log likelihood be monitored? The default is FALSE.
#' @param iter How many post-warmup samples? Defaults to 10000.
#' @param warmup How many warmup samples? Defaults to 1000.
#' @param adapt How many adaptation steps? Defaults to 2000.
#' @param chains How many chains? Defaults to 4.
#' @param thin Thinning interval. Defaults to 3.
#' @param method Defaults to "rjags" (single core run). For parallel, choose "rjparallel" or "parallel".
#' @param cl Use parallel::makeCluster(# clusters) to specify clusters for the parallel methods.
#' @param ... Other arguments to run.jags.
#'
#' @return A run.jags object
#' @export
#'
#' @examples
#' bayesEnet()
#'
bayesEnet  = function(formula, data, family = "gaussian", log_lik = FALSE, iter=10000, warmup=1000, adapt=2000, chains=4, thin=3, method = "rjags", cl = NULL, ...){

  X = model.matrix(formula, data)[,-1]
  y = model.frame(formula, data)[,1]

    jags_elastic_net = "model{

              tau ~ dgamma(.001, .001)
              sigma <- sqrt(1/tau)
              lambda1 ~ dunif(1e-6, 500)
              lambda2 ~ dunif(1e-6, 500)

              Intercept ~ dnorm(0, .01)

              for (p in 1:P){
                omega[p] ~ dgamma(.5, (8 * lambda2 * sigma^2) / lambda1^2) T(1,)
                beta_prec[p] <- (lambda2/sigma^2) * (omega[p]/(omega[p]-1))
                beta[p] ~ dnorm(0, beta_prec[p])
              }

              for (i in 1:N){
                 y[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 log_lik[i] <- logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 ySim[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
              }

              deviance <- -2 * sum(log_lik[1:N])
          }"

    P <- ncol(X)
    write_lines(jags_elastic_net, "jags_elastic_net.txt")
    jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X))
    monitor <- c("Intercept", "beta", "sigma", "lambda1", "lambda2", "deviance", "omega", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits <- lapply(1:chains, function(z) list("Intercept" = 0, "beta" = rep(0, P), "omega" = 1 + abs(jitter(rep(1, P), amount = .25)), "lambda1" = 50, "lambda2" = 15, "tau" = 1))
  out = run.jags(model = "jags_elastic_net.txt", modules = "glm", monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, ...)
  file.remove("jags_elastic_net.txt")
  return(out)
}
