#' Bayesian OLS for Heteroskedastic Residuals
#'
#' @description Just a modification of the \code{\link[Bayezilla]{glmFlat}} 
#' function to model heteroskedasticity. If the variance of the response variable
#' increases or decreases as a function of the predictor variable(s) then use this
#' function.
#' 
#' \cr
#' Model Specification:
#' \cr
#' \if{html}{\figure{hetLm.png}{}}
#' \if{latex}{\figure{hetLm.png}{}}
#' 
#' @param formula the model formula
#' @param data a data frame
#' @param log_lik Should the log likelihood be monitored? The default is FALSE.
#' @param iter How many post-warmup samples? Defaults to 10000.
#' @param warmup How many warmup samples? Defaults to 1000.
#' @param adapt How many adaptation steps? Defaults to 1000.
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
#' hetLm()
#'
hetLm  = function(formula, data, log_lik = FALSE, iter=10000, warmup=1000, adapt=1000, chains=4, thin=1, method = "parallel", cl = makeCluster(2), ...){
  
  X = as.matrix(model.matrix(formula, data)[,-1])
  y = model.frame(formula, data)[,1]

    jags_glm = "model{
              for (p in 1:P){
                beta[p] ~ dnorm(0, 1e-200)
              }
              
              Intercept ~ dnorm(0, 1e-200)
              sigmaIntercept ~ dnorm(0, 1e-200)
              
              for (i in 1:N){
                 mu[i] <- Intercept + sum(beta[1:P] * X[i,1:P])
                 log(sigma_hat[i]) <- sigmaIntercept + sum(beta[1:P] * X[i,1:P])
                 tau[i] <- 1 / pow(sigma_hat[i], 2)
                 y[i] ~ dnorm(mu[i], tau[i])
                 log_lik[i] <- logdensity.norm(y[i], mu[i], tau[i])
                 ySim[i] ~ dnorm(mu[i], tau[i])
              }
              
              Deviance <- -2 * sum(log_lik[1:N])
          }"
    
    P = ncol(X)
    write_lines(jags_glm, "jags_glm.txt")
    jagsdata = list(X = X, y = y,  N = length(y), P = ncol(X))
    monitor = c("Intercept", "beta", "sigma_hat", "Deviance", "ySim" ,"log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = as.vector(coef(lm(formula, data)))[1],
                                              "sigmaIntercept" = 0, 
                                              "beta" = as.vector(coef(lm(formula, data)))[-1], 
                                              "ySim" = y, 
                                              .RNG.name= "lecuyer::RngStream", 
                                              .RNG.seed = sample(1:20000, 1)))
    out = run.jags(model = "jags_glm.txt", modules = c("glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE,...)
    if (is.null(cl) == FALSE){
      parallel::stopCluster(cl = cl)
    }
    file.remove("jags_glm.txt")
    return(out)
}
