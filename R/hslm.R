#' Bayesian Ordinary Linear Regression with Error Variance Modeling 
#'
#' @description A modification of the \code{\link[Bayezilla]{glmBayes}} function to model 
#' heteroskedasticity as a function of all or some variables. 
#' 
#' \cr
#' Model Specification: \cr
#' \cr
#' \if{html}{\figure{hslm.png}{}}
#' \if{latex}{\figure{hslm.png}{}}
#' 
#' @param formula the model formula
#' @param sigma.formula the formula for the scale parameter. Defaults to ~ 1 (intercept only)
#' @param data a data frame
#' @param df degrees of freedom for prior.
#' @param s The desired prior scale. Defaults to 1. Is automatically squared within the model so
#' select a number here on the standard deviation scale.
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
#' hslm()
#'
hslm  = function(formula, sigma.formula = ~ 1, data, s = 1, df = 1, log_lik = FALSE, iter=10000, warmup=1000, adapt=2000, chains=4, thin=3, method = "parallel", cl = makeCluster(2), ...){
  
    jags_hslm = "model{
              
              Intercept ~ dnorm(0, 1e-10)
              eta ~ dgamma(df * 0.50, pow(s, 2) * (df * 0.50))
              
              for (p in 1:P){
                beta[p] ~ dnorm(0, eta)
              }
                            
              for (s in 1:S){
                sigma_beta[s] ~ dnorm(0, 1)
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
    write_lines(jags_hslm, "jags_hslm.txt")
    jagsdata = list(X = X, SX = SX, y = y,  N = length(y), P = P, S = S, s = s, df = df)
    monitor = c("Intercept", "beta", "sigma_beta", "deviance", "ySim" ,"log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = 0, 
                                              "beta" = jitter(rep(0, P), amount = 1), 
                                              "sigma_beta" = jitter(rep(0, S), amount = 1), 
                                              "ySim" = sample(y, length(y)), 
                                              "eta" = 1, 
                                              .RNG.name= "lecuyer::RngStream", 
                                              .RNG.seed = sample(1:10000, 1)))
  
  out = run.jags(model = "jags_hslm.txt", modules = c("glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, n.chains = chains, adapt = adapt, method = method, cl = cl, summarise = FALSE,...)
  if (is.null(cl) == FALSE){
    parallel::stopCluster(cl = cl)
  }
  file.remove("jags_hslm.txt")
  return(out)
}
