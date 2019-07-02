#' Conditional Bias Corrected Re-Estimation of Selected Variables
#'
#' @description At times the bias induced by continuous shrinkage prior methods such as the Bayesian LASSO may be too great to be
#' ignored and gains in predictive and inferential accuracy can be made by performing conditional re-estimation.  The simplest 
#' form of this is the OLS-LASSO / MLE-LASSO hybrid of Efron et al (2004), where the LASSO-selected variables are re-estimated 
#' using maximum likelihood. However, this ignores selection effects. Alternatives are the Relaxed LASSO procedure of Meinshausen N (2007)
#' which performs a second LASSO on the variables selected by the first LASSO. The second pass of the LASSO will select a smaller 
#' penalty due to the reduced variable set, hence providing estimates closer to the MLEs while retaining some small shrinkage to attenuate
#' optimistic selection effects. A general frequentist re-estimation procedure has been put forth by Taylor & Tibsharini (2015) that
#' is general enough to be used with various variable selection methods. Similarly, Sillanpää &  Mutshinda (2011) provided a means
#' of conditional re-estimation in their discussion of post-selection inference using the shape adaptive shrinkage prior (offered here
#' in the \code{\link[Bayezilla]{sasp}} function). The benefits of 
#' such two-stage estimation procedures are discussed by Belloni & Chernozhukov (2013). Belloni & Chernozhukov demonstrate that such
#' re-estimation dominates the naive MLE-LASSO (or MLE-your-favorite-variable-selection-method-here) procedure. \cr
#' 
#' \cr
#' The method of Sillanpää &  Mutshinda is adapted here. Their conditional re-estimation procedure involves supplying the full model
#' matrix that was passed to the LASSO or other selection procedure, along with a vector of indicator variables that take values of [0, 1]
#' to indicate the elimination or inclusion of the variable. They then estimate the raw coefficients using vague normal priors and multiply
#' these by the indicator variables. Here cauchy priors are used. This is done to account for the eliminated variables for purposes of avoiding selection bias, while
#' also only obtaining non-zero estimates for the selected variables.  
#' 
#' \cr
#' 
#' This procedure can be used following any variable selection process, Bayesian or not. All you need to do is provide the full model
#' formula and a vector of inclusion indicators to utilize this. This function assumes the variables are 
#' standardized just as the selection models require. \cr 
#' 
#' Standard gaussian, binomial, and poisson likelihood functions are available. \cr
#' \cr
#' Model Specification:\cr
#' \cr
#' \if{html}{\figure{cR.png}{}}
#' \if{latex}{\figure{cR.png}{}}
#' 
#' @param formula the model formula
#' @param data a data frame
#' @param delta a vector of inclusion indicators. Must only contain values of 0 and 1. 
#' @param family one of "gaussian", "binomial", or "poisson".
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
#' @references 
#' 
#' Belloni, A.; Chernozhukov, V. (2013) Least squares after model selection in high-dimensional sparse models. Bernoulli 19, no. 2, 521--547. doi:10.3150/11-BEJ410.  \cr
#' \cr
#' Efron B., Hastie T., Johnstone I., Tibshirani, R. (2004). Least angle regression. Ann Stat 32: 407–499. \cr
#' \cr
#' Meinshausen N (2007). Relaxed Lasso. Comp Stat Data Anal 52: 374–393. \cr
#' \cr
#' Sillanpää, S., &  Mutshinda, C., (2011) Bayesian shrinkage analysis of QTLs under shape-adaptive shrinkage priors, and accurate re-estimation of genetic effects. Heredity volume 107, pages 405–412. doi: 10.1038/hdy.2011.37 \cr
#' \cr
#' Taylor, J., & Tibshirani, R. (2015) Statistical learning and selective inference. Proceedings of the National Academy of Sciences, 112 (25) 7629-7634; doi: 10.1073/pnas.1507583112 \cr
#' \cr
#' 
#' @return A run.jags object
#' @export
#'
#' @examples
#' cR(formula, data, delta = c(0, 0, 1, 1, 1, 0, 1, 0))

cR  = function(formula, data, delta = NULL, family = "gaussian",log_lik = FALSE, iter=10000, warmup=1000, adapt=1000, chains=4, thin=1, method = "parallel", cl = makeCluster(2), ...){
  
  if (is.null(delta)) {
    stop("Please provide a vector of indicator variables")
  }

  X = as.matrix(model.matrix(formula, data)[,-1])
  y = model.frame(formula, data)[,1]
  
  
  if (family == "gaussian"){
    
    jags_glm = "model{
    
              tau ~ dgamma(0.01, 0.01)
              sigma <- sqrt(1/tau)
    
              for (j in 1:P){
                theta[j] ~ dt(0, 1, 1)
                beta[j] <- theta[j]*delta[j]
              }
              
              Intercept ~ dnorm(0, 1e-200)

              for (i in 1:N){
                 y[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 log_lik[i] <- logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 ySim[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
              }
              
              Deviance <- -2 * sum(log_lik[1:N])
          }"
    
    P = ncol(X)
    write_lines(jags_glm, "jags_glm.txt")
    jagsdata = list(X = X, y = y,  N = length(y), P = ncol(X), delta = delta, beta_hat= beta_hat)
    monitor = c("Intercept", "beta", "sigma", "Deviance", "ySim" ,"log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "gaussian", lambda = 0.025, alpha = 0, standardize = FALSE))[1,1]), 
                                              "theta" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "gaussian", lambda = 0.025, alpha = 0, standardize = FALSE))[-1,1]), 
                                              "tau" = 1, 
                                              "ySim" = y, 
                                              .RNG.name= "lecuyer::RngStream", 
                                              .RNG.seed = sample(1:20000, 1)))
  }
  
  if (family == "binomial" || family == "logistic"){
    
    jags_glm = "model{
      
              for (p in 1:P){
                theta[p] ~ dt(0, 1, 1)
                beta[p] <- theta[p]*delta[p]
              }

              Intercept ~ dnorm(0, 1e-200)

              for (i in 1:N){
                 logit(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                 y[i] ~ dbern(psi[i])
                 log_lik[i] <- logdensity.bern(y[i], psi[i])
                 ySim[i] ~ dbern(psi[i])
              }
             Deviance <- -2 * sum(log_lik[1:N])
          }"
    
    P = ncol(X)
    write_lines(jags_glm, "jags_glm.txt")
    jagsdata = list(X = X, y = y, N = length(y),  P = ncol(X), delta = delta)
    monitor = c("Intercept", "beta", "Deviance", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "binomial", lambda = 0.025, alpha = 0, standardize = FALSE))[1,1]), 
                                              "theta" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "binomial", lambda = 0.025, alpha = 0, standardize = FALSE))[-1,1]), 
                                              "ySim" = y, 
                                              .RNG.name= "lecuyer::RngStream", 
                                              .RNG.seed= sample(1:20000, 1)))
  }
  
  if (family == "poisson"){
    
    jags_glm = "model{

              for (p in 1:P){
                theta[p] ~ dt(0, 1, 1)
                beta[p] <- theta[p]*delta[p]
              }

              Intercept ~ dnorm(0, 1e-200)

              for (i in 1:N){
                 log(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                 y[i] ~ dpois(psi[i])
                 log_lik[i] <- logdensity.pois(y[i], psi[i])
                 ySim[i] ~ dpois(psi[i])
              }

              Deviance <- -2 * sum(log_lik[1:N])
          }"
    
    write_lines(jags_glm, "jags_glm.txt")
    P = ncol(X)
    jagsdata = list(X = X, y = y, N = length(y),  P = ncol(X), delta = delta)
    monitor = c("Intercept", "beta", "Deviance", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "poisson", lambda = 0.025, alpha = 0, standardize = FALSE))[1,1]), 
                                              "theta" = as.vector(coef(glmnet::glmnet(x = X, y = y, family = "poisson", lambda = 0.025, alpha = 0, standardize = FALSE))[-1,1]), 
                                              "ySim" = y,
                                              .RNG.name= "lecuyer::RngStream", 
                                              .RNG.seed= sample(1:20000, 1)))
  }
  
  out = run.jags(model = "jags_glm.txt", modules = c("glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE,...)
  if (is.null(cl) == FALSE){
    parallel::stopCluster(cl = cl)
  }
  file.remove("jags_glm.txt")
  return(out)
}
