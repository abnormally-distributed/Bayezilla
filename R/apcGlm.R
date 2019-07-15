#' Adaptive Powered Correlation Prior 
#'
#' This function implements the adaptive powered correlation prior for estimating a single general(ized) linear regression
#' model. 
#' \cr
#' The adaptive powered correlation prior extends the Zellner-Siow Cauchy g-prior by allowing the crossproduct of the 
#' model matrix to be raised to powers other than -1 (which gives the Fisher information matrix). The power here will 
#' be referred to as "lambda". A lambda of 0 results in an identity matrix, which results in a ridge-regression like
#' prior. Positive values of lambda adapt to collinearity by allowing correlated predictors to enter and exit the model 
#' together. Negative values of lambda on the other hand favor including only one of a set of correlated predictors. 
#' This can be understood as projecting the information matrix into a new space which leads to a model
#' similar in function to principal components regression (Krishna et al., 2009). In this implementation full Bayesian
#' inference is used for lambda, rather than searching via marginal likelihood maximization as Krishna et al. (2009) did. 
#' The reason for this is twofold. First, full Bayesian inference means the model has to be fit only once instead of
#' several times over a grid of candidate values for lambda. Second, this avoids any coherency problems such as those
#' that arise when using fixed-g priors.
#' \cr
#' \cr
#' The model specification is given below. Note that the model formulae have been adjusted to reflect the fact that JAGS
#' parameterizes the normal and multivariate normal distributions by their precision, rater than (co)variance. 
#' For generalized linear models plug-in pseudovariances are used. 
#' \cr
#' \cr
#' \if{html}{\figure{apc.png}{}}
#' \if{latex}{\figure{apc.png}{}}
#' \cr
#' \cr
#' Plugin Pseudo-Variances: \cr
#' \if{html}{\figure{pseudovar.png}{}}
#' \if{latex}{\figure{pseudovar.png}{}}
#'
#' @references 
#' 
#' Zellner, A. & Siow S. (1980). Posterior odds ratio for selected regression hypotheses. In Bayesian statistics. Proc. 1st int. meeting (eds J. M. Bernardo, M. H. DeGroot, D. V. Lindley & A. F. M. Smith), 585–603. University Press, Valencia. \cr 
#' \cr
#' Liang, Paulo, Molina, Clyde, & Berger (2008) Mixtures of g Priors for Bayesian Variable Selection, Journal of the American Statistical Association, 103:481, 410-423, DOI: 10.1198/016214507000001337 \cr
#' \cr
#' Krishna, A., Bondell, H. D., & Ghosh, S. K. (2009). Bayesian variable selection using an adaptive powered correlation prior. Journal of statistical planning and inference, 139(8), 2665–2674. doi:10.1016/j.jspi.2008.12.004
#' \cr
#' 
#' 
#' @param formula the model formula
#' @param data a data frame
#' @param lower lower limit on value of lambda. Is NULL by default and limits are set based on the minimum value that produces a positive definite covariance matrix.
#' @param uppper upper limit on value of lambda. Is NULL by default and limits are set based on the maximum value that produces a positive definite covariance matrix.
#' @param family one of "gaussian", "binomial", or "poisson".
#' @param log_lik Should the log likelihood be monitored? The default is FALSE.
#' @param iter How many post-warmup samples? Defaults to 15000.
#' @param warmup How many warmup samples? Defaults to 5000.
#' @param adapt How many adaptation steps? Defaults to 5000.
#' @param chains How many chains? Defaults to 4. 
#' @param thin Thinning interval. Defaults to 1.
#' @param method Defaults to "parallel". For an alternative parallel option, choose "rjparallel". Otherwise, "rjags" (single core run).
#' @param cl Use parallel::makeCluster(# clusters) to specify clusters for the parallel methods. Defaults to two cores.
#' @param ... Other arguments to run.jags.
#'
#' @return
#' A run.jags object.
#' @export
#'
#' @examples
#' apcGlm()
#' 
apcGlm = function(formula, data, family = "gaussian", lower = NULL, upper = NULL, log_lik = FALSE, iter = 15000, warmup = 5000, adapt = 5000, chains=4, thin=1, method = "parallel", cl = makeCluster(2), ...)
{
  data <- as.data.frame(data)
  y <- as.numeric(model.frame(formula, data)[, 1])
  X <- as.matrix(model.matrix(formula, data)[,-1])
  # Eigendecomposition
  cormat = cov2cor(fBasics::makePositiveDefinite(cor(X)))
  L = eigen(cormat)$vectors
  D = eigen(cormat)$values
  Trace = function(mat){sum(diag(mat))}
  P = ncol(X)
  
  
  apcLambda = function(formula, data){
    pdcheck = function(formula, data, lambda){
      X = model.matrix(formula, data)[,-1]
      cormat = cov2cor(fBasics::makePositiveDefinite(cor(X)))
      L = eigen(cormat)$vectors
      D = eigen(cormat)$values
      Trace = function(mat){sum(diag(mat))}
      Dpower = matrix(0, length(D), length(D))
      for (i in 1:length(D)){
        Dpower[i,i] <- pow(D[i], lambda)
      }
      fBasics::isPositiveDefinite(L %*% Dpower %*% t(L))
    }
    
    
    pd = as.numeric(sapply(seq(-20, 20, by = 1), function(l) pdcheck(formula, data, l)))
    l = seq(-20, 20, by = 1)[which(pd == 1)]
    c(lower.limit = min(l), upper.limit = max(l))
  }
  
  if (is.null(lower) || is.null(upper)){
    limits = apcLambda(formula, data)
    lower = limits[1]
    upper = limits[2]
  }
  
  if (family == "gaussian"){
    
    jags_apc = "model{
    
              tau ~ dgamma(.01, .01)
              g_inv ~ dgamma(.5, N * .5)
              g <- 1 / g_inv
              sigma <- sqrt(1/tau)
              lambda ~ dunif(lower, upper)
              
              for (i in 1:(P-1)) {
               for (j in (i+1):P) {
                 Dpower[i,j] <- 0
                 Dpower[j,i] <- Dpower[i,j]
                }
              }
              
              for (i in 1:P){
                Dpower[i,i] <- pow(D[i], lambda)
              }
              
              
              prior_cov_pre_raw <- L %*% Dpower %*% t(L)
              
              
              for (i in 1:P){
                for (j in 1:P){
                  prior_cov_raw[i,j] <- prior_cov_pre_raw[i,j] / N
                }
              }
              
              for (i in 1:P){
                d[i] <- prior_cov_raw[i, i]
              }
              
              trace <- sum(d[1:P])
              K <- t / trace
              
              for (i in 1:P){
                for (j in 1:P){
                    prior_cov[i, j] <- g * pow(sigma, 2)* (prior_cov_raw[i, j] * K)
                }
              }

              omega <- inverse(prior_cov) 
              beta[1:P] ~ dmnorm(rep(0,P), omega[1:P,1:P])
              Intercept ~ dnorm(0, 1e-10)
              
              for (i in 1:N){
                 y[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 log_lik[i] <- logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 ySim[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
              }
              Deviance <- -2 * sum(log_lik[1:N])
          }"
    
    P = ncol(X)
    write_lines(jags_apc, "jags_apc.txt")
    jagsdata = list(X = X, y = y,  N = length(y), P = ncol(X), t = Trace(XtXinv(X)), D=D, L=L, lower = lower, upper = upper)
    
    monitor = c("Intercept", "beta", "sigma", "g", "lambda", "Deviance", "ySim" , "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = lmSolve(formula, data)[1], 
                                              "beta" = lmSolve(formula, data)[-1], 
                                              "tau" = 1, 
                                              "g_inv" = 1/length(y), 
                                              "ySim" = y, 
                                              "lambda" = runif(1, lower, upper),
                                              .RNG.name= "lecuyer::RngStream", 
                                              .RNG.seed = sample(1:10000, 1)))
  }
  
  if (family == "binomial" || family == "logistic"){
    
    jags_apc = "model{

              g_inv ~ dgamma(.5, N * .5)
              g <- 1 / g_inv
              lambda ~ dunif(lower, upper)
              
              for (i in 1:(P-1)) {
               for (j in (i+1):P) {
                 Dpower[i,j] <- 0
                 Dpower[j,i] <- Dpower[i,j]
                }
              }
              
              for (i in 1:P){
                Dpower[i,i] <- pow(D[i], lambda)
              }
              
              
              prior_cov_pre_raw <- L %*% Dpower %*% t(L)
              
              
              for (i in 1:P){
                for (j in 1:P){
                  prior_cov_raw[i,j] <- prior_cov_pre_raw[i,j] / N
                }
              }
              
              for (i in 1:P){
                d[i] <- prior_cov_raw[i, i]
              }
              
              trace <- sum(d[1:P])
              K <- t / trace
              
              for (i in 1:P){
                for (j in 1:P){
                    prior_cov[i, j] <- g * sigma2 * (prior_cov_raw[i, j] * K)
                }
              }

              omega <- inverse(prior_cov) 
              beta[1:P] ~ dmnorm(rep(0,P), omega[1:P,1:P])
              Intercept ~ dnorm(0, 1e-10)
              for (i in 1:N){
                 logit(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                 y[i] ~ dbern(psi[i])
                 log_lik[i] <- logdensity.bern(y[i], psi[i])
                 ySim[i] ~ dbern(psi[i])
              }
             Deviance <- -2 * sum(log_lik[1:N])
          }"
    
    P = ncol(X)
    write_lines(jags_apc, "jags_apc.txt")
    jagsdata = list(X = X, y = y, N = length(y), P = ncol(X), t = Trace(XtXinv(X)), D=D, L=L, lower = lower, upper = upper, sigma2 = pow(mean(y), -1) * pow(1 - mean(y), -1), prior_cov = XtXinv(X))
    monitor = c("Intercept", "beta", "g", "Deviance", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = as.vector(coef(glm(formula, data, family = "binomial")))[1], 
                                              "beta" = as.vector(coef(glm(formula, data, family = "binomial")))[-1], 
                                              "g_inv" = 1/length(y), 
                                              "ySim" = y, 
                                              "lambda" = runif(1, lower, upper),
                                              .RNG.name= "lecuyer::RngStream", 
                                              .RNG.seed= sample(1:10000, 1)))
  }
  
  if (family == "poisson"){
    
    jags_apc = "model{
    
              g_inv ~ dgamma(.5, N * .5)
              g <- 1 / g_inv
              lambda ~ dunif(lower, upper)
              
              for (i in 1:(P-1)) {
               for (j in (i+1):P) {
                 Dpower[i,j] <- 0
                 Dpower[j,i] <- Dpower[i,j]
                }
              }
              
              for (i in 1:P){
                Dpower[i,i] <- pow(D[i], lambda)
              }
              
              
              prior_cov_pre_raw <- L %*% Dpower %*% t(L)
              
              
              for (i in 1:P){
                for (j in 1:P){
                  prior_cov_raw[i,j] <- prior_cov_pre_raw[i,j] / N
                }
              }
              
              for (i in 1:P){
                d[i] <- prior_cov_raw[i, i]
              }
              
              trace <- sum(d[1:P])
              K <- t / trace
              
              for (i in 1:P){
                for (j in 1:P){
                    prior_cov[i, j] <- g * sigma2 * (prior_cov_raw[i, j] * K)
                }
              }
              Intercept ~ dnorm(1e-10)
              omega <- inverse(prior_cov) 
              beta[1:P] ~ dmnorm(rep(0,P), omega[1:P,1:P])
              
              for (i in 1:N){
                 log(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                 y[i] ~ dpois(psi[i])
                 log_lik[i] <- logdensity.pois(y[i], psi[i])
                 ySim[i] ~ dpois(psi[i])
              }
              
              Deviance <- -2 * sum(log_lik[1:N])
          }"
    
    write_lines(jags_apc, "jags_apc.txt")
    P = ncol(X)
    jagsdata = list(X = X, y = y, N = length(y),  P = ncol(X), t = Trace(XtXinv(X)), D=D, L=L, lower = lower, upper = upper, sigma2 = pow(mean(y) , -1))
    monitor = c("Intercept", "beta", "g", "Deviance", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = as.vector(coef(glm(formula, data, family = "poisson")))[1], 
                                              "beta" = as.vector(coef(glm(formula, data, family = "poisson")))[1], 
                                              "g_inv" = 1/length(y), 
                                              "ySim" = y, 
                                              "lambda" = runif(1, lower, upper),
                                              .RNG.name= "lecuyer::RngStream", 
                                              .RNG.seed= sample(1:10000, 1)))
  }
  
  out = run.jags(model = "jags_apc.txt", modules = c("glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, n.chains = chains, cl = cl, summarise = FALSE,...)
  if (is.null(cl) == FALSE){
    parallel::stopCluster(cl = cl)
  }
  file.remove("jags_apc.txt")
  return(out)
}
