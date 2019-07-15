#' Adaptive Powered Correlation Prior Stochastic Search Variable Selection with design covariates
#'
#' @description
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
#' that arise when using fixed-g priors.\cr
#' \cr
#' In addition, this function allows for a set of covariates that are held constant across all models.  
#' For example, you may wish to keep variables such as age and gender constant in order to control for them,
#' so that the selected variables are chosen in light of the effects of age and gender on the outcome variable. \cr
#' \cr
#' The probability that a variable has no effect is 1 - mean(delta_i), where delta_i is an indicator variable that takes on values of 1 for
#' inclusion and 0 for exclusion. Averaging the number of 1s over the MCMC iterations gives the posterior inclusion probability (pip), hence, 
#' 1-pip gives the posterior exclusion probability. The overall rate of inclusion for all variables is controlled by the hyperparameter 
#' "phi". Phi is given a beta(1,1) prior which gives uniform probability to the inclusion rate. \cr
#' \cr
#' The posterior means of the coefficients give the Bayesian Model Averaged estimates, which are the expected values of each 
#' parameter averaged over all sampled models (Hoeting et al., 1999). \cr
#' \cr
#' The model specification is given below. Note that the model formulae have been adjusted to reflect the fact that JAGS
#' parameterizes the normal and multivariate normal distributions by their precision, rater than (co)variance. For generalized
#' linear models plug-in pseudovariances are used.
#' \cr
#' \cr
#' Model Specification:
#' \cr
#' \if{html}{\figure{apcSpikeDC.png}{}}
#' \if{latex}{\figure{apcSpikeDC.png}{}}
#' \cr
#' \cr
#' Plugin Pseudo-Variances: \cr
#' \if{html}{\figure{pseudovar.png}{}}
#' \if{latex}{\figure{pseudovar.png}{}}
#' 
#' 
#' @references 
#' Krishna, A., Bondell, H. D., & Ghosh, S. K. (2009). Bayesian variable selection using an adaptive powered correlation prior. Journal of statistical planning and inference, 139(8), 2665–2674. doi:10.1016/j.jspi.2008.12.004 \cr
#' \cr
#' Kuo, L., & Mallick, B. (1998). Variable Selection for Regression Models. Sankhyā: The Indian Journal of Statistics, Series B, 60(1), 65-81. \cr
#' \cr
#' Hoeting, J. , Madigan, D., Raftery, A. & Volinsky, C. (1999). Bayesian model averaging: a tutorial. Statistical Science 14 382–417. \cr
#'
#' @param formula the model formula
#' @param design.formula formula for the design covariates.
#' @param data a data frame
#' @param lower lower limit on value of lambda. Defaults to -10. If the model is failing due to a non-invertible
#' matrix, try adjusting this number.
#' @param uppper upper limit on value of lambda. Defaults to 10. If the model is failing due to a non-invertible
#' matrix, try adjusting this number.
#' @param family one of "gaussian", "binomial", or "poisson".
#' @param log_lik Should the log likelihood be monitored? The default is FALSE.
#' @param iter How many post-warmup samples? Defaults to 10000.
#' @param warmup How many warmup samples? Defaults to 1000.
#' @param adapt How many adaptation steps? Defaults to 2000.
#' @param chains How many chains? Defaults to 4. 
#' @param thin Thinning interval. Defaults to 1.
#' @param method Defaults to "rjparallel". For an alternative parallel option, choose "parallel". Otherwise, "rjags" (single core run).
#' @param cl Use parallel::makeCluster(# clusters) to specify clusters for the parallel methods. Defaults to two cores.
#' @param ... Other arguments to run.jags.
#'
#' @return
#' A run.jags object.
#' @export
#'
#' @examples
#' apcSpikeDC()
#' 
apcSpikeDC = function(formula, design.formula, data, family = "gaussian", lower = -10, upper = 10, log_lik = FALSE, 
                      iter = 10000, warmup=1000, adapt = 5000, chains=4, thin=1, method = "rjparallel", cl = makeCluster(2), ...)
{
  FX <- as.matrix(model.matrix(design.formula, data)[, -1])
  data <- as.data.frame(data)
  y <- as.numeric(model.frame(formula, data)[, 1])
  X <- as.matrix(model.matrix(formula, data)[,-1])
  # Eigendecomposition
  cormat = cov2cor(fBasics::makePositiveDefinite(cor(X)))
  L = eigen(cormat)$vectors
  D = eigen(cormat)$values
  Trace = function(mat){sum(diag(mat))}
  P = ncol(X)
  
  if (family == "gaussian"){
    
    jags_apc = "model{
              
              phi ~ dbeta(1, 1)
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
              
              Intercept ~ dnorm(0, 1e-10)
              
              omega <- inverse(prior_cov) 
              theta[1:P] ~ dmnorm(rep(0,P), omega[1:P,1:P])
              
              for (i in 1:P){
                delta[i] ~ dbern(phi)
                beta[i] <- delta[i] * theta[i]
              }
              
              # Design Variable Coefficients
                  for (f in 1:FP){
                  design_beta[f] ~ dnorm(0, 1e-200)
              }
              
              for (i in 1:N){
                 y[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]) + sum(design_beta[1:FP] * FX[i,1:FP]) , tau)
                 log_lik[i] <- logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]) + sum(design_beta[1:FP] * FX[i,1:FP]) , tau)
                 ySim[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]) + sum(design_beta[1:FP] * FX[i,1:FP]) , tau)
              }
              Deviance <- -2 * sum(log_lik[1:N])
              BIC <- (log(N) * (sum(delta[1:P]) + FP)) + Deviance
          }"
    
    P = ncol(X)
    FP <- ncol(FX)
    write_lines(jags_apc, "jags_apc.txt")
    jagsdata = list(X = X, y = y,  N = length(y), P = ncol(X), t = Trace(XtXinv(X)), D=D, L=L, lower = lower, upper = upper, FP = FP, FX = FX)
    monitor = c("Intercept", "beta", "design_beta", "sigma", "g", "lambda", "phi", "BIC", "Deviance", "delta", "ySim" ,"log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept"= as.vector(coef(lm(formula, data)))[1], 
                                              "design_beta" = as.vector(coef(lm(design.formula, data)))[-1], 
                                              "phi" = rbeta(1, 15, 15), 
                                              "lambda" = runif(1, lower, upper), 
                                              "delta" = rep(0, P), 
                                              "theta" = lmSolve(formula, data)[-1], 
                                              "tau" = 1, 
                                              "g_inv" = 1/length(y), 
                                              "ySim" = y, 
                                              .RNG.name= "lecuyer::RngStream", 
                                              .RNG.seed = sample(1:10000, 1)))
  }
  
  if (family == "binomial"){
    
    jags_apc = "model{
    
              phi ~ dbeta(1, 1)
              g_inv ~ dgamma(.5, N * .5)
              g <- 1 / g_inv
              
              for (j in 1:P){
                for (k in 1:P){
                  cov[j,k] = g * sigma2 * prior_cov[j,k]
                }
              }
              
              omega <- inverse(cov) 
              
              Intercept ~ dnorm(0, 1e-10)
              
              theta[1:P] ~ dmnorm(rep(0,P), omega[1:P,1:P])
              for (i in 1:P){
                delta[i] ~ dbern(phi)
                beta[i] <- delta[i] * theta[i]
              }
              
              
               # Design Variable Coefficients
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
             BIC <- (log(N) * (sum(delta[1:P]) + FP)) + Deviance
          }"
    
    P = ncol(X)
    FP <- ncol(FX)
    write_lines(jags_apc, "jags_apc.txt")
    jagsdata = list(X = X, y = y, N = length(y), P = ncol(X), t = Trace(XtXinv(X)), D=D, L=L, lower = lower, upper = upper, FP = FP, FX = FX, sigma2 = pow(mean(y), -1) * pow(1 - mean(y), -1))
    monitor = c("Intercept", "beta", "design_beta", "g", "lambda", "phi", "BIC", "Deviance", "delta", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = as.vector(coef(glm(formula, data, family = "binomial")))[1], 
                                              "design_beta" = rep(0, FP),  
                                              "phi" = rbeta(1, 15, 15), 
                                              "delta" = rep(0, P), 
                                              "lambda" = runif(1, lower, upper), 
                                              "theta" = as.vector(coef(glm(formula, data, family = "binomial")))[-1], 
                                              "g_inv" = 1/length(y), 
                                              "ySim" = y, 
                                              .RNG.name= "lecuyer::RngStream", 
                                              .RNG.seed= sample(1:10000, 1)))
  }
  
  if (family == "poisson"){
    
    jags_apc = "model{
    
              phi ~ dbeta(1, 1)

              g_inv ~ dgamma(.5, N * .5)
              g <- 1 / g_inv
              
              for (j in 1:P){
                for (k in 1:P){
                  cov[j,k] = g * sigma2 * prior_cov[j,k]
                }
              }
              
              omega <- inverse(cov) 
              
              Intercept ~ dnorm(0, 1e-10)
              
              theta[1:P] ~ dmnorm(rep(0,P), omega[1:P,1:P])
              for (i in 1:P){
                delta[i] ~ dbern(phi)
                beta[i] <- delta[i] * theta[i]
              }
             
             # Design Variable Coefficients
                  for (f in 1:FP){
                  design_beta[f] ~ dnorm(0, 1e-200)
              }
              
              for (i in 1:N){
                 log(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                 y[i] ~ dpois(psi[i])
                 log_lik[i] <- logdensity.pois(y[i], psi[i]) + sum(design_beta[1:FP] * FX[i,1:FP]) 
                 ySim[i] ~ dpois(psi[i])
              }
              
              Deviance <- -2 * sum(log_lik[1:N])
              BIC <- (log(N) * (sum(delta[1:P]) + FP)) + Deviance
          }"
    
    write_lines(jags_apc, "jags_apc.txt")
    P = ncol(X)
    FP <- ncol(FX)
    jagsdata = list(X = X, y = y, N = length(y),  P = ncol(X), t = Trace(XtXinv(X)), D=D, L=L, lower = lower, upper = upper, FP = FP, FX = FX, sigma2 = pow(mean(y) , -1))
    monitor = c("Intercept", "beta", "design_beta", "g", "lambda", "phi", "BIC", "Deviance", "delta", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = as.vector(coef(glm(formula, data, family = "poisson")))[1], 
                                              "design_beta" = rep(0, FP), 
                                              "g_inv" = 1/length(y), 
                                              "ySim" = y, 
                                              .RNG.name= "lecuyer::RngStream", 
                                              .RNG.seed= sample(1:10000, 1),
                                              "phi" = rbeta(1, 15, 15), 
                                              "lambda" = runif(1, lower, upper), 
                                              "delta" = rep(0, P), 
                                              "theta" = as.vector(coef(glm(formula, data, family = "poisson")))[-1]))
  }
  
  
  out = run.jags(model = "jags_apc.txt", 
                 modules = c("glm on", "dic off"), 
                 monitor = monitor, 
                 data = jagsdata, 
                 inits = inits, 
                 burnin = warmup, 
                 sample = iter, 
                 n.chains = chains, 
                 thin = thin, 
                 adapt = adapt, 
                 method = method, 
                 cl = cl, 
                 summarise = FALSE,
                 ...)
  
  if (is.null(cl) == FALSE){
    parallel::stopCluster(cl = cl)
  }
  file.remove("jags_apc.txt")
  return(out)
}
