#' Bernoulli-Normal Adaptive Powered Correlation Prior for Variable Selection
#'
#' @description IMPORTANT: I suggest not using any factor predictor variables, only numeric. In my experience the inclusion of categorical
#' predictors tends to lead to odd results in calculating the prior scale.
#' 
#' This function impements the adaptive powered correlation g-prior model described by Krishna, Bondell, and Ghosh (2009). Typically, in regression the cross-product XtX is inverted in the process of calculating the coefficients. In addition, 
#' the Zellner-Siow cauchy g-prior utilizes the inverse crossproduct is used as an empirical Bayesian method of determining the proper scale of the coefficient
#' priors by treating the inverse crossproduct as a covariance matrix, which is scaled by the parameter "g". 
#' 
#' The adaptive powered correlation prior simply extends this to allow using other powers besides -1. The power here will be referred to as "lambda".
#' Setting lambda to 0 results in a ridge-regression like prior. Setting lambda to a positive value adapts to collinearity by allowing
#' correlated predictors to enter and exit the model together. Negative values of lambda on the other hand favor including only one
#' of a set of correlated predictors. Of course, setting lambda to -1 is just the Zellner-Siow cauchy g-prior. This is designed to deal with collinearity in a more adaptive way than even ridge regression
#' by allowing the analyst to use model comparison techniques to choose an optimal value of lambda, and then using the best model for inference.
#' An analysts beliefs about which type of selection is preferred, or the goals of the particular analysis, can also inform the choice of lambda.
#' 
#' The difference between this and the apcSpike function is that there is a hyperprior on the overall inclusion probability prior.
#' This can be useful to encourage better sampling if the apcSpike function is having troubles, or you may wish
#' to use it if you do not wish to assume a beta prior on phi with fixed shape parameters. The structure of this hyperprior is as follows:
#' 
#' w ~ beta(1, 1) \cr
#' k_raw ~ uniform(0, 2) \cr
#' k <- k_raw + 2 \cr
#' a <- w * (K - 2) + 1 \cr
#' b <- (1-w) * (K - 2) + 1 \cr
#' phi ~ beta(a, b) \cr
#' 
#' for (i in 1:P){delta[i] ~ bernoulli(phi)} \cr
#' 
#' The specification of the k_raw as having a lower bound of zero, then adding 2 to k_raw only to subtract 2 from
#' k in the formula for the shape parameters a and b of the beta distribution might seem unnecessary, but coding it
#' this way prevents compilation errors.
#' 
#' Note, however, that this prior is designed to deal with collinearity but not necessarily P > N scenarios. For that you may wish to take a look
#' at the \code{\link[Bayezilla]{extLASSO}} or \code{\link[Bayezilla]{bayesEnet}} functions. 
#'
#'
#' @references 
#' Krishna, A., Bondell, H. D., & Ghosh, S. K. (2009). Bayesian variable selection using an adaptive powered correlation prior. Journal of statistical planning and inference, 139(8), 2665–2674. doi:10.1016/j.jspi.2008.12.004
#' 
#' Kuo, L., & Mallick, B. (1998). Variable Selection for Regression Models. Sankhyā: The Indian Journal of Statistics, Series B, 60(1), 65-81.
#' 
#' @param formula the model formula
#' @param data a data frame
#' @param lambda the power to use in the adaptive correlation prior. Default is -1, which gives the Zellner-Siow g prior. I suggest 
#' fitting multiple values of lambda and selecting which is best via LOO-IC or WAIC. 
#' @param family one of "gaussian", "binomial", or "poisson".
#' @param log_lik Should the log likelihood be monitored? The default is FALSE.
#' @param iter How many post-warmup samples? Defaults to 8000.
#' @param warmup How many warmup samples? Defaults to 2000.
#' @param adapt How many adaptation steps? Defaults to 1500.
#' @param chains How many chains? Defaults to 4. 
#' @param thin Thinning interval. Defaults to 1.
#' @param method Defaults to "rjparallel". For an alternative parallel option, choose "parallel". Otherwise, "rjags" (single core run).
#' @param cl Use parallel::makeCluster(# clusters) to specify clusters for the parallel methods. Defaults to three cores.
#' @param ... Other arguments to run.jags.
#'
#' @return
#' A run.jags object.
#' @export
#'
#' @examples
#' apcSpikeHS()
#' 
apcSpikeHS = function(formula, data, family = "gaussian", lambda = -1, log_lik = FALSE, 
                    iter = 8000, warmup = 2000, adapt = 1500, chains= 4, thin=1, method = "rjparallel", cl = makeCluster(3), ...)
{
  data <- as.data.frame(data)
  y <- as.numeric(model.frame(formula, data)[, 1])
  X <- model.matrix(formula, data)[, -1]
  ## Ensure that the correlation matrix is positive definite.
  cormat = cor(X)
  diag(cormat) <- diag(cormat) + 1e-2
  cormat = cov2cor(solve(pseudoinverse(cov2cor(cormat))))
  # Eigendecomposition
  L = eigen(cormat)$vectors
  D = eigen(cormat)$values
  Trace = function(mat){sum(diag(mat))}
  P = ncol(X)
  Dpower = rep(0, P)
  t = XtXinv(X, tol=1e-2)
  for(i in 1:P) {
    Dpower[i] <- (D[i]^lambda);
  }
  prior_cov = (L %*% diag(Dpower) %*% t(L)) / length(y)
  K = Trace(t) / Trace(prior_cov)
  prior_cov = K * (prior_cov)
  
  if (family == "gaussian"){
    
    jags_apc = "model{
              
              # Priors on the overall inclusion probability 
              w ~ dbeta(1, 1)
              k_raw ~ dunif(0, 4)
              k <- k_raw + 2
              a <- w * (k-2) + 1
              b <- (1-w) * (k-2) + 1
              phi ~ dbeta(a, b)
              
              # Priors on precision and g
              tau ~ dgamma(.0001, .0001)
              g_inv ~ dgamma(.5, N * .5)
              g <- 1 / g_inv
              sigma <- sqrt(1/tau)
              
              # Construct Covariance Matrix for Prior
              for (j in 1:P){
                for (k in 1:P){
                  cov[j,k] = g * pow(sigma, 2) * prior_cov[j,k]
                }
              }
              
              # Convert Cov. Matrix to Precision Matrix
              omega <- inverse(cov) 
              
              # Multivariate prior on the raw coefficients
              theta[1:P] ~ dmnorm(rep(0,P), omega[1:P,1:P])
              
              # Inclusion indicator variable and censored coefficients 
              for (i in 1:P){
                delta[i] ~ dbern(phi)
                beta[i] <- delta[i] * theta[i]
              }
              
              
              Intercept ~ dnorm(0, .001)
              
              # Likelihood Function
              for (i in 1:N){
                 y[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 log_lik[i] <- logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 ySim[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
              }
              
              Deviance <- -2 * sum(log_lik[1:N])
          }"
    
    P = ncol(X)
    write_lines(jags_apc, "jags_apc.txt")
    jagsdata = list(X = X, y = y,  N = length(y), P = ncol(X), prior_cov = prior_cov)
    monitor = c("Intercept", "beta", "sigma", "g", "Deviance", "w", "k", "a", "b", "phi", "delta", "ySim" ,"log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = 0, "phi" = .5, "k_raw" = 2.5, "w" = .5, "delta" = sample(c(0,1), P, replace = TRUE), "theta" = jitter(rep(0, P), amount = 1), "tau" = 1, "g_inv" = 1/length(y), "ySim" = y, .RNG.name= "lecuyer::RngStream", .RNG.seed = sample(1:10000, 1)))
  }
  
  if (family == "binomial" || family == "logistic"){
    
    jags_apc = "model{
    
              
              # Priors on the overall inclusion probability 
              w ~ dbeta(1, 1)
              k_raw ~ dunif(0, 4)
              k <- k_raw + 2
              a <- w * (k-2) + 1
              b <- (1-w) * (k-2) + 1
              phi ~ dbeta(a, b)
              
              # Priors on g
              g_inv ~ dgamma(.5, N * .5)
              g <- 1 / g_inv
              
              # Construct Covariance Matrix for Prior
              for (j in 1:P){
                for (k in 1:P){
                  cov[j,k] = g * 1.0 * prior_cov[j,k]
                }
              }
              
              # Convert Cov. Matrix to Precision Matrix
              omega <- inverse(cov) 
              
              # Multivariate prior on the raw coefficients
              theta[1:P] ~ dmnorm(rep(0,P), omega[1:P,1:P])
              
              # Inclusion indicator variable and censored coefficients 
              for (i in 1:P){
                delta[i] ~ dbern(phi)
                beta[i] <- delta[i] * theta[i]
              }
              
              Intercept ~ dnorm(0, .001)
              
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
    jagsdata = list(X = X, y = y, N = length(y), P = ncol(X), prior_cov = prior_cov)
    monitor = c("Intercept", "beta", "g", "Deviance", "w", "k", "a", "b", "phi", "delta", "ySim" ,"log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = 0, "phi" = .5, "k_raw" = 2.5, "w" = .5, "delta" = sample(c(0,1), P, replace = TRUE), "theta" = jitter(rep(0, P), amount = 1), "g_inv" = 1/length(y), "ySim" = y, .RNG.name= "lecuyer::RngStream", .RNG.seed = sample(1:10000, 1)))
  }
  
  if (family == "poisson"){
    
    jags_apc = "model{
    
              # Priors on the overall inclusion probability 
              w ~ dbeta(1, 1)
              k_raw ~ dunif(0, 4)
              k <- k_raw + 2
              a <- w * (k-2) + 1
              b <- (1-w) * (k-2) + 1
              phi ~ dbeta(a, b)
              
              # Priors on g
              g_inv ~ dgamma(.5, N * .5)
              g <- 1 / g_inv
              
              # Construct Covariance Matrix for Prior
              for (j in 1:P){
                for (k in 1:P){
                  cov[j,k] = g * 1.0 * prior_cov[j,k]
                }
              }
              
              # Convert Cov. Matrix to Precision Matrix
              omega <- inverse(cov) 
              
              # Multivariate prior on the raw coefficients
              theta[1:P] ~ dmnorm(rep(0,P), omega[1:P,1:P])
              
              # Inclusion indicator variable and censored coefficients 
              for (i in 1:P){
                delta[i] ~ dbern(phi)
                beta[i] <- delta[i] * theta[i]
              }
              
              Intercept ~ dnorm(0, .001)
              
              theta[1:P] ~ dmnorm(rep(0,P), omega[1:P,1:P])
              for (i in 1:P){
                delta[i] ~ dbern(phi)
                beta[i] <- delta[i] * theta[i]
              }
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
    jagsdata = list(X = X, y = y, N = length(y),  P = ncol(X), prior_cov = prior_cov)
    monitor = c("Intercept", "beta", "g", "Deviance", "w", "k", "a", "b", "phi", "delta", "ySim" ,"log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = 0 , "phi" = .5, "k_raw" = 2.5, "w" = .5, "delta" = sample(c(0,1), P, replace = TRUE), "theta" = jitter(rep(0, P), amount = 1), "g_inv" = 1/length(y), "ySim" = y, .RNG.name= "lecuyer::RngStream", .RNG.seed = sample(1:10000, 1)))
  }
  
  out = run.jags(model = "jags_apc.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE,...)
  if (is.null(cl) == FALSE){
    parallel::stopCluster(cl = cl)
  }
  file.remove("jags_apc.txt")
  return(out)
}
