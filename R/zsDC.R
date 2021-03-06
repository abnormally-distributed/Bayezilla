#' Zellner-Siow g-prior with design covariates
#'
#' The Zellner-Siow cauchy g-prior utilizes the inverse crossproduct is to determine the proper scale of the coefficient priors 
#' by treating the inverse crossproduct of the model matrix as a covariance matrix for a multivariate normal prior distribution 
#' for the coefficients, which is scaled by the parameter "g". The logic is that variables which carry the most information will
#' consequently have a more dispersed prior, while variables that carry less information will have priors more concentrated about zero. 
#' While the joint prior is multivariate normal, the implied independent marginal priors are Cauchy distributions. The approach 
#' here is to let g be a random variable estimated as part of the model, rather than fixed values of g=N. This avoids several
#' problems associated with fixed-g priors. For more information, see Liang et al. (2008). \cr
#' \cr
#' In addition, this function allows for a set of covariates that are held constant across all models.  
#' For example, you may wish to keep variables such as age and gender constant in order to control for them,
#' so that the selected variables are chosen in light of the effects of age and gender on the outcome variable. \cr
#' \cr
#' \if{html}{\figure{zsDC.png}{}}
#' \if{latex}{\figure{zsDC.png}{}}
#' \cr
#' \cr
#' Plugin Pseudo-Variances: \cr
#' \if{html}{\figure{pseudovar.png}{}}
#' \if{latex}{\figure{pseudovar.png}{}}
#'
#' 
#' @references 
#' Zellner, A. & Siow S. (1980). Posterior odds ratio for selected regression hypotheses. In Bayesian statistics. Proc. 1st int. meeting (eds J. M. Bernardo, M. H. DeGroot, D. V. Lindley & A. F. M. Smith), 585–603. University Press, Valencia. \cr 
#' \cr
#' Zellner, A. (1986) On assessing prior distributions and Bayesian regression analysis with g-prior distributions. In P. K. Goel and A. Zellner, editors, Bayesian Inference and Decision Techniques: Essays in Honor of Bruno de Finetti, 233–243.  \cr
#' \cr
#' Liang, Paulo, Molina, Clyde, & Berger (2008) Mixtures of g Priors for Bayesian Variable Selection, Journal of the American Statistical Association, 103:481, 410-423, DOI: 10.1198/016214507000001337 \cr
#' \cr
#' @param formula the model formula
#' @param design.formula the formula for the design covariates
#' @param data a data frame
#' @param family one of "gaussian", "binomial", or "poisson".
#' @param log_lik Should the log likelihood be monitored? The default is FALSE.
#' @param iter How many post-warmup samples? Defaults to 10000.
#' @param warmup How many warmup samples? Defaults to 1000.
#' @param adapt How many adaptation steps? Defaults to 2000.
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
#' zsDC()
#' 
zsDC = function(formula, design.formula, data, family = "gaussian", log_lik = FALSE, iter=10000, warmup=1000, adapt=5000, chains=4, thin=1, method = "parallel", cl = makeCluster(2), ...)
{
  data <- as.data.frame(data)
  FX <- as.matrix(model.matrix(design.formula, data)[, -1])
  y <- as.numeric(model.frame(formula, data)[, 1])
  X <- as.matrix(model.matrix(formula, data)[,-1])
  prior_cov = XtXinv(X)
  
  if (family == "gaussian"){
    
    jags_zs = "model{
    
              tau ~ dgamma(.01, .01)
              g_inv ~ dgamma(.5, N * .5)
              g <- 1 / g_inv
              sigma <- sqrt(1/tau)
              
              for (j in 1:P){
                for (k in 1:P){
                  cov[j,k] = g * pow(sigma, 2) * prior_cov[j,k]
                }
              }
              
              # Design Variable Coefficients
                  for (f in 1:FP){
                  design_beta[f] ~ dnorm(0, 1e-200)
              }
              
              omega <- inverse(cov) 
              beta[1:P] ~ dmnorm(rep(0,P), omega[1:P,1:P])
              Intercept ~ dnorm(0, 1e-10)
       
              for (i in 1:N){
                 y[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]) + sum(design_beta[1:FP] * FX[i,1:FP]) , tau)
                 log_lik[i] <- logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]) + sum(design_beta[1:FP] * FX[i,1:FP]), tau)
                 ySim[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]) + sum(design_beta[1:FP] * FX[i,1:FP]) , tau)
              }
              
              Deviance <- -2 * sum(log_lik[1:N])
          }"
    

    write_lines(jags_zs, "jags_zs.txt")
    P = ncol(X)
    FP = ncol(FX)
    jagsdata = list(X = X, y = y, N = length(y),  P = ncol(X), prior_cov = prior_cov, FP = FP, FX = FX)
    monitor = c("Intercept", "beta", "design_beta", "sigma", "g", "Deviance", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = 0, "beta" = jitter(rep(0, P), amount = 1), "design_beta" = rep(0, FP), "tau" = 1,  "g_inv" = 1/length(y), "ySim" = sample(y, length(y)), .RNG.name= "lecuyer::RngStream", .RNG.seed= sample(1:10000, 1)))
    out = run.jags(model = "jags_zs.txt", modules = c("glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE,...)
}
  if (family == "binomial" || family == "logistic"){
    
    jags_zs = "model{

              g_inv ~ dgamma(.5, N * .5)
              g <- 1 / g_inv
              
              for (j in 1:P){
                for (k in 1:P){
                  cov[j,k] = g * sigma2 * prior_cov[j,k]
                }
              }
              
              # Design Variable Coefficients
                  for (f in 1:FP){
                  design_beta[f] ~ dnorm(0, 1e-200)
              }
              
              omega <- inverse(cov) 
              beta[1:P] ~ dmnorm(rep(0,P), omega[1:P,1:P])
              Intercept ~ dlogis(0, 1e-10)
              for (i in 1:N){
                 logit(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P]) + sum(design_beta[1:FP] * FX[i,1:FP]) 
                 y[i] ~ dbern(psi[i])
                 log_lik[i] <- logdensity.bern(y[i], psi[i])
                 ySim[i] ~ dbern(psi[i])
              }
             Deviance <- -2 * sum(log_lik[1:N])
          }"

    write_lines(jags_zs, "jags_zs.txt")
    P = ncol(X)
    FP <- ncol(FX)
    jagsdata = list(X = X, y = y, N = length(y),  P = ncol(X), prior_cov = prior_cov, FP = FP, FX = FX, sigma2 = pow(mean(y), -1) * pow(1 - mean(y), -1))
    monitor = c("Intercept", "beta", "design_beta", "g", "Deviance", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = 0, "beta" = jitter(rep(0, P), amount = 1), "design_beta" = rep(0, FP), "g_inv" = 1/length(y), "ySim" = sample(y, length(y)), .RNG.name= "lecuyer::RngStream", .RNG.seed= sample(1:10000, 1)))
    out = run.jags(model = "jags_zs.txt", modules = c("glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE,...)
}
  
  if (family == "poisson"){
    
    jags_zs = "model{

              g_inv ~ dgamma(.5, N * .5)
              g <- 1 / g_inv
              
              for (j in 1:P){
                for (k in 1:P){
                  cov[j,k] = g * sigma2 * prior_cov[j,k]
                }
              }

              # Design Variable Coefficients
                  for (f in 1:FP){
                  design_beta[f] ~ dnorm(0, 1e-200)
              }
              
              omega <- inverse(cov) 
              Intercept ~ dnorm(0, 1e-10)
              beta[1:P] ~ dmnorm(rep(0,P), omega[1:P,1:P])
            
              for (i in 1:N){
                 log(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                 y[i] ~ dpois(psi[i])
                 log_lik[i] <- logdensity.pois(y[i], psi[i]) + sum(design_beta[1:FP] * FX[i,1:FP]) 
                 ySim[i] ~ dpois(psi[i])
              }
              
              Deviance <- -2 * sum(log_lik[1:N])
          }"
    
    write_lines(jags_zs, "jags_zs.txt")
    P = ncol(X)
    FP <- ncol(FX)
    jagsdata = list(X = X, y = y, N = length(y),  P = ncol(X), prior_cov = prior_cov, FP = FP, FX = FX, sigma2 = pow(mean(y) , -1))
    monitor = c("Intercept", "beta", "design_beta", "g", "Deviance", "ySim", "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = 0, "beta" = jitter(rep(0, P), amount = 1), "design_beta" = rep(0, FP), "g_inv" = 1/length(y), "ySim" = sample(y, length(y)), .RNG.name= "lecuyer::RngStream", .RNG.seed= sample(1:10000, 1)))
    out = run.jags(model = "jags_zs.txt", modules = c("glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE,...)
  }
  if (is.null(cl) == FALSE){
    parallel::stopCluster(cl = cl)
  }
  file.remove("jags_zs.txt")
  return(out)
}
