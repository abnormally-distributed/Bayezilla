#' Group Extended Bayesian LASSO
#'
#' @description This is the extended Bayesian LASSO presented by
#' Crispin M. Mutshinda and Mikko J. Sillanpää (2010) which is an improvement on the Baysian LASSO
#' of Park & Casella (2008).  It is further adapted such that coefficients can be assigned groups in the 
#' spirit of the Group LASSO (Yuan & Lin, 2006) and Group Bayesian LASSO (Kyung et al., 2010). 
#' \cr
#' Two variants are provided here. \cr
#' \cr
#' The first version of the model is based on the original specification of
#' the extended Bayesian LASSO in Mutshinda and Sillanpää (2010) and the Group Bayesian LASSO described
#' in and Group Bayesian LASSO (Kyung et al., 2010), labeled
#' "classic" in the options here. This requires you to choose upper limits for the uniform priors on both 
#' the top-level shrinkage hyperparameter as well as the local shrinkage parameters. These can be tuned through
#' model comparison if neccessary. 
#' \cr
#' \cr
#' Model Specification:
#' \cr
#' \if{html}{\figure{groupExtLASSO.png}{}}
#' \if{latex}{\figure{groupExtLASSO.png}{}}
#' \cr
#' \cr
#' The second version is the "gamma" prior. This places a gamma(0.50 , 0.01) prior on the
#' top level shrinkage hyperparameter. The individual shrinkage parameters are still given independent uniform(0, local_u) 
#' priors just as in the classic version. 
#' \cr
#' \cr
#' Model Specification:
#' \cr
#' \if{html}{\figure{groupExtLASSO2.png}{}}
#' \if{latex}{\figure{groupExtLASSO2.png}{}}
#' \cr
#' \cr
#' The author of the extended Bayesian Lasso (Sillanpää, personal communication) confirmed that gamma prior on the top
#' level shrinkage parameter 
#' does work well in some settings, which is why I opted to include the "gamma" variant.
#' \cr
#' \cr
#' One really nice feature of the extended Bayesian Lasso is that inclusion probabilities and Bayes Factors for the 
#' coefficients are easily obtained. The prior inclusion probability is given by 1/local_u, so for example
#' uniform(0, 2) priors on the shrinkage parameters indicate a 50% prior inclusion
#' probability. Common in Bayesian variable selection is to use a 20% probability if
#' dealing with a high dimensional problem, so for this choose local_u = 5. If you have genuine prior
#' information you can and should use this to guide your choice. If you are unsure, use model comparison
#' to select which value of u to choose. Inclusion indicators are given by a step function based on the 
#' marginal individual shrinkage parameter, delta = step(1 - eta). Inlcusion probabilities are then given as the number of
#' 1s that appear in the vector of monte carlo samples out of the total number of iterations. This will appear as
#' the mean for each inclusion indicator in the summary. \cr
#' \cr
#' Bayes Factors for each cofficient can then be manually derived using the formula below (Mutshinda, & Sillanpää, 2010; 2012). 
#' \cr
#' \cr
#' \if{html}{\figure{extLASSO_BF.png}{}}
#' \if{latex}{\figure{extLASSO_BF.png}{}}
#' \cr
#' \cr
#' @references 
#' Kyung, M., Gill, J., Ghosh, M., and Casella, G. (2010). Penalized regression, standard errors, and bayesian lassos. Bayesian Analysis, 5(2):369–411. \cr
#' \cr
#' Li, Z.,and Sillanpää, M. J. (2012) Overview of LASSO-related penalized regression methods for quantitative trait mapping and genomic selection. Theoretical and Applied Genetics 125: 419-435. \cr
#' \cr
#' Mutshinda, C.M., & Sillanpää, M.J. (2010). Extended Bayesian LASSO for multiple quantitative trait loci mapping and unobserved phenotype prediction. Genetics, 186 3, 1067-75 . \cr 
#' \cr
#' Mutshinda, C. M., and M. J. Sillanpää (2012) A decision rule for quantitative trait locus detection under the extended Bayesian LASSO model. Genetics 192: 1483-1491. \cr
#' \cr
#' Yuan, Ming; Lin, Yi (2006). Model Selection and Estimation in Regression with Grouped Variables. Journal of the Royal Statistical Society. Series B (statistical Methodology). Wiley. 68 (1): 49–67. doi:10.1111/j.1467-9868.2005.00532.x \cr
#' 
#' 
#' @param X the model matrix. Construct this manually with model.matrix()[,-1]
#' @param y the outcome variable
#' @param idx the group labels. Should be of length = to ncol(model.matrix()[,-1]) with the group assignments
#' for each covariate. Please ensure that you start numbering with 1, and not 0.
#' @param family one of "gaussian", "binomial", or "poisson".
#' @param eta_prior one of "classic (default)", or "gamma".
#' @param local_u This must be assigned a value. Default is 2. 
#' @param top_u If using eta_prior = "classic" this must be assigned a value. Default is 50. 
#' @param log_lik Should the log likelihood be monitored? The default is FALSE.
#' @param iter How many post-warmup samples? Defaults to 10000.
#' @param warmup How many warmup samples? Defaults to 10000.
#' @param adapt How many adaptation steps? Defaults to 15000.
#' @param chains How many chains? Defaults to 4.
#' @param thin Thinning interval. Defaults to 3.
#' @param method Defaults to "parallel". For an alternative parallel option, choose "rjparallel" or. Otherwise, "rjags" (single core run).
#' @param cl Use parallel::makeCluster(# clusters) to specify clusters for the parallel methods. Defaults to two cores.
#' @param ... Other arguments to run.jags.
#'
#' @return A run.jags object
#' @export
#'
#' @examples
#' groupExtLASSO()
#'
groupExtLASSO  = function(X, y, idx, data, family = "gaussian", eta_prior = "classic", local_u = 2, top_u = 50, log_lik = FALSE, iter=10000, warmup = 10000, adapt=15000, chains=4, thin = 3, method = "parallel", cl = makeCluster(2), ...){
  
  nG = max(idx)
  
  if (family == "gaussian" || family == "normal") {
    
    if (eta_prior == "gamma"){
      
      if(is.na(local_u)) stop("Please select an upper limit for the uniform prior on eta.")
      
      jags_group_extended_LASSO = "model{

              # Precision
              tau ~ dgamma(.01, .01)

              # Shrinkage top-level-hyperparameter
              Omega ~ dgamma(0.5 , 0.20)

              Intercept ~ dnorm(0, 1e-10)

              for (g in 1:nG){
              
                # Group Level shrinkage hyperparameter
                eta[g] ~ dunif(0, local_u)
                lambda[g] <- Omega * eta[g]

                # Beta Precision
                w[g] <- pow(lambda[g],2)/2
                beta_var[g] ~ dgamma( (k[g] + 1) * 0.50 , w[g] * 0.50)
                beta_prec[g] <- 1 / beta_var[g]
                
              }
              
              for (p in 1:P){
              
                # Coefficient
                beta[p] ~ dnorm(0, beta_prec[idx[p]])
                
                # Indicator Function
                delta[p] <- step(1 - eta[idx[p]])
                
              }

              for (i in 1:N){
                 y[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 log_lik[i] <- logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 ySim[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
              }
              
              sigma <- sqrt(1/tau)
              Deviance <- -2 * sum(log_lik[1:N])
          }"
      
      P = ncol(X)
      monitor = c("Intercept", "beta", "sigma",  "Omega", "Deviance",   "lambda", "delta", "ySim" ,"log_lik")
      if (log_lik == FALSE){
        monitor = monitor[-(length(monitor))]
      }
      jagsdata = list(X = X, y = y, N = length(y), P = ncol(X), local_u = local_u, idx = idx, nG = max(idx), k = as.vector(table(idx)))
      inits = lapply(1:chains, function(z) list("Intercept" = lmSolve(y ~ ., data = data.frame(y = y, X))[1], 
                                                "ySim" = y, 
                                                "Omega" = 2, 
                                                "beta" = lmSolve(y ~ ., data = data.frame(y = y, X))[-1], 
                                                "eta" = rep(1, nG), 
                                                "beta_var" = abs(jitter(rep(.5, nG), amount = .25)), 
                                                "tau" = 1, 
                                                .RNG.name="lecuyer::RngStream", 
                                                .RNG.seed= sample(1:10000, 1)))
      
      write_lines(jags_group_extended_LASSO, "jags_group_extended_LASSO.txt")
    }
    
    
    if (eta_prior == "classic"){
      
      if(is.na(local_u)) stop("Please select an upper limit for the uniform prior on eta.")
      if (is.na(top_u)) stop("Please select an upper limit for the uniform prior on Omega.")
      jags_group_extended_LASSO = "model{

              # Precision
              tau ~ dgamma(.01, .01)

              # Shrinkage top-level-hyperparameter
              Omega ~ dunif(0, top_u)

              Intercept ~ dnorm(0, 1e-10)

              for (g in 1:nG){
              
                # Group Level shrinkage hyperparameter
                eta[g] ~ dunif(0, local_u)
                lambda[g] <- Omega * eta[g]

                # Beta Precision
                w[g] <- pow(lambda[g],2)/2
                beta_var[g] ~ dgamma( (k[g] + 1) * 0.50 , w[g] * 0.50)
                beta_prec[g] <- 1 / beta_var[g]
                
              }
              
              for (p in 1:P){
              
                # Coefficient
                beta[p] ~ dnorm(0, beta_prec[idx[p]])
                
                # Indicator Function
                delta[p] <- step(1 - eta[idx[p]])
                
              }

              for (i in 1:N){
                 y[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 log_lik[i] <- logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 ySim[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
              }
              sigma <- sqrt(1/tau)
              Deviance <- -2 * sum(log_lik[1:N])
          }"
      
      P = ncol(X)
      monitor = c("Intercept", "beta", "sigma",  "Omega",  "Deviance",  "lambda", "delta", "ySim" ,"log_lik")
      if (log_lik == FALSE){
        monitor = monitor[-(length(monitor))]
      }
      jagsdata = list(X = X, y = y, N = length(y), P = ncol(X), local_u = local_u, top_u = top_u, idx = idx, nG = max(idx), k = as.vector(table(idx)))
      
      inits = lapply(1:chains, function(z) list("Intercept" = lmSolve(y ~ ., data = data.frame(y = y, X))[1], 
                                                "ySim" = y, 
                                                "Omega" = 2, 
                                                "beta" = lmSolve(y ~ ., data = data.frame(y = y, X))[-1], 
                                                "eta" = rep(1, nG), 
                                                "beta_var" = abs(jitter(rep(.5, nG), amount = .25)), 
                                                "tau" = 1, 
                                                .RNG.name="lecuyer::RngStream", 
                                                .RNG.seed= sample(1:10000, 1)))      
      
      write_lines(jags_group_extended_LASSO, "jags_group_extended_LASSO.txt")
    }
  }
  
  if (family == "binomial" || family == "logistic") {
    
    
    if (eta_prior == "gamma"){
      
      if(is.na(local_u)) stop("Please select an upper limit for the uniform prior on eta.")
      
      jags_group_extended_LASSO = "model{

              # Shrinkage top-level-hyperparameter
              Omega ~ dgamma(0.50 , 0.20)

              Intercept ~ dnorm(0, 1e-10)

              for (g in 1:nG){
              
                # Group Level shrinkage hyperparameter
                eta[g] ~ dunif(0, local_u)
                lambda[g] <- Omega * eta[g]

                # Beta Precision
                w[g] <- pow(lambda[g],2)/2
                beta_var[g] ~ dgamma( (k[g] + 1) * 0.50 , w[g] * 0.50)
                beta_prec[g] <- 1 / beta_var[g]
                
              }
              
              for (p in 1:P){
              
                # Coefficient
                beta[p] ~ dnorm(0, beta_prec[idx[p]])
                
                # Indicator Function
                delta[p] <- step(1 - eta[idx[p]])
                
              }

              for (i in 1:N){
                 logit(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                  y[i] ~ dbern(psi[i])
                  log_lik[i] <- logdensity.bern(y[i], psi[i])
                  ySim[i] ~ dbern(psi[i])
              }
              Deviance <- -2 * sum(log_lik[1:N])
            }"
      
      P = ncol(X)
      monitor = c("Intercept", "beta", "Omega", "Deviance",  "lambda", "delta", "ySim" ,"log_lik")
      if (log_lik == FALSE){
        monitor = monitor[-(length(monitor))]
      }
      jagsdata = list(X = X, y = y, N = length(y), P = ncol(X), local_u = local_u, idx = idx, nG = max(idx), k = as.vector(table(idx)))
      inits = lapply(1:chains, function(z) list("Intercept" = 0, 
                                                "ySim" = y, 
                                                "Omega" = 2, 
                                                "beta" = rep(0, P), 
                                                "eta" = rep(1, nG), 
                                                "beta_var" = abs(jitter(rep(.5, nG), amount = .25)), 
                                                .RNG.name="lecuyer::RngStream", 
                                                .RNG.seed= sample(1:10000, 1)))
      write_lines(jags_group_extended_LASSO, "jags_group_extended_LASSO.txt")
    }
    
    else if (eta_prior == "classic"){
      
      if(is.na(local_u)) stop("Please select an upper limit for the uniform prior on eta.")
      if (is.na(top_u)) stop("Please select an upper limit for the uniform prior on Omega.")
      jags_group_extended_LASSO = "model{

              # Shrinkage top-level-hyperparameter
              Omega ~ dunif(0, top_u)

              Intercept ~ dnorm(0, 1e-10)

              for (g in 1:nG){
              
                # Group Level shrinkage hyperparameter
                eta[g] ~ dunif(0, local_u)
                lambda[g] <- Omega * eta[g]

                # Beta Precision
                w[g] <- pow(lambda[g],2)/2
                beta_var[g] ~ dgamma( (k[g] + 1) * 0.50 , w[g] * 0.50)
                beta_prec[g] <- 1 / beta_var[g]
                
              }
              
              for (p in 1:P){
              
                # Coefficient
                beta[p] ~ dnorm(0, beta_prec[idx[p]])
                
                # Indicator Function
                delta[p] <- step(1 - eta[idx[p]])
                
              }

              for (i in 1:N){
                 logit(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                  y[i] ~ dbern(psi[i])
                  log_lik[i] <- logdensity.bern(y[i], psi[i])
                  ySim[i] ~ dbern(psi[i])
              }
              Deviance <- -2 * sum(log_lik[1:N])
            }"
      
      P = ncol(X)
      monitor = c("Intercept", "beta", "Omega", "Deviance",   "lambda", "delta", "ySim" ,"log_lik")
      if (log_lik == FALSE){
        monitor = monitor[-(length(monitor))]
      }
      jagsdata = list(X = X, y = y, N = length(y), P = ncol(X), local_u = local_u, top_u = top_u, idx = idx, nG = max(idx), k = as.vector(table(idx)))
      inits = lapply(1:chains, function(z) list("Intercept" = 0, 
                                                "ySim" = y, 
                                                "Omega" = 2, 
                                                "beta" = rep(0, P), 
                                                "eta" = rep(1, nG), 
                                                "beta_var" = abs(jitter(rep(.5, nG), amount = .25)), 
                                                .RNG.name="lecuyer::RngStream", 
                                                .RNG.seed= sample(1:10000, 1)))      
      write_lines(jags_group_extended_LASSO, "jags_group_extended_LASSO.txt")
    }
    
  }
  
  else if (family == "poisson") {
    
    if (eta_prior == "gamma"){
      
      if (is.na(local_u)) stop("Please select an upper limit for the uniform prior on eta.")
      
      jags_group_extended_LASSO = "model{

              # Shrinkage top-level-hyperparameter
              Omega ~ dgamma(0.50 , 0.20)

              Intercept ~ dnorm(0, 1e-10)

              for (g in 1:nG){
              
                # Group Level shrinkage hyperparameter
                eta[g] ~ dunif(0, local_u)
                lambda[g] <- Omega * eta[g]

                # Beta Precision
                w[g] <- pow(lambda[g],2)/2
                beta_var[g] ~ dgamma( (k[g] + 1) * 0.50 , w[g] * 0.50)
                beta_prec[g] <- 1 / beta_var[g]
                
              }
              
              for (p in 1:P){
              
                # Coefficient
                beta[p] ~ dnorm(0, beta_prec[idx[p]])
                
                # Indicator Function
                delta[p] <- step(1 - eta[idx[p]])
                
              }

              for (i in 1:N){
                 log(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                 y[i] ~ dpois(psi[i])
                 log_lik[i] <- logdensity.pois(y[i], psi[i])
                 ySim[i] ~ dpois(psi[i])
              }
              Deviance <- -2 * sum(log_lik[1:N])
          }"
      
      P = ncol(X)
      monitor = c("Intercept", "beta", "Omega",  "Deviance",  "lambda", "delta", "ySim" ,"log_lik")
      if (log_lik == FALSE){
        monitor = monitor[-(length(monitor))]
      }
      jagsdata = list(X = X, y = y, N = length(y), P = ncol(X), local_u = local_u, idx = idx, nG = max(idx), k = as.vector(table(idx)))
      inits = lapply(1:chains, function(z) list("Intercept" = 0, 
                                                "ySim" = y, 
                                                "Omega" = 2, 
                                                "beta" = rep(0, P), 
                                                "eta" = rep(1, nG), 
                                                "beta_var" = abs(jitter(rep(.5, nG), amount = .25)), 
                                                .RNG.name="lecuyer::RngStream", 
                                                .RNG.seed= sample(1:10000, 1)))      
      write_lines(jags_group_extended_LASSO, "jags_group_extended_LASSO.txt")
    }
    
    if (eta_prior == "classic"){
      
      if(is.na(local_u)) stop("Please select an upper limit for the uniform prior on eta.")
      if (is.na(top_u)) stop("Please select an upper limit for the uniform prior on Omega.")
      jags_group_extended_LASSO = "model{

              # Shrinkage top-level-hyperparameter
              Omega ~ dunif(0, top_u)

              Intercept ~ dnorm(0, 1e-10)
              
              for (g in 1:nG){
              
                # Group Level shrinkage hyperparameter
                eta[g] ~ dunif(0, local_u)
                lambda[g] <- Omega * eta[g]

                # Beta Precision
                w[g] <- pow(lambda[g],2)/2
                beta_var[g] ~ dgamma( (k[g] + 1) * 0.50 , w[g] * 0.50)
                beta_prec[g] <- 1 / beta_var[g]
                
              }
              
              for (p in 1:P){
              
                # Coefficient
                beta[p] ~ dnorm(0, beta_prec[idx[p]])
                
                # Indicator Function
                delta[p] <- step(1 - eta[idx[p]])
                
              }
              
              for (i in 1:N){
                 log(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                 y[i] ~ dpois(psi[i])
                 log_lik[i] <- logdensity.pois(y[i], psi[i])
                 ySim[i] ~ dpois(psi[i])
              }
              Deviance <- -2 * sum(log_lik[1:N])
          }"
      
      P = ncol(X)
      monitor = c("Intercept", "beta", "Omega", "Deviance",  "lambda", "delta", "ySim" ,"log_lik")
      if (log_lik == FALSE){
        monitor = monitor[-(length(monitor))]
      }
      jagsdata = list(X = X, y = y, N = length(y), P = ncol(X), local_u = local_u, top_u = top_u, idx = idx, nG = max(idx), k = as.vector(table(idx)))
      inits = lapply(1:chains, function(z) list("Intercept" = 0, 
                                                "ySim" = y, 
                                                "Omega" = 2, 
                                                "beta" = rep(0, P), 
                                                "eta" = rep(1, nG), 
                                                "beta_var" = abs(jitter(rep(.5, nG), amount = .25)), 
                                                .RNG.name="lecuyer::RngStream", 
                                                .RNG.seed= sample(1:10000, 1)))      
      write_lines(jags_group_extended_LASSO, "jags_group_extended_LASSO.txt")
    }
  }
  
  out = run.jags(model = "jags_group_extended_LASSO.txt", modules = c("bugs on", "glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, n.chains = chains, method = method, cl = cl, summarise = FALSE, ...)
  file.remove("jags_group_extended_LASSO.txt")
  if (is.null(cl) == FALSE){
    parallel::stopCluster(cl = cl)
  }
  return(out)
}
