#' Bernoulli-Normal Mixture Group Selection for GLMs
#'
#' @description 
#' 
#' IMPORTANT NOTICE: This model works best on smaller to medium sized data sets with a small number of variables (less than 20).
#' If you experience difficulty with running times or obtaining a good effective sample size consider using the extended LASSO.
#' 
#' IMPORTANT NOTICE: Center and scale your predictors before using this function.
#' If you do not scale and center your numeric predictors, this will likely not give reasonable results. 
#'
#' This is the most basic type of Bayesian variable selection. This models the
#' regression coefficients as coming from either a null distribution represented
#' by a probability mass of 100% at zero (the "spike"), or as coming from a Normal( mu=0 , sigma = sqrt(2) ) distribution. This variant of the
#' Bernoulli-Normal mixture prior models the selection of parameters as groups, akin to the group LASSO.
#' 
#'
#' Some Suggestions for Priors on phi (the overall inclusion probability):
#'
#' beta(0.5, 0.5) Jeffrey's Prior (Truly uninformative)
#'
#' beta(1.0, 1.0) Laplace's Uniform Prior
#'
#' beta(1.0, 4.0) Weakly Informative, mean probability = 0.20 [This is the default prior used here]
#'
#' beta(2.0, 2.0) Weakly Informative, mean probability = 0.50 
#'
#' beta(4.0, 1.0) Weakly Informative, mean probability = 0.80
#' 
#' beta(2.0, 8.0) Moderately Informative, Regularizing, mean probability = 0.20 
#' 
#' beta(4.0, 4.0) Moderately Informative, Regularizing, mean probability = 0.50
#'
#' beta(8.0, 2.0) Moderately Informative, Regularizing, mean probability = 0.80
#'
#' @param X the model matrix. Construct this manually with model.matrix()[,-1]
#' @param y the outcome variable
#' @param idx the group labels. Should be of length = to ncol(model.matrix()[,-1]) with the group assignments
#' for each covariate. Please ensure that you start numbering with 1, and not 0.
#' @param family one of "gaussian", "binomial", or "poisson".
#' @param phi_prior Default is c(1, 4).
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
#' @return A run.jags object
#' @export
#'
#' @examples
#' groupSpike()
#'
groupSpike  = function(formula, data, family = "gaussian", phi_prior = c(1, 4), log_lik = FALSE, iter=10000, warmup=1000, adapt=2000, chains=4, thin=1, method = "parallel", cl = makeCluster(2), ...){

  if (family == "gaussian"){

    jags_group_glm_spike = "model{
              tau ~ dgamma(.001, .001)
              phi ~ dbeta(a, b)

              # Indicator Variables For Groups
              for (r in 1:nG){
                 delta[r] ~ dbern(phi)
              }

              # Coefficients
              for (p in 1:P){
                theta[p] ~ dnorm(0, .5)
                beta[p] <- delta[idx[p]] * theta[p]
              }

              Intercept ~ dnorm(0, .01)

              for (i in 1:N){
                 y[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 log_lik[i] <- logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 ySim[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
              }

              sigma <- sqrt(1/tau)
              Deviance <- -2 * sum(log_lik[1:N])
          }"

    P = ncol(X)
    nG <- length(unique(idx))
    write_lines(jags_group_glm_spike, "jags_group_glm_spike.txt")
    jagsdata = list(X = X, y = y, N = length(y), P = ncol(X), a = phi_prior[1], b = phi_prior[2], idx = idx, nG = nG)
    monitor = c("Intercept", "beta", "sigma", "phi", "Deviance", "delta",  "theta" , "ySim" , "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = 0, .RNG.name= "lecuyer::RngStream", .RNG.seed = sample(1:10000, 1),   "ySim" = y, "delta"=rep(1, nG), "phi" = .20 , "theta" = jitter(rep(0, P), amount = .25), "tau" = 1))
    out = run.jags(model = "jags_group_glm_spike.txt", modules = c("glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE,...)
    return(out)
  }

  if (family == "binomial" || family == "logistic"){

    jags_group_glm_spike = "model{
    
              phi ~ dbeta(a, b)
              
              # Indicator Variables For Groups
              for (r in 1:nG){
                 delta[r] ~ dbern(phi)
              }

              # Coefficients
              for (p in 1:P){
                theta[p] ~ dnorm(0, .5)
                beta[p] <- delta[idx[p]] * theta[p]
              }

              Intercept ~ dnorm(0, .01)
              for (i in 1:N){
                 logit(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                 y[i] ~ dbern(psi[i])
                 log_lik[i] <- logdensity.bern(y[i], psi[i])
                 ySim[i] ~ dbern(psi[i])
              }
             Deviance <- -2 * sum(log_lik[1:N])
          }"

    P = ncol(X)
    write_lines(jags_group_glm_spike, "jags_group_glm_spike.txt")
    y = as.numeric(as.factor(y)) - 1
    nG <- length(unique(idx))
    jagsdata = list(X = X, y = y, N = length(y), P = ncol(X), a = phi_prior[1], b = phi_prior[2], idx = idx, nG = nG)
    monitor = c("Intercept", "beta", "phi", "Deviance","delta",  "theta" ,"ySim" , "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = 0, .RNG.name= "lecuyer::RngStream", .RNG.seed = sample(1:10000, 1),  "ySim" = y, "delta" = rep(1, nG), "phi" = .20 , "theta" = jitter(rep(0, P), amount = .25)))
    out = run.jags(model = "jags_group_glm_spike.txt", modules = c("glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE,...)
    return(out)
  }

  if (family == "poisson"){

    jags_group_glm_spike = "model{
              phi ~ dbeta(a, b)

              # Indicator Variables For Groups
              for (r in 1:nG){
                 delta[r] ~ dbern(phi)
              }

              # Coefficients
              for (p in 1:P){
                theta[p] ~ dnorm(0, .5)
                beta[p] <- delta[idx[p]] * theta[p]
              }

              Intercept ~ dnorm(0, .01)

              for (i in 1:N){
                 log(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                 y[i] ~ dpois(psi[i])
                 log_lik[i] <- logdensity.pois(y[i], psi[i])
                 ySim[i] ~ dpois(psi[i])
              }
              Deviance <- -2 * sum(log_lik[1:N])
          }"

    write_lines(jags_group_glm_spike, "jags_group_glm_spike.txt")
    P = ncol(X)
    nG <- length(unique(idx))
    jagsdata = list(X = X, y = y, N = length(y),  P = ncol(X), a = phi_prior[1], b = phi_prior[2], idx = idx, nG = nG)
    monitor = c("Intercept", "beta", "phi", "Deviance","delta",  "theta" , "ySim" , "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = 0, .RNG.name= "lecuyer::RngStream", .RNG.seed = sample(1:10000, 1), "ySim" = y, "delta"=rep(1, nG), "phi" = .20 , "theta" = jitter(rep(0, P), amount = .25)))
    out = run.jags(model = "jags_group_glm_spike.txt", modules = c("glm on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE,...)
    file.remove("jags_group_glm_spike.txt")
    if (!is.null(cl)) {
      parallel::stopCluster(cl = cl)
    }
    return(out)
  }
}


