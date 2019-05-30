#' Bernoulli-Normal Mixture Group Selection for GLMs
#'
#' @description This is the most basic type of Bayesian variable selection, however, it is
#' very intuitive to understand and often performs well**. This models the
#' regression coefficients as coming from either a null distribution represented
#' by a probability mass of 100% at zero, or as coming from a cauchy distribution parameterized
#' as a normal-gamma(.5,.5) mixture with the mixing parameter given the name "omega". This variant of the
#' Bernoulli-Normal mixture model models the selection of parameters as groups, akin to the group LASSO.
#'
#' Note that you should center and scale your variables before using this function.
#' If you do not scale and center your numeric predictors, this will likely not perform well or
#' give reasonable results. The mixing hyperparameter omega assumes all covariates are on the same scale.
#'
#'
#' ** This model works best on smaller to medium sized data sets. If you experience difficulty with
#' running times or obtaining a good effective sample size consider using the extended LASSO. Another tip that
#' may improve performance is using a beta(1, 1) or beta(1.5, 1.5) prior on phi. If all coefficients are set to
#' zero and there is truly good reason to believe this is not correct (ie, you aren't just hunting for "statistical
#' significance", are you?) a prior such as beta(4, 2) may help. If there is no sparsity, and you have genuine reason
#' to believe there should be try a beta(2,8) prior (if not, ask yourself why are you using variable selection?
#' If it is to deal with collinearity or other foibles for least squares, consider trying the glm_bayes function instead).
#'
#'
#' @param X the model matrix. Construct this manually with model.matrix()[,-1]
#' @param y the outcome variable
#' @param idx the group labels. Should be of length = to ncol(model.matrix()[,-1]) with the group assignments
#' for each covariate. Please ensure that you start numbering with 1, and not 0.
#' @param family one of "gaussian", "binomial", or "poisson".
#' @param phi_prior The beta distribution parameters on the inclusion probabilities. Default is Jeffrey's prior, c(0.5, 0.5).
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
#' groupSpike()
#'
groupSpike  = function(formula, data, family = "gaussian", phi_prior = c(.5, .5), log_lik = FALSE, iter=10000, warmup=1000, adapt=2000, chains=4, thin=3, method = "rjags", cl = NA, ...){

  if (family == "gaussian"){

    jags_group_glm_spike = "model{
              tau ~ dgamma(.001, .001)
              phi ~ dbeta(a, b)
              omega ~ dgamma(.5, .5)

              # Indicator Variables For Groups
              for (r in 1:nG){
                 delta[r] ~ dbern(phi)
              }

              # Coefficients
              for (p in 1:P){
                theta[p] ~ dnorm(0, omega)
                beta[p] <- delta[idx[p]] * theta[p]
              }

              Intercept ~ dnorm(0, .01)

              for (i in 1:N){
                 y[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 log_lik[i] <- logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 ySim[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
              }

              sigma <- sqrt(1/tau)
              deviance <- -2 * sum(log_lik[1:N])
          }"

    P = ncol(X)
    nG <- length(unique(idx))
    write_lines(jags_group_glm_spike, "jags_group_glm_spike.txt")
    jagsdata = list(X = X, y = y, N = length(y), P = ncol(X), a = phi_prior[1], b = phi_prior[2], idx = idx, nG = nG)
    monitor = c("Intercept", "beta", "sigma", "phi","omega", "delta", "deviance", "ySim" , "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = 0, "omega" = .001,  "ySim" = y, "delta"=rep(1, nG), "phi" = .20 , "theta" = jitter(rep(0, P), amount = .25), "tau" = 1))
    out = run.jags(model = "jags_group_glm_spike.txt", modules = "glm", monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, ...)
    return(out)
  }

  if (family == "binomial" || family == "logistic"){

    jags_group_glm_spike = "model{
              omega ~ dgamma(.5, .5)
              phi ~ dbeta(a, b)

              # Indicator Variables For Groups
              for (r in 1:nG){
                 delta[r] ~ dbern(phi)
              }

              # Coefficients
              for (p in 1:P){
                theta[p] ~ dnorm(0, omega)
                beta[p] <- delta[idx[p]] * theta[p]
              }

              Intercept ~ dnorm(0, .01)
              for (i in 1:N){
                 logit(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                 y[i] ~ dbern(psi[i])
                 log_lik[i] <- logdensity.bern(y[i], psi[i])
                 ySim[i] ~ dbern(psi[i])
              }
             deviance <- -2 * sum(log_lik[1:N])
          }"

    P = ncol(X)
    write_lines(jags_group_glm_spike, "jags_group_glm_spike.txt")
    y = as.numeric(as.factor(y)) - 1
    nG <- length(unique(idx))
    jagsdata = list(X = X, y = y, N = length(y), P = ncol(X), a = phi_prior[1], b = phi_prior[2], idx = idx, nG = nG)
    monitor = c("Intercept", "beta", "phi","omega", "delta", "deviance", "ySim" , "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = 0, "omega" = .001,  "ySim" = y, "delta" = rep(1, nG), "phi" = .20 , "theta" = jitter(rep(0, P), amount = .25)))
    out = run.jags(model = "jags_group_glm_spike.txt", modules = "glm", monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, ...)
    return(out)
  }

  if (family == "poisson"){

    jags_group_glm_spike = "model{
              omega ~ dgamma(.5, .5)
              phi ~ dbeta(a, b)

              # Indicator Variables For Groups
              for (r in 1:nG){
                 delta[r] ~ dbern(phi)
              }

              # Coefficients
              for (p in 1:P){
                theta[p] ~ dnorm(0, omega)
                beta[p] <- delta[idx[p]] * theta[p]
              }

              Intercept ~ dnorm(0, .01)

              for (i in 1:N){
                 log(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                 y[i] ~ dpois(psi[i])
                 log_lik[i] <- logdensity.pois(y[i], psi[i])
                 ySim[i] ~ dpois(psi[i])
              }
              deviance <- -2 * sum(log_lik[1:N])
          }"

    write_lines(jags_group_glm_spike, "jags_group_glm_spike.txt")
    P = ncol(X)
    nG <- length(unique(idx))
    jagsdata = list(X = X, y = y, N = length(y),  P = ncol(X), a = phi_prior[1], b = phi_prior[2], idx = idx, nG = nG)
    monitor = c("Intercept", "beta", "phi","omega", "delta", "deviance", "ySim" , "log_lik")
    if (log_lik == FALSE){
      monitor = monitor[-(length(monitor))]
    }
    inits = lapply(1:chains, function(z) list("Intercept" = 0, "omega" = .001,  "ySim" = y, "delta"=rep(1, nG), "phi" = .20 , "theta" = jitter(rep(0, P), amount = .25)))
    out = run.jags(model = "jags_group_glm_spike.txt", modules = "glm", monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, ...)
    return(out)
  }
}


