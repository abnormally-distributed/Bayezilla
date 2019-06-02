#' Bernoulli-Normal Mixture Selection for GLMs
#'
#' @description
#'
#' IMPORTANT NOTICE: This model works best on smaller to medium sized data sets with a small number of variables (less than 20).
#' If you experience difficulty with running times or obtaining a good effective sample size consider using the extended LASSO.
#'
#' IMPORTANT NOTICE: Center and scale your predictors before using this function.
#' If you do not scale and center your numeric predictors, this will likely not give reasonable results. 
#' 
#' 
#' This is the most basic type of Bayesian variable selection. This models the
#' regression coefficients as coming from either a null distribution represented
#' by a probability mass of 100% at zero (the "spike"), or as coming from a Normal( mu=0 , sigma = sqrt(2) ) distribution. 
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
#'
#' @param formula the model formula
#' @param data a data frame
#' @param family one of "gaussian", "binomial", or "poisson".
#' @param phi_prior Default is c(1, 4).
#' @param log_lik Should the log likelihood be monitored? The default is FALSE.
#' @param iter How many post-warmup samples? Defaults to 10000.
#' @param warmup How many warmup samples? Defaults to 1000.
#' @param adapt How many adaptation steps? Defaults to 2000.
#' @param chains How many chains? Defaults to 4. Max allowed is 4.
#' @param thin Thinning interval. Defaults to 3.
#' @param method Defaults to "parallel". For an alternative parallel option, choose "rjparallel". Otherwise, "rjags" (single core run).
#' @param cl Use parallel::makeCluster(# clusters) to specify clusters for the parallel methods. Defaults to two cores.
#' @param ... Other arguments to run.jags.
#'
#' @return A run.jags object
#' @export
#'
#' @examples
#' glmSpike()
Spike <- function(formula, data, family = "gaussian", phi_prior = c(1, 4), log_lik = FALSE, iter = 10000, warmup = 1000, adapt = 2000, chains = 4, thin = 3, method = "parallel", cl = makeCluster(2),summarise = FALSE, ...) {
  X <- model.matrix(formula, data)[, -1]
  y <- model.frame(formula, data)[, 1]

  RNGlist <- c("base::Wichmann-Hill", "base::Marsaglia-Multicarry", "base::Super-Duper", "base::Mersenne-Twister")
  if (chains > 4) {
    chains <- 4
  }

  if (family == "gaussian") {
    jags_glm_spike <- "model{
              tau ~ dgamma(.001, .001)
              phi ~ dbeta(a, b)
              omega ~ dgamma(3, 3)
              
              for (p in 1:P){
                theta[p] ~ dnorm(0, omega)
                delta[p] ~ dbern(phi)
                beta[p] <- delta[p] * theta[p]
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

    P <- ncol(X)
    write_lines(jags_glm_spike, "jags_glm_spike.txt")
    jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), a = phi_prior[1], b = phi_prior[2])
    monitor <- c("Intercept", "beta", "sigma", "phi", "deviance",  "delta", "theta" ,"ySim", "log_lik")
    if (log_lik == FALSE) {
      monitor <- monitor[-(length(monitor))]
    }
    inits <- lapply(1:chains, function(z) list("Intercept" = 0, .RNG.name = RNGlist[z], .RNG.seed = sample(1:10000, 1),   "ySim" = y, "delta" = rep(1, P), "phi" = .20, "theta" = jitter(rep(0, P), amount = .25), "tau" = 1))
    out <- run.jags(model = "jags_glm_spike.txt", modules = "glm", monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl,summarise = FALSE, ...)
    return(out)
  }

  if (family == "binomial" || family == "logistic") {
    jags_glm_spike <- "model{
              phi ~ dbeta(a, b)
              omega ~ dgamma(3, 3)
              for (p in 1:P){
                theta[p] ~ dnorm(0, omega)
                delta[p] ~ dbern(phi)
                beta[p] <- delta[p] * theta[p]
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

    P <- ncol(X)
    write_lines(jags_glm_spike, "jags_glm_spike.txt")
    y <- as.numeric(as.factor(y)) - 1
    jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), a = phi_prior[1], b = phi_prior[2])
    monitor <- c("Intercept", "beta", "phi",  "delta", "deviance", "theta" , "ySim", "log_lik")
    if (log_lik == FALSE) {
      monitor <- monitor[-(length(monitor))]
    }
    inits <- lapply(1:chains, function(z) list("Intercept" = 0, .RNG.name = RNGlist[z], .RNG.seed = sample(1:10000, 1),   "ySim" = y, "delta" = rep(1, P), "phi" = .20, "theta" = jitter(rep(0, P), amount = .25)))
    out <- run.jags(model = "jags_glm_spike.txt", modules = "glm", monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE,...)
    return(out)
  }

  if (family == "poisson") {
    jags_glm_spike <- "model{
              phi ~ dbeta(a, b)
              omega ~ dgamma(3, 3)
              for (p in 1:P){
                theta[p] ~ dnorm(0, omega)
                delta[p] ~ dbern(phi)
                beta[p] <- delta[p] * theta[p]
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

    write_lines(jags_glm_spike, "jags_glm_spike.txt")
    P <- ncol(X)
    jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), a = phi_prior[1], b = phi_prior[2])
    monitor <- c("Intercept", "beta", "phi",  "delta", "deviance", "theta" , "ySim", "log_lik")
    if (log_lik == FALSE) {
      monitor <- monitor[-(length(monitor))]
    }
    inits <- lapply(1:chains, function(z) list("Intercept" = 0, .RNG.name = RNGlist[z], .RNG.seed = sample(1:10000, 1), "ySim" = y, "delta" = rep(1, P), "phi" = .20, "theta" = jitter(rep(0, P), amount = .25)))
    out <- run.jags(model = "jags_glm_spike.txt", modules = "glm", monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
    file.remove("jags_glm_spike.txt")
    if (!is.null(cl)) {
      parallel::stopCluster(cl = cl)
    }
    return(out)
  }
}