#' Stochastic Search Variable Selection (Bernoulli-Normal Mixture)
#'
#' @description
#'
#' IMPORTANT NOTICE: This model works best on smaller to medium sized data sets with a small number of variables (less than 20).
#' If you experience difficulty with running times or obtaining a good effective sample size consider using the extended LASSO. \cr
#' \cr
#' IMPORTANT NOTICE: Center and scale your predictors before using this function.
#' If you do not scale and center your numeric predictors, this will likely not give reasonable results. 
#' \cr
#' \cr
#' This is the most basic type of Bayesian variable selection. This is inspired by the method presented by Kuo and Mallick (1998), 
#' with some improvements. This models the regression coefficients as coming from either a null distribution represented
#' by a probability mass of 100% at zero (the "spike") or unpenalized. The probability that a coefficient 
#' comes from the null-spike is controlled by a hyperparameter "phi" which estimates the overall probability of inclusion, 
#' i.e., the proportion of the P-number of predictors that are non-zero. This hyperparameter is given a Jeffrey's prior, 
#' beta(1/2, 1/2) which is non-informative and objective. The specification of the prior is pictured below. \cr
#' \cr
#' \if{html}{\figure{spike.png}{}}
#' \if{latex}{\figure{spike.png}{}}
#' \cr
#' Standard gaussian, binomial, and poisson likelihood functions are available. 
#' \cr 
#' For an alternative variable selection method, see the \code{\link[Bayezilla]{apcSpike}} and \code{\link[Bayezilla]{IAt}} 
#' functions. This package also implements several variants of the Bayesian LASSO.
#' \cr
#'
#' @references  
#' Kuo, L., & Mallick, B. (1998). Variable Selection for Regression Models. Sankhyā: The Indian Journal of Statistics, Series B, 60(1), 65-81.
#'
#'
#' @param formula the model formula
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
#' @return A run.jags object
#' @export
#'
#' @examples
#' glmSpike()
#' 
#' @seealso 
#' \code{\link[Bayezilla]{apcSpike}} 
#' \code{\link[Bayezilla]{extLASSO}}
#' \code{\link[Bayezilla]{negLASSO}}
#' \code{\link[Bayezilla]{bayesEnet}}
#' 
Spike <- function(formula, data, family = "gaussian", log_lik = FALSE, iter = 10000, warmup = 1000, adapt = 2000, chains = 4, thin = 1, method = "parallel", cl = makeCluster(2), ...) {
  X <- model.matrix(formula, data)[, -1]
  y <- model.frame(formula, data)[, 1]


  if (family == "gaussian") {
    
    jags_glm_spike <- "model{
              tau ~ dgamma(.01, .01)
              phi ~ dbeta(.5, .5) 
              
              for (p in 1:P){
                theta[p] ~ dnorm(0, 0.01)
                delta[p] ~ dbern(phi)
                beta[p] <- delta[p] * theta[p]
              }

              Intercept ~ dnorm(0, 1e-10)

              for (i in 1:N){
                 y[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 log_lik[i] <- logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 ySim[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
              }
              sigma <- sqrt(1/tau)
              Deviance <- -2 * sum(log_lik[1:N])
          }"

    P <- ncol(X)
    write_lines(jags_glm_spike, "jags_glm_spike.txt")
    jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X))
    monitor <- c("Intercept", "beta", "sigma", "phi", "Deviance",  "delta", "theta" ,"ySim", "log_lik")
    if (log_lik == FALSE) {
      monitor <- monitor[-(length(monitor))]
    }
    inits <- lapply(1:chains, function(z) list("Intercept" = 0, .RNG.name = "lecuyer::RngStream", .RNG.seed = sample(1:10000, 1),   "ySim" = y, "delta" = rep(1, P), "phi" = .20, "theta" = jitter(rep(0, P), amount = .25), "tau" = 1))
    out <- run.jags(model = "jags_glm_spike.txt", modules = c("glm on", "bugs on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl,summarise = FALSE, ...)
    return(out)
  }
  
  if (family == "binomial" || family == "logistic") {
    jags_glm_spike <- "model{
              phi ~ dbeta(.5, .5) 
              
              for (p in 1:P){
                theta[p] ~ dlogis(0, 0.01)
                delta[p] ~ dbern(phi)
                beta[p] <- delta[p] * theta[p]
              }
              
              Intercept ~ dnorm(0, 1e-10)
              
              for (i in 1:N){
                 logit(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                 y[i] ~ dbern(psi[i])
                 log_lik[i] <- logdensity.bern(y[i], psi[i])
                 ySim[i] ~ dbern(psi[i])
              }
             Deviance <- -2 * sum(log_lik[1:N])
          }"

    P <- ncol(X)
    write_lines(jags_glm_spike, "jags_glm_spike.txt")
    y <- as.numeric(as.factor(y)) - 1
    jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X))
    monitor <- c("Intercept", "beta", "phi", "Deviance", "delta",  "theta" , "ySim", "log_lik")
    if (log_lik == FALSE) {
      monitor <- monitor[-(length(monitor))]
    }
    inits <- lapply(1:chains, function(z) list("Intercept" = 0, .RNG.name = "lecuyer::RngStream", .RNG.seed = sample(1:10000, 1),   "ySim" = y, "delta" = rep(1, P), "phi" = .20, "theta" = jitter(rep(0, P), amount = .25)))
    out <- run.jags(model = "jags_glm_spike.txt", modules = c("glm on", "bugs on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE,...)
    return(out)
  }

  if (family == "poisson") {
    jags_glm_spike <- "model{
              phi ~ dbeta(.5, .5) 
              
              for (p in 1:P){
                theta[p] ~ dnorm(0, 0.01)
                delta[p] ~ dbern(phi)
                beta[p] <- delta[p] * theta[p]
              }
              
              Intercept ~ dnorm(0, 1e-10)
              
              for (i in 1:N){
                 log(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                 y[i] ~ dpois(psi[i])
                 log_lik[i] <- logdensity.pois(y[i], psi[i])
                 ySim[i] ~ dpois(psi[i])
              }
              Deviance <- -2 * sum(log_lik[1:N])
          }"

    write_lines(jags_glm_spike, "jags_glm_spike.txt")
    P <- ncol(X)
    jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X))
    monitor <- c("Intercept", "beta", "phi", "Deviance", "delta",  "theta" , "ySim", "log_lik")
    if (log_lik == FALSE) {
      monitor <- monitor[-(length(monitor))]
    }
    inits <- lapply(1:chains, function(z) list("Intercept" = 0, .RNG.name = "lecuyer::RngStream", .RNG.seed = sample(1:10000, 1), "ySim" = y, "delta" = rep(1, P), "phi" = .20, "theta" = jitter(rep(0, P), amount = .25)))
    out <- run.jags(model = "jags_glm_spike.txt", modules = c("glm on", "bugs on", "dic off"), monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
    file.remove("jags_glm_spike.txt")
    if (!is.null(cl)) {
      parallel::stopCluster(cl = cl)
    }
    return(out)
  }
}