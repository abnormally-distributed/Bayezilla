#' Extended Group Bayesian LASSO
#'
#' @description An adaptation of the extended Bayesian LASSO presented by Crispin M. Mutshinda and Mikko J. Sillanpää (2010)
#' for the purposes of group variable selection. The model specification is identical to the extended
#' Bayesian Lasso in all regards except the individual shrinkage hyperparameters & associated inclusion
#' indicators are replaced by group level parameters.
#'
#' Three versions of the model are presented. The default is "gamma" which places a gamma(.01, .01) prior
#' on the top-level shrinkage hyperparameter and independent gamma(1.007137, 0.7000579) priors on the
#' individual group shrinkage hyperparameters. This prior has a median of exactly 1, making the prior
#' inclusion probability 50%. This leads to a denominator in the Bayes Factor for each coefficient of .50/.50 = 1,
#' so that the posterior odds are equivalent to the Bayes Factor. In other words the posterior is fully
#' data driven. This prior also can be a little better with convergence and sampling effeciency on certain data sets.
#'
#' The second version is the "fixed_u" prior. While this retain the gamma(.01, .01) prior on the
#' top level shrinkage hyperparameter, the individaul shrinkage priors are of uniform(0, u) just as
#' in the original extended Bayesian LASSO. This requires you to choose an upper limit to the uniform
#' prior. Note that the prior inclusion probability is given by 1/u, so if you want a 50% inclusion
#' probability choose u=2. Common in Bayesian variable selection is to use a 20% probability if
#' dealing with a high dimensional problem, so for this choose u=5. If you have genuine prior
#' information you can and should use this to guide your choice. If you are unsure, use model comparison
#' to select which value of u to choose.
#'
#' The third version of the model is the original specification in Mutshinda and Sillanpää (2010), albeit with
#' the group selection modification, labeled
#' "classic" in the options here. The only differing feature from the above is that the top level shrinkage
#' hyperparameter is given a uniform(0, 100) prior. The advice for the "fixed_u" version applies here as well.
#'
#'
#' @param X the model matrix. Construct this manually with model.matrix()[,-1]
#' @param y the outcome variable
#' @param idx the group labels. Should be of length = to ncol(model.matrix()[,-1]) with the group assignments
#' for each covariate. Please ensure that you start numbering with 1, and not 0.
#' @param family one of "gaussian", "binomial", or "poisson".
#' @param eta_prior one of "gamma" (default), "fixed_u", or "classic".
#' @param fixed_u if using "fixed_u" or "classic" this must be assigned a value. It is empty by
#' default to force you to choose a good value through model comparison.
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
#' groupExtLASSO
groupExtLASSO <- function(X, y, idx, family = "normal", eta_prior = "gamma", fixed_u = NA, log_lik = FALSE, iter = 10000, warmup = 1000, adapt = 2000, chains = 4, thin = 3, method = "rjags", cl = NULL, ...) {
  if (family == "gaussian" || family == "normal") {
    if (eta_prior == "gamma") {
      jags_grp_extended_LASSO <- "model{

              # Precision
              tau ~ dgamma(.001, .001)

              # Shrinkage top-level-hyperparameter
              Omega ~ dgamma(.01, .01)

              Intercept ~ dnorm(0, .01)

              for (g in 1:nG){
                # Group Level shrinkage hyperparameter
                eta[g] ~ dgamma(1.007137, 0.7000579)
                lambda[g] <- Omega * eta[g]
                # Indicator Function
                delta[g] <- step(1-eta[g])
                w[g]<- pow(lambda[g],2)/2
              }

              for (p in 1:P){
                # Beta Precision
                beta_var[p] ~ dexp(w[idx[p]])
                beta_prec[p] <- 1 / beta_var[p]

                # Coefficient
                beta[p] ~ dnorm(0, beta_prec[p])

              }

              for (i in 1:N){
                 y[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 log_lik[i] <- logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 ySim[i] ~  dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
              }
              sigma <- sqrt(1/tau)
              deviance <- -2 * sum(log_lik[1:N])
          }"

      P <- ncol(X)
      nG <- length(unique(idx))
      monitor <- c("Intercept", "beta", "sigma", "Omega", "eta", "lambda", "delta", "deviance", "ySim", "log_lik")
      if (log_lik == FALSE) {
        monitor <- monitor[-(length(monitor))]
      }
      jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), idx = idx, nG = nG)
      inits <- lapply(1:chains, function(z) list("Intercept" = 0, "ySim" = y, "Omega" = 1, "beta" = jitter(rep(0, P), amount = .025), "eta" = rep(1, max(idx)), "beta_var" = abs(jitter(rep(.5, P), amount = .25)), "tau" = 1))
      write_lines(jags_grp_extended_LASSO, "jags_grp_extended_LASSO.txt")
    }

    if (eta_prior == "fixed_u") {
      if (is.na(fixed_u)) stop("Please select an upper limit for the uniform prior on eta.")

      jags_grp_extended_LASSO <- "model{

              # Precision
              tau ~ dgamma(.001, .001)

              # Shrinkage top-level-hyperparameter
              Omega ~ dgamma(.01, .01)

              Intercept ~ dnorm(0, .01)

              for (g in 1:nG){
                # Group Level shrinkage hyperparameter
                eta[g] ~ dunif(0, u)
                lambda[g] <- Omega * eta[g]
                # Indicator Function
                delta[g] <- step(1-eta[g])
                w[g]<- pow(lambda[g],2)/2
              }

              for (p in 1:P){
                # Beta Precision
                beta_var[p] ~ dexp(w[idx[p]])
                beta_prec[p] <- 1 / beta_var[p]

                # Coefficient
                beta[p] ~ dnorm(0, beta_prec[p])
              }

              for (i in 1:N){
                 y[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 log_lik[i] <- logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 ySim[i] ~  dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
              }
              sigma <- sqrt(1/tau)
              deviance <- -2 * sum(log_lik[1:N])
          }"

      monitor <- c("Intercept", "beta", "sigma", "Omega", "eta", "lambda", "delta", "deviance", "ySim", "log_lik")
      P <- ncol(X)
      nG <- length(unique(idx))
      jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), u = fixed_u, idx = idx, nG = nG)
      inits <- lapply(1:chains, function(z) list("Intercept" = 0, "ySim" = y, "Omega" = 1, "beta" = rnorm(P, 0, 1), "eta" = rep(1, nG), "beta_var" = abs(jitter(rep(.5, P), amount = .25)), "tau" = 1))
      write_lines(jags_grp_extended_LASSO, "jags_grp_extended_LASSO.txt")
    }


    if (eta_prior == "classic") {
      if (is.na(fixed_u)) stop("Please select an upper limit for the uniform prior on eta.")

      jags_grp_extended_LASSO <- "model{

              # Precision
              tau ~ dgamma(.001, .001)

              # Shrinkage top-level-hyperparameter
              Omega ~ dunif(0, 100)

              Intercept ~ dnorm(0, .01)

              for (g in 1:nG){
                # Group Level shrinkage hyperparameter
                eta[g] ~ dunif(0, u)
                lambda[g] <- Omega * eta[g]
                # Indicator Function
                delta[g] <- step(1-eta[g])
                w[g]<- pow(lambda[g],2)/2
              }

              for (p in 1:P){
                # Beta Precision
                beta_var[p] ~ dexp(w[idx[p]])
                beta_prec[p] <- 1 / beta_var[p]

                # Coefficient
                beta[p] ~ dnorm(0, beta_prec[p])

              }

              for (i in 1:N){
                 y[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 log_lik[i] <- logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 ySim[i] ~  dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
              }
              sigma <- sqrt(1/tau)
              deviance <- -2 * sum(log_lik[1:N])
          }"


      monitor <- c("Intercept", "beta", "sigma", "Omega", "eta", "lambda", "delta", "deviance", "ySim", "log_lik")
      if (log_lik == FALSE) {
        monitor <- monitor[-(length(monitor))]
      }
      P <- ncol(X)
      nG <- length(unique(idx))
      jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), u = fixed_u, idx = idx, nG = nG)
      inits <- lapply(1:chains, function(z) list("Intercept" = 0, "ySim" = y, "Omega" = 1, "beta" = rnorm(P, 0, 1), "eta" = rep(1, nG), "beta_var" = abs(jitter(rep(.5, P), amount = .25)), "tau" = 1))
      write_lines(jags_grp_extended_LASSO, "jags_grp_extended_LASSO.txt")
    }
  }

  if (family == "binomial" || family == "logistic") {
    if (eta_prior == "gamma") {
      jags_grp_extended_LASSO <- "model{

              # Shrinkage top-level-hyperparameter
              Omega ~ dgamma(.01, .01)

              Intercept ~ dnorm(0, .01)

              for (g in 1:nG){
                # Group Level shrinkage hyperparameter
                eta[g] ~ dgamma(1.007137, 0.7000579)
                lambda[g] <- Omega * eta[g]
                # Indicator Function
                delta[g] <- step(1-eta[g])
                w[g]<- pow(lambda[g],2)/2
              }

              for (p in 1:P){
                # Beta Precision
                beta_var[p] ~ dexp(w[idx[p]])
                beta_prec[p] <- 1 / beta_var[p]

                # Coefficient
                beta[p] ~ dnorm(0, beta_prec[p])

              }

              for (i in 1:N){
                 logit(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                 y[i] ~ dbern(psi[i])
                 log_lik[i] <- logdensity.bern(y[i], psi[i])
                 ySim[i] ~ dbern(psi[i])
              }
              deviance <- -2 * sum(log_lik[1:N])
          }"

      P <- ncol(X)
      nG <- length(unique(idx))
      monitor <- c("Intercept", "beta", "Omega", "eta", "lambda", "delta", "deviance", "ySim", "log_lik")
      if (log_lik == FALSE) {
        monitor <- monitor[-(length(monitor))]
      }
      jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), nG = nG)
      inits <- lapply(1:chains, function(z) list("Intercept" = 0, "ySim" = y, "Omega" = 1, "beta" = rnorm(P, 0, 1), "eta" = rep(1, nG), "beta_var" = abs(jitter(rep(.5, P), amount = .25))))
      write_lines(jags_grp_extended_LASSO, "jags_grp_extended_LASSO.txt")
    }

    else if (eta_prior == "fixed_u") {
      if (is.na(fixed_u)) stop("Please select an upper limit for the uniform prior on eta.")

      jags_grp_extended_LASSO <- "model{

              # Shrinkage top-level-hyperparameter
              Omega ~ dgamma(.01, .01)

              Intercept ~ dnorm(0, .01)
              for (g in 1:nG){
                # Group Level shrinkage hyperparameter
                eta[g] ~ dunif(0, u)
                lambda[g] <- Omega * eta[g]
                # Indicator Function
                delta[g] <- step(1-eta[g])
                w[g]<- pow(lambda[g],2)/2
              }

              for (p in 1:P){
                # Beta Precision
                beta_var[p] ~ dexp(w[idx[p]])
                beta_prec[p] <- 1 / beta_var[p]

                # Coefficient
                beta[p] ~ dnorm(0, beta_prec[p])

              }

              for (i in 1:N){
                 logit(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                  y[i] ~ dbern(psi[i])
                  log_lik[i] <- logdensity.bern(y[i], psi[i])
                  ySim[i] ~ dbern(psi[i])
              }
              deviance <- -2 * sum(log_lik[1:N])
            }"

      P <- ncol(X)
      monitor <- c("Intercept", "beta", "Omega", "eta", "lambda", "delta", "deviance", "ySim", "log_lik")
      if (log_lik == FALSE) {
        monitor <- monitor[-(length(monitor))]
      }
      jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), u = fixed_u, idx = idx, nG = max(idx))
      inits <- lapply(1:chains, function(z) list("Intercept" = 0, "ySim" = y, "Omega" = 1, "beta" = jitter(rep(0, P), amount = .25), "eta" = rep(1, P), "beta_var" = abs(jitter(rep(.5, P), amount = .25))))
      write_lines(jags_grp_extended_LASSO, "jags_grp_extended_LASSO.txt")
    }

    else if (eta_prior == "classic") {
      if (is.na(fixed_u)) stop("Please select an upper limit for the uniform prior on eta.")

      jags_grp_extended_LASSO <- "model{

              # Shrinkage top-level-hyperparameter
              Omega ~ dunif(0, 100)

              Intercept ~ dnorm(0, .01)

              for (g in 1:nG){
                # Group Level shrinkage hyperparameter
                eta[g] ~ dunif(0, u)
                lambda[g] <- Omega * eta[g]
                # Indicator Function
                delta[g] <- step(1-eta[g])
                w[g]<- pow(lambda[g],2)/2
              }

              for (p in 1:P){
                # Beta Precision
                beta_var[p] ~ dexp(w[idx[p]])
                beta_prec[p] <- 1 / beta_var[p]

                # Coefficient
                beta[p] ~ dnorm(0, beta_prec[p])
              }

              for (i in 1:N){
                 logit(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                  y[i] ~ dbern(psi[i])
                  log_lik[i] <- logdensity.bern(y[i], psi[i])
                  ySim[i] ~ dbern(psi[i])
              }
              deviance <- -2 * sum(log_lik[1:N])
            }"

      P <- ncol(X)
      monitor <- c("Intercept", "beta", "Omega", "eta", "lambda", "delta", "deviance", "ySim", "log_lik")
      if (log_lik == FALSE) {
        monitor <- monitor[-(length(monitor))]
      }
      jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), u = fixed_u, idx = idx, nG = max(idx))
      inits <- lapply(1:chains, function(z) list("Intercept" = 0, "ySim" = y, "Omega" = 1, "beta" = jitter(rep(0, P), amount = .25), "eta" = rep(1, P), "beta_var" = abs(jitter(rep(.5, P), amount = .25))))
      write_lines(jags_grp_extended_LASSO, "jags_grp_extended_LASSO.txt")
    }
  }

  else if (family == "poisson") {
    if (eta_prior == "gamma") {
      jags_grp_extended_LASSO <- "model{

              # Shrinkage top-level-hyperparameter
              Omega ~ dgamma(.01, .01)

              Intercept ~ dnorm(0, .01)

              for (g in 1:nG){
                # Group Level shrinkage hyperparameter
                eta[g] ~ dgamma(1.007137, 0.7000579)
                lambda[g] <- Omega * eta[g]
                # Indicator Function
                delta[g] <- step(1-eta[g])
                w[g]<- pow(lambda[g],2)/2
              }

              for (p in 1:P){
                # Beta Precision
                beta_var[p] ~ dexp(w[idx[p]])
                beta_prec[p] <- 1 / beta_var[p]

                # Coefficient
                beta[p] ~ dnorm(0, beta_prec[p])
              }

              for (i in 1:N){
                 log(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                 y[i] ~ dpois(psi[i])
                 log_lik[i] <- logdensity.pois(y[i], psi[i])
                 ySim[i] ~ dpois(psi[i])
              }
              deviance <- -2 * sum(log_lik[1:N])
          }"

      P <- ncol(X)
      monitor <- c("Intercept", "beta", "Omega", "eta", "lambda", "delta", "deviance", "ySim", "log_lik")
      if (log_lik == FALSE) {
        monitor <- monitor[-(length(monitor))]
      }
      jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X))
      inits <- lapply(1:chains, function(z) list("Intercept" = 0, "ySim" = y, "Omega" = 1, "beta" = jitter(0, amount = .025), "eta" = rep(1, P), "beta_var" = abs(jitter(rep(.5, P), amount = .25))))
      write_lines(jags_grp_extended_LASSO, "jags_grp_extended_LASSO.txt")
    }

    if (eta_prior == "fixed_u") {
      if (is.na(fixed_u)) stop("Please select an upper limit for the uniform prior on eta.")

      jags_grp_extended_LASSO <- "model{

              # Shrinkage top-level-hyperparameter
              Omega ~ dgamma(.01, .01)

              Intercept ~ dnorm(0, .01)

              for (g in 1:nG){
                # Group Level shrinkage hyperparameter
                eta[g] ~ dunif(0, u)
                lambda[g] <- Omega * eta[g]
                # Indicator Function
                delta[g] <- step(1-eta[g])
                w[g]<- pow(lambda[g],2)/2
              }

              for (p in 1:P){

                # Beta Precision
                beta_var[p] ~ dexp(w[idx[p]])
                beta_prec[p] <- 1 / beta_var[p]

                # Coefficient
                beta[p] ~ dnorm(0, beta_prec[p])
              }

              for (i in 1:N){
                 log(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                 y[i] ~ dpois(psi[i])
                 log_lik[i] <- logdensity.pois(y[i], psi[i])
                 ySim[i] ~ dpois(psi[i])
              }
              deviance <- -2 * sum(log_lik[1:N])
          }"

      P <- ncol(X)
      monitor <- c("Intercept", "beta", "Omega", "eta", "lambda", "delta", "deviance", "ySim", "log_lik")
      if (log_lik == FALSE) {
        monitor <- monitor[-(length(monitor))]
      }
      jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), u = fixed_u, idx = idx, nG = max(idx))
      inits <- lapply(1:chains, function(z) list("Intercept" = 0, "ySim" = y, "Omega" = 1, "beta" = jitter(0, amount = .025), "eta" = rep(1, P), "beta_var" = abs(jitter(rep(.5, P), amount = .25))))
      write_lines(jags_grp_extended_LASSO, "jags_grp_extended_LASSO.txt")
    }

    if (eta_prior == "classic") {
      if (is.na(fixed_u)) stop("Please select an upper limit for the uniform prior on eta.")

      jags_grp_extended_LASSO <- "model{

              # Shrinkage top-level-hyperparameter
              Omega ~ dunif(0, 100)

              Intercept ~ dnorm(0, .01)

              for (g in 1:nG){
                # Group Level shrinkage hyperparameter
                eta[g] ~ dunif(0, u)
                lambda[g] <- Omega * eta[g]
                # Indicator Function
                delta[g] <- step(1-eta[g])
                w[g]<- pow(lambda[g],2)/2
              }

              for (p in 1:P){
                # Beta Precision
                beta_var[p] ~ dexp(w[idx[p]])
                beta_prec[p] <- 1 / beta_var[p]

                # Coefficient
                beta[p] ~ dnorm(0, beta_prec[p])
              }

              for (i in 1:N){
                 log(psi[i]) <- Intercept + sum(beta[1:P] * X[i,1:P])
                 y[i] ~ dpois(psi[i])
                 log_lik[i] <- logdensity.pois(y[i], psi[i])
                 ySim[i] ~ dpois(psi[i])
              }
              deviance <- -2 * sum(log_lik[1:N])
          }"

      P <- ncol(X)
      monitor <- c("Intercept", "beta", "Omega", "eta", "lambda", "delta", "deviance", "ySim", "log_lik")
      if (log_lik == FALSE) {
        monitor <- monitor[-(length(monitor))]
      }
      jagsdata <- list(X = X, y = y, N = length(y), P = ncol(X), u = fixed_u, idx = idx, nG = max(idx))
      inits <- lapply(1:chains, function(z) list("Intercept" = 0, "Omega" = 1, "beta" = jitter(0, amount = .025), "eta" = rep(1, P), "beta_var" = abs(jitter(rep(.5, P), amount = .25))))
      write_lines(jags_grp_extended_LASSO, "jags_grp_extended_LASSO.txt")
    }
  }
  out <- run.jags(model = "jags_grp_extended_LASSO.txt", modules = "glm", monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, n.chains = chains, method = method, cl = cl, ...)
  file.remove("jags_grp_extended_LASSO.txt")
  if (!is.null(cl)) {
    parallel::stopCluster(cl = cl)
  }
  return(out)
}
