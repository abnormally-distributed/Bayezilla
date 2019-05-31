#' Adaptive Powered Correlation Prior for Gaussian likelihoods
#'
#' @param formula the model formula
#' @param data a data frame
#' @param lambda the power to use in the adaptive correlation prior. Default is -1, which gives the Zellner-Siow g prior. Setting
#' lambda to 0 results in a ridge-regression like prior. Setting lambda to a positive value adapts to collinearity by shrinking 
#' collinear variables towards a common value. Negative values of lambda pushes collinear variables further apart. I suggest 
#' fitting multiple values of lambda and selecting which is best via LOO-IC or WAIC. 
#' @param family one of "gaussian", "binomial", or "poisson".
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
#' @return
#' A run.jags object.
#' @export
#'
#' @examples
#' apcLm()
apcLm = function(formula, data, lambda = -1, log_lik = FALSE, iter=10000, warmup=1000, adapt=5000, chains=4, thin=3, method = "parallel", cl = makeCluster(2), ...)
{
  
  RNGlist = c("base::Wichmann-Hill", "base::Marsaglia-Multicarry", "base::Super-Duper", "base::Mersenne-Twister")
  
  if (chains > 4){
    chains = 4
  }
  
  data <- as.data.frame(data)
  y <- as.numeric(model.frame(formula, data)[, 1])
  X <- model.matrix(formula, data)[, -1]
  cormat = cor(X)
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
  
  apcLm  = " 
          model{
              tau ~ dgamma(.001, .001)
              g_inv ~ dgamma(.5, .5 * N)
              g <- 1 / g_inv
              for (j in 1:P){
                for (k in 1:P){
                   prior_scale[j,k] = g * (1/tau) * prior_cov[j,k]
                }
              }
              beta[1:P] ~ dmnorm.vcov(zeros[1:P], prior_scale[1:P,1:P])
              Intercept ~ dnorm(0, .01)            
              for (i in 1:N){
                 y[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 log_lik[i] <- logdensity.norm(y[i], Intercept + sum(beta[1:P] * X[i,1:P]), tau)
                 ySim[i] ~ dnorm(Intercept + sum(beta[1:P] * X[i,1:P]), tau)
              }
              sigma <- sqrt(1/tau)
              deviance <- -2 * sum(log_lik)
          }"
  
  write_lines(apcLm , "jags_apcLm.txt")
  jagsdata = list(X = X, y = y, N = length(y), P = ncol(X), prior_cov = prior_cov, zeros = rep(0, P))
  monitor = c("Intercept", "beta", "sigma", "g", "deviance", "ySim", "log_lik")
  if (log_lik == FALSE){
    monitor = monitor[-(length(monitor))]
  }
  inits = lapply(1:chains, function(z) list("Intercept" = 0, "beta" = jitter(rep(0, P), amount = .25), "g_inv" = .001, "tau" = 3, "ySim" = y,
                                            .RNG.name=RNGlist[z], .RNG.seed= sample(1:10000, 1)))
  out = run.jags(model = "jags_apcLm.txt", modules = "glm", monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, ...)
  if (!is.null(cl)){
    parallel::stopCluster(cl = cl)
  }
  file.remove("jags_apcLm.txt")
  return(out)
}
