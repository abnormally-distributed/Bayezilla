apcLmBetaMix = function(formula, data, phi_prior = .5 , lambda = -1, log_lik = FALSE, iter=10000, warmup=1000, adapt=5000, chains=4, thin=3,method = "rjags", cl = NULL, ...)
{
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
              a ~ dunif(1, 10)
              b ~ dunif(1, 10)
              for (p in 1:P){
                delta[p] ~ dbeta(a, b)
              }
              for (j in 1:P){
                for (k in 1:P){
                   prior_scale[j,k] = g * (1/tau) * prior_cov[j,k]
                }
              }
              theta[1:P] ~ dmnorm.vcov(zeros[1:P], prior_scale[1:P,1:P])
              for (i in 1:P){
                  beta[i] <- delta[i] * theta[i]
              }
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
  monitor = c("Intercept", "beta", "sigma", "g", "deviance", "a", "b", "theta", "delta", "ySim", "log_lik")
  if (log_lik == FALSE){
    monitor = monitor[-(length(monitor))]
  }
  inits = lapply(1:chains, function(z) list("Intercept" = 0, "a" =1 , "b" = 1, "delta" = sample(seq(0, 1, by = .01), P), "theta" = jitter(rep(0, P), amount = .25), "g_inv" = .001, "tau" = 3, "ySim" = y))
  out = run.jags(model = "jags_apcLm.txt", modules = "glm", monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, ...)
  if (!is.null(cl)){
    parallel::stopCluster(cl = cl)
  }
  file.remove("jags_apcLm.txt")
  return(out)
}
