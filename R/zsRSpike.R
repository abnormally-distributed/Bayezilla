#' Robust Zellner-Siow g-prior Stochastic Search Variable Selection
#'
#' @description The same model as the \code{\link[Bayezilla]{zsSpike}} Gaussian model, but with the likelihood
#' function replaced by a "contaminated normal" model which consists of a convolution
#' of two Gaussians with the same mean but different variances. The Gaussian with the lower variance models the 
#' "uncontaminated" data while the Gaussian with the higher variance models the "contaminated"
#' portion of the data. Tukey, 1960; Box & Tiao, 1968; Gleason, 1993). The mixture proportions are obtained using either  Huber's
#' Psi, Hampel's Psi, or Tukey's Bisquare functions to generate weights. Information about the weight functions
#' is given in the figure below. (Huber, 1964; Hampel et al., 1986; Hoaglin, Mosteller, & Tukey, 1983) 
#'
#' The Fisher information matrix used as part of the g-prior (see the documentation in the \code{\link[Bayezilla]{zsSpike}}
#' function for more information) is calculated here using a minimum volume ellipsoid method to obtain a robust
#' covariance matrix, which is subjected to eigendecomposition and reconstructed as an inverse using the formula
#' L * D^-1 * L'. The matrix is then scaled according to the ratio of the trace of the non-robust Fisher information matrix
#' to the robust matrix. Note this particular model is meant to be an empirical Bayes or frequentist method estimated using
#' full Bayesian 'machinery' as a means to an end, and not so much a fully Bayesian model as most other functions in this
#' package. 
#'
#' \cr
#' \cr
#' \if{html}{\figure{robust.png}{}}
#' \if{latex}{\figure{robust.png}{}}
#' \cr
#'
#' @references 
#' Box, G., & Tiao, G. (1968). "A Bayesian Approach to Some Outlier Problems."" Biometrika, 55(1), 119-129. doi:10.2307/2334456 \cr
#' \cr
#' Gleason, J. R. (1993). "Understanding Elongation: The Scale Contaminated Normal Family" JASA 88(421) \cr
#' \cr
#' Hampel, F., E. Ronchetti, P. Rousseeuw, and W. Stahel (1986). "Robust Statistics: The Approach Based on Influence Functions." \cr
#' \cr
#' Hoaglin, D., Mosteller, F., & Tukey, J. (1983). "Understanding Robust and Exploratory Data Analysis." John Wiley and Sons, Inc., New York. \cr
#' \cr
#' Huber, Peter J. (1964). "Robust Estimation of a Location Parameter". Annals of Statistics. 53 (1): 73–101. doi:10.1214/aoms/1177703732 \cr
#' \cr
#' Liang, Paulo, Molina, Clyde, & Berger (2008). Mixtures of g Priors for Bayesian Variable Selection, Journal of the American Statistical Association, 103:481, 410-423, DOI: 10.1198/016214507000001337 \cr
#' \cr
#' Maronna, R. A., Martin, R. D., Yohai, V. J., & Salibian-Barrera, M. (2019). "Robust statistics: Theory and methods (with R)." Hoboken, NJ: John Wiley & Sons. \cr
#' \cr
#' Tukey, J. W. (1960). "A Survey of Sampling from Contaminated Distributions" in I. Olkin, ed., Contributions to Probability and Statistics \cr
#' \cr
#' Zellner, A. & Siow S. (1980). "Posterior odds ratio for selected regression hypotheses." In Bayesian statistics. Proc. 1st int. meeting (eds J. M. Bernardo, M. H. DeGroot, D. V. Lindley & A. F. M. Smith), 585?603. University Press, Valencia. \cr 
#' \cr
#' Zellner, A. (1986). "On assessing prior distributions and Bayesian regression analysis with g-prior distributions." In P. K. Goel and A. Zellner, editors, Bayesian Inference and Decision Techniques: Essays in Honor of Bruno de Finetti, 233?243.  \cr
#' \cr
#' @param formula the model formula
#' @param data a data frame
#' @param robfun "huber" for Huber's Psi, "tukey" for Tukey's Bisquare, or "hampel" for Hampel's Psi (the default).
#' @param c the tuning constant for the Huber weights function. Defaults to 1.345, which gives 95\% the efficiency of OLS when there are no outliers. 
#' @param t the tuning constant for the Tukey's bisquare weights function. Defaults to 4.685, which gives 95\% the efficiency of OLS when there are no outliers.
#' @param k the tuning constant controlling the asymptotic relative efficiency of Hampel's psi. The default is 0.9016085, which gives 95\% the efficiency of OLS when there are no outliers.
#' @param log_lik Should the log likelihood be monitored? The default is FALSE.
#' @param iter How many post-warmup samples? Defaults to 10000.
#' @param warmup How many warmup samples? Defaults to 1000.
#' @param adapt How many adaptation steps? Defaults to 2000.
#' @param chains How many chains? Defaults to 4. 
#' @param thin Thinning interval. Defaults to 1.
#' @param method Defaults to "rjparallel". For an alternative parallel option, choose "parallel". Otherwise, "rjags" (single core run).
#' @param cl Use parallel::makeCluster(# clusters) to specify clusters for the parallel methods. Defaults to two cores.
#' @param ... Other arguments to run.jags.
#'
#' @return
#' A run.jags object.
#' @export
#'
#' @examples
#' zsRSpike()
#' 
zsRSpike = function(formula, data, robfun = "hampel", c = 1.345, t = 4.685, k = 0.9016085, log_lik = FALSE, iter=10000, warmup=1000, adapt=5000, chains=4, thin=1, method = "rjparallel", cl = makeCluster(2), ...)
{
  
  robXtXinv = function (X) {
    
    data = MASS::mvrnorm(n = nrow(X), 
                         mu = colMeans(X), 
                         Sigma = as.matrix(corpcor::cov.shrink(X, verbose = FALSE)), 
                         empirical = TRUE)
    
    cor = MASS::cov.rob(data, cor = TRUE, method = "mve", nsamp = 2500)$cor
    cormat = cov2cor(fBasics::makePositiveDefinite(cor(X)))
    L = eigen(cormat)$vectors
    D = eigen(cormat)$values
    Trace = function(mat){sum(diag(mat))}
    P = ncol(X)
    Dpower =  pow(D, -1)
    xtxinv = L %*% diag(Dpower) %*% t(L)
    t = Trace(XtXinv(X))
    trace = Trace(xtxinv)
    K = t / trace
    xtxinv * K
  }
  
  huber.wts = function (r, c = 1.345){
    R <- abs(r)
    Rgtc <- R > c
    w <- r
    w[!Rgtc] <- 1
    w[Rgtc] <- c/R[Rgtc]
    return(w)
  }
  
  tukey.wts = function (r, t = 4.685){
    return((1 - pmin(1, abs(r/t))^2)^2)
  }
  
  
  hampel.wts = function(r, k = 0.9016085){
    
    psi.hampell = function(x, a = 1.5 * k, b = 2.8 * k, r = 6 * k){
      
      if (abs(x) <= a){
        x
      } else if (a < abs(x) & abs(x) <= b){
        a*sign(x)
      } else if (b < abs(x) & abs(x) <= r){
        a * sign(x) * ((r - abs(x)) / (r - b))
      } else if (r < abs(x)){
        0
      }
    }
    p  = sapply(r, function(x) psi.hampell(x) / x)
    if (any(is.nan(p))){
      wch = which(is.nan(p))
      p[wch] <- 1
    }
    return(p)
  }
  
  huberScale = function (y, k = 1.345, tol = 1e-06) {
    
    huberScale0 = function (y, k = 1.345, tol = 1e-06) 
    {
      y <- y[!is.na(y)]
      n <- length(y)
      mu <- median(y)
      s <- mad(y)
      if (s == 0) 
        s = 1e-6
      repeat {
        yy <- pmin(pmax(mu - k * s, y), mu + k * s)
        mu1 <- sum(yy)/n
        if (abs(mu - mu1) < tol * s) 
          break
        mu <- mu1
      }
      list(mu = mu, s = s)
    }
    
    
    huberScale1 = function (y, k = 1.345, mu, s, initmu = median(y), tol = 1e-06) 
    {
      mmu <- missing(mu)
      ms <- missing(s)
      y <- y[!is.na(y)]
      n <- length(y)
      if (mmu) {
        mu0 <- initmu
        n1 <- n - 1
      }
      else {
        mu0 <- mu
        mu1 <- mu
        n1 <- n
      }
      if (ms) {
        s0 <- mad(y)
        if (s0 == 0) 
          return(list(mu = mu0, s = 0))
      }
      else {
        s0 <- s
        s1 <- s
      }
      th <- 2 * pnorm(k) - 1
      beta <- th + k^2 * (1 - th) - 2 * k * dnorm(k)
      for (i in 1:30) {
        yy <- pmin(pmax(mu0 - k * s0, y), mu0 + k * s0)
        if (mmu) 
          mu1 <- sum(yy)/n
        if (ms) {
          ss <- sum((yy - mu1)^2)/n1
          s1 <- sqrt(ss/beta)
        }
        if ((abs(mu0 - mu1) < tol * s0) && abs(s0 - s1) < tol * 
            s0) 
          break
        mu0 <- mu1
        s0 <- s1
      }
      list(mu = mu0, s = s0)
    }
    
    
    estimate = huberScale0(y, k = k, tol = tol)
    estimate2 = huberScale1(y, k = k, initmu = estimate$mu, s = estimate$s, tol = tol)
    estimate3 = huberScale1(y, k = k, mu = estimate2$mu, tol = tol)
    estimate3$s
    
  }
  
  data <- as.data.frame(data)
  y <- as.numeric(model.frame(formula, data)[, 1])
  X <- as.matrix(model.matrix(formula, data)[,-1])
  prior_cov = robXtXinv(X)
  
  resids = y - as.vector(lmSolve(formula , data) %*% t(model.matrix(formula, data)))
  sigmaA = huberScale(resids)
  resids = (resids - mean(resids)) / sd(resids)
  
  if (robfun == "Tukey" || robfun == "tukey"){
    w = tukey.wts(resids, t = t)
  }
  
  if (robfun == "Huber" || robfun == "huber"){
    w = huber.wts(resids, c = c)
  }
  
  if (robfun == "Hampel" || robfun == "hampel"){
    w = hampel.wts(resids, k = k)
  }
  
  
  jags_apc = "model{
              phi ~ dbeta(1, 1)
              g_inv ~ dgamma(.5, N * .5)
              g <- 1 / g_inv
              
              for (j in 1:P){
                for (k in 1:P){
                  cov[j,k] = g * pow(sigmaA, 2) * prior_cov[j,k]
                }
              }
              
              omega <- inverse(cov) 
              theta[1:P] ~ dmnorm(rep(0,P), omega[1:P, 1:P])
              
              for (i in 1:P){
                delta[i] ~ dbern(phi)
                beta[i] <- delta[i] * theta[i]
              }
              
              Intercept ~ dnorm(0, 1e-10)
              for (i in 1:N){
                 mu[i] <- Intercept + sum(beta[1:P] * X[i, 1:P])
                 y[i] ~ dnormmix(rep(mu[i], 2), tau, c(w[i], 1-w[i]))
                 ySim[i] ~ dnormmix(rep(mu[i], 2), tau, c(w[i], 1 - w[i]))
                 log_lik[i] <- logdensity.normmix(y[i], rep(mu[i], 2), tau, c(w[i], 1 - w[i]))
                 residuals[i] <- pow(mu[i] - y[i], 2)
              }
               sigma <- sqrt(sum(residuals[1:N]) / (N - 1))
               Deviance <- -2 * sum(log_lik[1:N])
          }"
  
  P = ncol(X)
  write_lines(jags_apc, "jags_apc.txt")
  jagsdata = list(X = X, y = y,  N = length(y), P = ncol(X), prior_cov = prior_cov, w = sapply(w, function(x) max(c(x * 0.999, 0.001))), sigmaA = sigmaA, tau = c(1 / square(sigmaA), 0.01))
  monitor = c("Intercept", "beta", "g", "sigma", "phi", "delta", "Deviance", "ySim" ,"log_lik")
  if (log_lik == FALSE){
    monitor = monitor[-(length(monitor))]
  }
  inits = lapply(1:chains, function(z) list("Intercept" = coef(MASS::rlm(formula, data, method = "MM"))[1], 
                                            "theta" = coef(MASS::rlm(formula, data, method = "MM"))[-1], 
                                            "g_inv" = 1/length(y), 
                                            "ySim" = sample(y, length(y)), 
                                            "phi" = rbeta(1, 2, 2), 
                                            "delta" = sample(c(0, 1), replace = TRUE, size = P),
                                            .RNG.name= "lecuyer::RngStream", 
                                            .RNG.seed = sample(1:10000, 1)))
  
  
  
  out = run.jags(model = "jags_apc.txt", modules = c("mix on", "glm on", "dic off"), n.chains = chains, monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, factories = "mix::TemperedMix sampler off", ...)
  if (is.null(cl) == FALSE){
    parallel::stopCluster(cl = cl)
  }
  file.remove("jags_apc.txt")
  return(out)
  
}