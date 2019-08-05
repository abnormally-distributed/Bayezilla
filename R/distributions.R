
###########################################################################
#########  Power Exponential (Generalized Normal) Distribution)  ##########
###########################################################################


#' Power Exponential random number generator
#'
#' @param n number of observations
#' @param mu vector of location parameter values
#' @param sigma vector of scale parameter values
#' @param kappa vector of shape parameter values
#'
#' @return vector
#' @export
#'
#' @examples
#' rpowexp()
rpowexp = function (n, mu = 0, sigma = 1, kappa = 2) 
{
  if (any(sigma <= 0)) 
    stop(paste("the sigma parameter must be positive", "\n", ""))
  if (any(n <= 0)) 
    stop(paste("n must be a positive integer", "\n", ""))
  n <- ceiling(n)
  p <- runif(n)
  r <- qpowexp(p, mu = mu, sigma = sigma, kappa = kappa)
  r
}


#' Power Exponential quantile function
#'
#' @param p vector of probabilities.
#' @param mu vector of location parameter values
#' @param sigma vector of scale parameter values
#' @param kappa vector of shape parameter values
#' @param lower.tail if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param log if TRUE, probabilities p are given as log(p).
#'
#' @return vector
#' @export
#'
#' @examples
#' qpowexp()
qpowexp = function (p, mu = 0, sigma = 1, kappa = 2, lower.tail = TRUE, log = FALSE) 
{
  if (any(sigma < 0)) 
    stop(paste("the sigma parameter must be positive", "\n", ""))
  if (any(kappa < 0)) 
    stop(paste("the kappa parameter must be positive", "\n", ""))
  if (log == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  log.c <- 0.5 * (-(2/kappa) * log(2) + lgamma(1/kappa) - lgamma(3/kappa))
  c <- exp(log.c)
  suppressWarnings(s <- qgamma((2 * p - 1) * sign(p - 0.5), 
                               shape = (1/kappa), scale = 1))
  z <- sign(p - 0.5) * ((2 * s)^(1/kappa)) * c
  ya <- mu + sigma * z
  ya
}

#' Power Exponential probability density function
#'
#' @param x vector of quantiles
#' @param mu vector of location parameter values
#' @param sigma vector of scale parameter values
#' @param kappa vector of shape parameter values
#' @param log if TRUE, probabilities p are given as log(p).
#'
#' @return vector
#' @export
#'
#' @examples
#' dpowexp()
dpowexp = function (x, mu = 0, sigma = 1, kappa = 2, log = FALSE) {
  if (any(sigma < 0)) 
    stop(paste("the sigma parameter must be positive", "\n", ""))
  if (any(kappa < 0)) 
    stop(paste("the kappa parameter must be positive", "\n", ""))
  log.c <- 0.5 * (-(2/kappa) * log(2) + lgamma(1/kappa) - lgamma(3/kappa))
  c <- exp(log.c)
  z <- (x - mu)/sigma
  log.lik <- -log(sigma) + log(kappa) - log.c - (0.5 * (abs(z/c)^kappa)) - 
    (1 + (1/kappa)) * log(2) - lgamma(1/kappa)
  if (log == FALSE) 
    fy <- exp(log.lik)
  else fy <- log.lik
  fy
}



#' Power Exponential cumulative probability function
#'
#' @param q vector of quantiles
#' @param mu vector of location parameter values
#' @param sigma vector of scale parameter values
#' @param kappa vector of shape parameter values
#' @param lower.tail if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param log if TRUE, probabilities p are given as log(p).
#'
#' @return vector
#' @export
#'
#' @examples
#' ppowexp()
ppowexp = function (q, mu = 0, sigma = 1, kappa = 2, lower.tail = TRUE, log = FALSE) {
  if (any(sigma < 0)) 
    stop(paste("the sigma parameter must be positive", "\n", ""))
  if (any(kappa < 0)) 
    stop(paste("the kappa parameter must be positive", "\n", ""))
  log.c <- 0.5 * (-(2/kappa) * log(2) + lgamma(1/kappa) - lgamma(3/kappa))
  c <- exp(log.c)
  z <- (q - mu)/sigma
  s <- 0.5 * ((abs(z/c))^kappa)
  cdf <- 0.5 * (1 + pgamma(s, shape = 1/kappa, scale = 1) * sign(z))
  if (length(kappa) > 1) 
    cdf <- ifelse(kappa > 10000, (q - (mu - sqrt(3) * sigma))/(sqrt(12) * sigma), cdf)
  else cdf <- if (kappa > 10000) 
    (q - (mu - sqrt(3) * sigma))/(sqrt(12) * sigma)
  else cdf
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}


###########################################################################
#########             Student's T Distribution                   ##########
###########################################################################

#' Student-T random number generator
#'
#' @param n number of observations
#' @param nu vector of normality parameter values
#' @param mu vector of location parameter values
#' @param scale vector of scale parameter values
#'
#' @return vector
#' @export
#'
#' @examples
#' rst()
rst = function (n, nu = 3, mu = 0, scale = 1) {
  
  if (any(scale <= 0)) 
    stop("the scale parameter must be positive.")
  if (any(nu <= 0)) 
    stop("the nu parameter must be positive.")
  
  rGamma = function(n, shape = 1, rate = 1) {
    return(qgamma(seq(1/n, 1 - 1/n, length.out = n), 
                  shape, rate))
  }
  
  gamma.samps = rGamma(n, nu * 0.50, (scale^2) * (nu * 0.50))
  sigmas = sqrt(1 / gamma.samps)
  x = rnorm(n, mu, sigmas)
  return(x)
  
}

#' Student-T probability density function
#'
#' @param x vector of quantiles
#' @param nu vector of normality parameter values
#' @param mu vector of location parameter values
#' @param scale vector of scale parameter values
#' @param log if TRUE, probabilities p are given as log(p).
#'
#' @return vector
#' @export
#'
#' @examples
#' dst()
dst = function (x, nu = 3, mu = 0, scale = 1, log = FALSE) 
{
  if (any(scale <= 0)) {
    stop("the scale parameter must be positive.")
  }
  if (any(nu <= 0))  {
    stop("the nu parameter must be positive.")
  }
  if (log) {
    dt((x - mu)/scale, df = nu, log = TRUE) - log(scale)
  }
  else {
    dt((x - mu)/scale, df = nu)/scale
  }
}


#' Student-T quantile function
#'
#' @param p vector of probabilities.
#' @param nu vector of normality parameter values
#' @param mu vector of location parameter values
#' @param scale vector of scale parameter values
#' @param lower.tail if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param log if TRUE, probabilities p are given as log(p).
#'
#' @return vector
#' @export
#'
#' @examples
#' qst()
qst = function (p, nu = 3, mu = 0, scale = 1, lower.tail = TRUE, log = FALSE) 
{
  if (any(p < 0) || any(p > 1)) 
    stop("p must be in [0,1].")
  if (any(scale <= 0)) 
    stop("the scale parameter must be positive.")
  if (any(nu <= 0)) 
    stop("the nu parameter must be positive.")
  NN <- max(length(p), length(mu), length(scale), length(nu))
  p <- rep(p, len = NN)
  mu <- rep(mu, len = NN)
  scale <- rep(scale, len = NN)
  nu <- rep(nu, len = NN)
  q <- mu + scale * qt(p, df = nu, lower.tail = lower.tail)
  temp <- which(nu > 1e+06)
  q[temp] <- qnorm(p[temp], mu[temp], scale[temp], lower.tail = lower.tail, 
                   log.p = log)
  return(q)
}

#' Student-T cumulative probability function
#'
#' @param q vector of quantiles
#' @param kappa vector of normality parameter values
#' @param mu vector of location parameter values
#' @param scale vector of scale parameter values
#' @param lower.tail if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param log if TRUE, probabilities p are given as log(p).
#'
#' @return vector
#' @export
#'
#' @examples
#' pst()
pst = function (q, nu = 3, mu = 0, scale = 1, lower.tail = TRUE, log = FALSE) 
{
  if (any(scale <= 0)) {
    stop("the scale parameter must be positive.")
  }
  if (any(nu <= 0))  {
    stop("the nu parameter must be positive.")
  }
  pt((q - mu)/scale, df = nu, lower.tail = lower.tail, log.p = log)
}

###########################################################################
#########             Asymmetric Laplace Distribution            ##########
###########################################################################


#' Asymmetric Laplace random number generator
#'
#' @param n number of observations
#' @param mu vector of location parameter values
#' @param scale vector of scale parameter values
#' @param q vector of quantiles
#'
#' @return vector
#' @export
#'
#' @examples
#' ral()
ral  = function (n, mu = 0, scale = 1, k = 0.5) 
{
  if (any(scale <= 0)) 
    stop("the scale parameter must be positive.")
  u <- runif(n)
  qal(u, mu = mu, scale = scale, k = k)
}

#' Asymmetric Laplace quantile function
#'
#' @param p vector of probabilities.
#' @param mu vector of location parameter values
#' @param scale vector of scale parameter values
#' @param k vector of quantile locations
#' @param lower.tail if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param log.p if TRUE, probabilities p are given as log(p).
#'
#' @return vector
#' @export
#'
#' @examples
#' qal()
qal  = function (p, mu = 0, scale = 1, k = 0.5, lower.tail = TRUE, log.p = FALSE) 
{
  if (any(scale <= 0)) 
    stop("the scale parameter must be positive.")
  
  if (log.p) {
    p <- exp(p)
  }
  if (!lower.tail) {
    p <- 1 - p
  }
  if (length(k) == 1L) {
    k <- rep(k, length(mu))
  }
  ifelse(p < k, yes = mu + ((scale * log(p/k))/(1 - k)), no = mu - ((scale * log((1 - p)/(1 - k)))/k))
}


#' Asymmetric Laplace cumulative probability function
#'
#' @param q vector of quantiles
#' @param mu vector of location parameter values
#' @param scale vector of scale parameter values
#' @param k vector of quantile locations
#' @param lower.tail if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param log.p if TRUE, probabilities p are given as log(p).
#'
#' @return vector
#' @export
#'
#' @examples
#' pal()
pal  = function (q, mu = 0, scale = 1, k = 0.5, lower.tail = TRUE, 
                 log.p = FALSE) 
{
  if (any(scale <= 0)) 
    stop("the scale parameter must be positive.")
  
  out <- ifelse(q < mu, k * exp((1 - k) * (q - mu)/scale), 1 - (1 - k) * exp(-k * (q - mu)/scale))
  if (!lower.tail){
    out <- 1 - out
  }
  if (log.p) {
    out <- log(out)
  }
  out
}


#' Asymmetric Laplace probability density function
#'
#' @param x vector of observations
#' @param mu vector of location parameter values
#' @param scale vector of scale parameter values
#' @param k vector of quantile locations
#' @param log if TRUE, probabilities p are given as log(p).
#'
#' @return vector
#' @export
#'
#' @examples
#' dal()
dal  = function (x, mu = 0, scale = 1, k = 0.5, log = FALSE) 
{
  
  if (any(scale <= 0)) 
    stop("the scale parameter must be positive.")
  
  out <- ifelse(x < mu, yes = (k * (1 - k)/scale) * 
                  exp((1 - k) * (x - mu)/scale), no = (k * (1 - k)/scale) * exp(-k * (x - mu)/scale))
  if (log) {
    out <- log(out)
  }
  out
}


###########################################################################
#########           Sinh-Arcsinh (SHASH) Distribution            ##########
###########################################################################

#' Sinh-Arcsinh (SHASH) random number generator
#'
#' @param n number of observations
#' @param mu vector of location parameter values
#' @param scale vector of scale parameter values
#' @param alpha vector of skewness parameter
#' @param kappa vector of kurtosis parameter
#'
#' @return vector
#' @export
#'
#' @examples
#' rshash()
rshash = function (n, mu = 0, scale = 1, alpha = 0, kappa = 1) 
{
  if (any(n <= 0)) 
    stop(paste("n must be a positive integer", "\n", ""))
  n <- ceiling(n)
  p <- runif(n)
  r <- qshash(p, mu = mu, scale = scale, alpha = alpha, kappa = kappa)
  r
}

#' Sinh-Arcsinh (SHASH) probability density function
#'
#' @param x vector of quantiles
#' @param mu vector of location parameter values
#' @param scale vector of scale parameter values
#' @param alpha vector of skewness parameter
#' @param kappa vector of kurtosis parameter
#' @param log if TRUE, probabilities p are given as log(p).
#'
#' @return vector
#' @export
#'
#' @examples
#' dshash()
dshash = function (x, mu = 0, scale = 1, alpha = 0, kappa = 1, log = FALSE){
  
  if (any(scale < 0)) 
    stop(paste("scale must be positive", "\n", ""))
  if (any(kappa < 0)) 
    stop(paste("kappa must be positive", "\n", ""))
  z <- (x - mu)/(scale * kappa)
  c <- cosh(kappa * asinh(z) - alpha)
  r <- sinh(kappa * asinh(z) - alpha)
  loglik <- -log(scale) - 0.5 * log(2 * pi) - 0.5 * log(1 +  (z^2)) + log(c) - 0.5 * (r^2)
  if (log == FALSE) 
    fy <- exp(loglik)
  else fy <- loglik
  fy
}

#' Sinh-Arcsinh (SHASH) cumulative probability function
#'
#' @param v vector of quantiles
#' @param mu vector of location parameter values
#' @param scale vector of scale parameter values
#' @param alpha vector of skewness parameter
#' @param kappa vector of kurtosis parameter
#' @param lower.tail if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param log.p if TRUE, probabilities p are given as log(p).
#'
#' @return vector
#' @export
#'
#' @examples
#' pshash()
pshash = function (q, mu = 0, scale = 1, alpha = 0, kappa = 1, lower.tail = TRUE, 
                   log.p = FALSE) 
{
  if (any(scale < 0)) 
    stop(paste("scale must be positive", "\n", ""))
  if (any(kappa < 0)) 
    stop(paste("kappa must be positive", "\n", ""))
  z <- (q - mu)/(scale * kappa)
  c <- cosh(kappa * asinh(z) - alpha)
  r <- sinh(kappa * asinh(z) - alpha)
  p <- pNO(r)
  if (lower.tail == TRUE) 
    p <- p
  else p <- 1 - p
  if (log.p == FALSE) 
    p <- p
  else p <- log(p)
  p
}

#' Sinh-Arcsinh (SHASH) quantile function
#'
#' @param p vector of probabilities.
#' @param mu vector of location parameter values
#' @param scale vector of scale parameter values
#' @param alpha vector of skewness parameter
#' @param kappa vector of kurtosis parameter
#' @param lower.tail if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param log.p if TRUE, probabilities p are given as log(p).
#'
#' @return vector
#' @export
#'
#' @examples
#' qshash()
qshash = function (p, mu = 0, scale = 1, alpha = 0, kappa = 1, lower.tail = TRUE, log.p = FALSE) 
{
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (any(p <= 0) | any(p >= 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  if (lower.tail == TRUE) 
    p <- p
  else p <- 1 - p
  y <- mu + (kappa * scale) * sinh((1/kappa) * asinh(qnorm(p)) + 
                                     (alpha/kappa))
  y
}

