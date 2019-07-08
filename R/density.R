#' Modified version of density() with additional kernels
#'
#' @description This modifies the base R density() function to add three new kernels - the Laplace
#' , Student-T, and Huber (with fixed k = 1.58) densities. The default number of points in the estimated density curve have
#' also been changed from the paltry 512 to 2048 in order to give better density plots.
#'
#' @param x the data from which the estimate is to be computed. For the default method 
#' a numeric vector: long vectors are not supported.
#' @param bw the smoothing bandwidth to be used. The kernels are scaled such that 
#' this is the standard deviation of the smoothing kernel. (Note this differs from 
#' the reference books cited below, and from S-PLUS.) \cr \cr
#' bw can also be a character string giving a rule to choose the bandwidth. 
#' See bw.nrd. \cr \cr "SJ" is the default. "nrd0" has remained the default in the unmodified version of this function
#' for historical reasons, rather than as a the best method, where "SJ" better fits. Hence this package takes the progressive path
#' and just uses the better option as the default.
#' @param adjust The specified (or computed) value of bw is multiplied by adjust.
#' @param df The degrees of freedom for the student t kernel.
#' @param kernel a character string giving the smoothing kernel to be used. 
#' This must partially match one of "normal", "laplacian", "student", "huber", "rectangular", "triangular", "epanechnikov", "biweight",
#' "cosine" or "optcosine", with default "normal", and may be abbreviated to a unique prefix (single letter). \cr \cr
#' "cosine" is smoother than "optcosine", which is the usual ‘cosine’ kernel in the literature and almost MSE-efficient.
#' However, "cosine" is the version used by S. 
#' @param weights numeric vector of non-negative observation weights, hence of same length as x. 
#' The default NULL is equivalent to weights = rep(1/nx, nx) where nx is the length of (the finite entries of) x[].
#' @param window a character string giving the smoothing kernel to be used. 
#' This must partially match one of "normal", "laplacian", "student", "huber", "rectangular", "triangular", "epanechnikov", "biweight",
#' "cosine" or "optcosine", with default "normal", and may be abbreviated to a unique prefix (single letter). \cr \cr
#' "cosine" is smoother than "optcosine", which is the usual ‘cosine’ kernel in the literature and almost MSE-efficient.
#' However, "cosine" is the version used by S. 
#' @param width this exists for compatibility with S; if given, and bw is not, will set bw to width if this is a character string, or to a kernel-dependent multiple of width if this is numeric.
#' @param give.Rkern logical; if true, no density is estimated, and the ‘canonical bandwidth’ of the chosen kernel is returned instead.
#' @param n the number of equally spaced points at which the density is to be estimated. When n > 512, it is rounded up to a power of 2 during the calculations (as fft is used) and the final result is interpolated by approx. So it almost always makes sense to specify n as a power of two.
#' @param from the left and right-most points of the grid at which the density is to be estimated; the defaults are cut * bw outside of range(x).
#' @param to the left and right-most points of the grid at which the density is to be estimated; the defaults are cut * bw outside of range(x).
#' @param cut by default, the values of from and to are cut bandwidths beyond the extremes of the data. This allows the estimated density to drop to approximately zero at the extremes.
#' @param na.rm logical; if TRUE, missing values are removed from x. If FALSE any missing values cause an error.
#' @param ... further arguments for (non-default) methods.
#'
#' @return If give.Rkern is true, the number R(K), otherwise an object with class "density" whose underlying structure is a list containing the following components.
#' @export
#'
#' @examples
#' require(graphics)
#' plot(density(c(-20, rep(0,98), 20)), xlim = c(-4, 4))  # IQR = 0 
#' # The Old Faithful geyser data 
#' d <- density(faithful$eruptions, bw = "sj")
#' d
#' plot(d)
#' plot(d, type = "n")
#' polygon(d, col = "wheat")
#' 
#' 
density <- function (x, 
                     bw = "sj", 
                     adjust = 1, 
                     df = 6, 
                     kernel = c("normal", 
                                "laplacian", 
                                "laplace" , 
                                "huber",
                                "student", 
                                "epanechnikov", 
                                "rectangular", 
                                "triangular", 
                                "biweight", 
                                "cosine", 
                                "optcosine"), 
                     weights = NULL, window = kernel, width, give.Rkern = FALSE, n = 2048, from, to, cut = 3, na.rm = FALSE, ...) {
  
        UseMethod("density")
}

#' Modified version of density() with laplacian and student-t kernels
#'
#' @description This modifies the base R density() function to add three new kernels - the Laplace
#' , Student-T, and Huber (with fixed k = 1.58) densities. The default number of points in the estimated density curve have
#' also been changed from the paltry 512 to 2048 in order to give better density plots.
#'
#' @param x the data from which the estimate is to be computed. For the default method 
#' a numeric vector: long vectors are not supported.
#' @param bw the smoothing bandwidth to be used. The kernels are scaled such that 
#' this is the standard deviation of the smoothing kernel. (Note this differs from 
#' the reference books cited below, and from S-PLUS.) \cr \cr
#' bw can also be a character string giving a rule to choose the bandwidth. 
#' See bw.nrd. \cr \cr "SJ" is the default. "nrd0" has remained the default in the unmodified version of this function
#' for historical reasons, rather than as a the best method, where "SJ" better fits. Hence this package takes the progressive path
#' and just uses the better option as the default.
#' @param adjust The specified (or computed) value of bw is multiplied by adjust.
#' @param df The degrees of freedom for the student t kernel. 
#' @param kernel a character string giving the smoothing kernel to be used. 
#' This must partially match one of "normal", "laplacian", "student", "huber", "rectangular", "triangular", "epanechnikov", "biweight",
#' "cosine" or "optcosine", with default "normal", and may be abbreviated to a unique prefix (single letter). \cr \cr
#' "cosine" is smoother than "optcosine", which is the usual ‘cosine’ kernel in the literature and almost MSE-efficient.
#' However, "cosine" is the version used by S. 
#' @param weights numeric vector of non-negative observation weights, hence of same length as x. 
#' The default NULL is equivalent to weights = rep(1/nx, nx) where nx is the length of (the finite entries of) x[].
#' @param window a character string giving the smoothing kernel to be used. 
#' This must partially match one of "normal", "laplacian", "student", "huber","rectangular", "triangular", "epanechnikov", "biweight",
#' "cosine" or "optcosine", with default "normal", and may be abbreviated to a unique prefix (single letter). \cr \cr
#' "cosine" is smoother than "optcosine", which is the usual ‘cosine’ kernel in the literature and almost MSE-efficient.
#' However, "cosine" is the version used by S. 
#' @param width this exists for compatibility with S; if given, and bw is not, will set bw to width if this is a character string, or to a kernel-dependent multiple of width if this is numeric.
#' @param give.Rkern logical; if true, no density is estimated, and the ‘canonical bandwidth’ of the chosen kernel is returned instead.
#' @param n the number of equally spaced points at which the density is to be estimated. When n > 512, it is rounded up to a power of 2 during the calculations (as fft is used) and the final result is interpolated by approx. So it almost always makes sense to specify n as a power of two.
#' @param from the left and right-most points of the grid at which the density is to be estimated; the defaults are cut * bw outside of range(x).
#' @param to the left and right-most points of the grid at which the density is to be estimated; the defaults are cut * bw outside of range(x).
#' @param cut by default, the values of from and to are cut bandwidths beyond the extremes of the data. This allows the estimated density to drop to approximately zero at the extremes.
#' @param na.rm logical; if TRUE, missing values are removed from x. If FALSE any missing values cause an error.
#' @param ... further arguments for (non-default) methods.
#'
#' @return If give.Rkern is true, the number R(K), otherwise an object with class "density" whose underlying structure is a list containing the following components.
#' @export
#'
#' @examples
#' require(graphics)
#' plot(density(c(-20, rep(0,98), 20)), xlim = c(-4, 4))  # IQR = 0 
#' # The Old Faithful geyser data 
#' d <- density(faithful$eruptions, bw = "sj")
#' d
#' plot(d)
#' plot(d, type = "n")
#' polygon(d, col = "wheat")
#' 
#' 
density.default <- function (x, bw = "sj", adjust = 1, df = 6, kernel = c("normal", "laplacian", "laplace" , "student", "huber",
                                                  "epanechnikov", "rectangular", "triangular", "biweight", 
                                                 "cosine", "optcosine"), weights = NULL, window = kernel, 
          width, give.Rkern = FALSE, n = 2048, from, to, cut = 3, na.rm = FALSE, 
          ...) 
{

  chkDots(...)
  if (!missing(window) && missing(kernel)) 
    kernel <- window
  kernel <- match.arg(kernel)
  if (give.Rkern) 
    return(switch(kernel, 
                  normal = 1/(2 * sqrt(pi)),
                  huber = 1/(1.5 * sqrt(pi)),
                  student = 1/(sqrt(((df-2) / df) * pi)), 
                  laplacian = 1/(2 * pi), 
                  laplace = 1/(2 * pi), 
                  rectangular = sqrt(3)/6, 
                  triangular = sqrt(6)/9, 
                  epanechnikov = 3/(5 * sqrt(5)), 
                  biweight = 5 * sqrt(7)/49, 
                  cosine = 3/4 * sqrt(1/3 - 2/pi^2), 
                  optcosine = sqrt(1 - 8/pi^2) * pi^2/16))
  
  if (!is.numeric(x)) 
    stop("argument 'x' must be numeric")
  name <- deparse(substitute(x))
  x <- as.vector(x)
  x.na <- is.na(x)
  if (any(x.na)) {
    if (na.rm) 
      x <- x[!x.na]
    else stop("'x' contains missing values")
  }
  N <- nx <- as.integer(length(x))
  if (is.na(N)) 
    stop(gettextf("invalid value of %s", "length(x)"), domain = NA)
  x.finite <- is.finite(x)
  if (any(!x.finite)) {
    x <- x[x.finite]
    nx <- length(x)
  }
  if (is.null(weights)) {
    weights <- rep.int(1/nx, nx)
    totMass <- nx/N
  }
  else {
    if (length(weights) != N) 
      stop("'x' and 'weights' have unequal length")
    if (!all(is.finite(weights))) 
      stop("'weights' must all be finite")
    if (any(weights < 0)) 
      stop("'weights' must not be negative")
    wsum <- sum(weights)
    if (any(!x.finite)) {
      weights <- weights[x.finite]
      totMass <- sum(weights)/wsum
    }
    else totMass <- 1
    if (!isTRUE(all.equal(1, wsum))) 
      warning("sum(weights) != 1  -- will not get true density")
  }
  n.user <- n
  n <- max(n, 512)
  if (n > 512) 
    n <- 2^ceiling(log2(n))
  if (missing(bw) && !missing(width)) {
    if (is.numeric(width)) {
      fac <- switch(kernel, normal = 4, huber = 4, student = 4, laplace = 4, laplacian = 4, rectangular = 2 * 
                      sqrt(3), triangular = 2 * sqrt(6), epanechnikov = 2 * 
                      sqrt(5), biweight = 2 * sqrt(7), cosine = 2/sqrt(1/3 -  2/pi^2), optcosine = 2/sqrt(1 - 8/pi^2))
      bw <- width/fac
    }
    if (is.character(width)) 
      bw <- width
  }
  if (is.character(bw)) {
    if (nx < 2) 
      stop("need at least 2 points to select a bandwidth automatically")
    bw <- switch(tolower(bw), 
                 nrd0 = bw.nrd0(x), 
                 nrd = bw.nrd(x), 
                 ucv = bw.ucv(x), 
                 bcv = bw.bcv(x), 
                 sj = , `sj-ste` = bw.SJ(x, method = "ste"), 
                 `sj-dpi` = bw.SJ(x, method = "dpi"), 
                 stop("unknown bandwidth rule"))
  }
  if (!is.finite(bw)) 
    stop("non-finite 'bw'")
  bw <- adjust * bw
  if (bw <= 0) 
    stop("'bw' is not positive.")
  if (missing(from)) 
    from <- min(x) - cut * bw
  if (missing(to)) 
    to <- max(x) + cut * bw
  if (!is.finite(from)) 
    stop("non-finite 'from'")
  if (!is.finite(to)) 
    stop("non-finite 'to'")
  lo <- from - 4 * bw
  up <- to + 4 * bw
  y <- .Call(stats:::C_BinDist, x, weights, lo, up, n) * totMass
  kords <- seq.int(0, 2 * (up - lo), length.out = 2L * n)
  kords[(n + 2):(2 * n)] <- -kords[n:2]
  
  dhuber <- function(x, bw, k = 1.5)
  {
    z = x - mean(x) / sd(x)
    cnorm <- (2*pnorm(k)-1) + 2*dnorm(k)/k
    ifelse(abs(z) < k, dnorm(z, sd = bw), 
            dlaplace(z, scale = 2 * sqrt(bw)))/cnorm
  }
  
  dlaplace = function (x, location = 0, scale = 1, log = FALSE) 
  {
    x <- as.vector(x)
    location <- as.vector(location)
    scale <- as.vector(scale)
    if (any(scale <= 0)) 
      stop("The scale parameter must be positive.")
    NN <- max(length(x), length(location), length(scale))
    x <- rep(x, len = NN)
    location <- rep(location, len = NN)
    scale <- rep(scale, len = NN)
    dens <- (-abs(x - location)/scale) - log(2 * scale)
    if (log == FALSE) 
      dens <- exp(dens)
    return(dens)
  }
  
  dstudent = function(x, df = 6, mu = 0, sigma = 1){
    dt((x - mu)/sigma, df = df)/sigma
  }
  
  kords <- switch(kernel, 
                  normal = dnorm(kords, sd = bw), 
                  laplacian = dlaplace(kords, scale = bw), 
                  laplace = dlaplace(kords, scale = bw), 
                  huber = dhuber(kords, bw = bw, k = 1.58),
                  student = dstudent(kords, df = df, mu = 0, sigma = bw), 
                  rectangular = {
                    a <- bw * sqrt(3)
                    ifelse(abs(kords) < a, 0.5/a, 0)
                  }, triangular = {
                    a <- bw * sqrt(6)
                    ax <- abs(kords)
                    ifelse(ax < a, (1 - ax/a)/a, 0)
                  }, epanechnikov = {
                    a <- bw * sqrt(5)
                    ax <- abs(kords)
                    ifelse(ax < a, 3/4 * (1 - (ax/a)^2)/a, 0)
                  }, biweight = {
                    a <- bw * sqrt(7)
                    ax <- abs(kords)
                    ifelse(ax < a, 15/16 * (1 - (ax/a)^2)^2/a, 0)
                  }, cosine = {
                    a <- bw/sqrt(1/3 - 2/pi^2)
                    ifelse(abs(kords) < a, (1 + cos(pi * kords/a))/(2 * 
                                                                      a), 0)
                  }, optcosine = {
                    a <- bw/sqrt(1 - 8/pi^2)
                    ifelse(abs(kords) < a, pi/4 * cos(pi * kords/(2 * 
                                                                    a))/a, 0)
                  })
  kords <- fft(fft(y) * Conj(fft(kords)), inverse = TRUE)
  kords <- pmax.int(0, Re(kords)[1L:n]/length(y))
  xords <- seq.int(lo, up, length.out = n)
  x <- seq.int(from, to, length.out = n.user)
  structure(list(x = x, y = approx(xords, kords, x)$y, bw = bw, 
                 n = N, call = match.call(), data.name = name, has.na = FALSE), 
            class = "density")
}