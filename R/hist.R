#' Modified Version of hist() with additional functionality
#'
#' @param x a vector of values for which the histogram is desired.
#' @param breaks 	one of: a vector giving the breakpoints between histogram cells, a function to compute the vector of breakpoints,
#' a single number giving the number of cells for the histogram, or a character string naming an algorithm to compute the number of cells (see ‘Details’), a function to compute the number of cells. Defaults to "dhist".
#' @param freq 	logical; if TRUE, the histogram graphic is a representation of frequencies, the counts component of the result; if FALSE, probability densities, component density, are plotted (so that the histogram has a total area of one). Defaults to FALSE.
#' @param probability 	an alias for !freq, for S compatibility.
#' @param include.lowest logical; if TRUE, an x[i] equal to the breaks value will be included in the first (or last, for right = FALSE) bar. This will be ignored (with a warning) unless breaks is a vector.
#' @param right logical; if TRUE, the histogram cells are right-closed (left open) intervals.
#' @param density the density of shading lines, in lines per inch. The default value of NULL means that no shading lines are drawn. Non-positive values of density also inhibit the drawing of shading lines.
#' @param angle the slope of shading lines, given as an angle in degrees (counter-clockwise).
#' @param col a colour to be used to fill the bars. The default is "#71D795CC".
#' @param border the color of the border around the bars. The default is to use the standard foreground color.
#' @param main 	these arguments to title have useful defaults here.
#' @param xlim the range of x and y values with sensible defaults. Note that xlim is not used to define the histogram (breaks), but only for plotting (when plot = TRUE).
#' @param ylim the range of x and y values with sensible defaults. Note that xlim is not used to define the histogram (breaks), but only for plotting (when plot = TRUE).
#' @param xlab these arguments to title have useful defaults here.
#' @param ylab these arguments to title have useful defaults here.
#' @param axes logical. If TRUE (default), axes are draw if the plot is drawn.
#' @param rug Should a rug be plotted under the histogram? Defaults to TRUE. 
#' @param rug.col The rug color. 
#' @param plot 	logical. If TRUE (default), a histogram is plotted, otherwise a list of breaks and counts is returned. In the latter case, a warning is used if (typically graphical) arguments are specified that only apply to the plot = TRUE case.
#' @param labels logical or character string. Additionally draw labels on top of bars, if not FALSE; see plot.histogram.
#' @param nclass 	numeric (integer). For S(-PLUS) compatibility only, nclass is equivalent to breaks for a scalar or character argument.
#' @param warn.unused 	logical. If plot = FALSE and warn.unused = TRUE, a warning will be issued when graphical parameters are passed to hist.default().
#' @param ... further arguments and graphical parameters passed to plot.histogram and thence to title and axis (if plot = TRUE).
#'
#' @return a histogram
#' @export
#'
#' @examples
#' x <- rgamma(1000, 2, .25)
#' par(mfrow = c(2, 2))
#' hist(x, breaks = "dhist")
#' hist(x, breaks = "fd")
#' par(mfrow = c(2, 2))
#' hist(x, breaks = "dhist", main = "'dhist' breaks method")
#' hist(x, breaks = "fd", main = "'Freedman-Diaconis' breaks method")
#' hist(x, breaks = "scott", main = "'Scott' breaks method")
#' hist(x, breaks = "sturges", main = "'Sturges' breaks method")
#' 
#' 
#' @details 
#' The definition of histogram differs by source (with country-specific biases). 
#' R's default with equi-spaced breaks (also the default) is to plot the counts in the cells defined by breaks. 
#' Thus the height of a rectangle is proportional to the number of points falling into the cell, as is the area 
#' provided the breaks are equally-spaced. The modified hist() function included in the Bayezilla package however
#' utilizes non-equal sized breaks. More details are given below.
#' If right = TRUE (default), the histogram cells are intervals of the form (a, b], i.e., 
#' they include their right-hand endpoint,
#' but not their left one, with the exception of the first cell when include.lowest is TRUE.
#' For right = FALSE, the intervals are of the form [a, b), and include.lowest means ‘include highest’.
#' A numerical tolerance of 1e-7 times the median bin size (for more than four bins, otherwise the median is substituted) 
#' is applied when counting entries on the edges of bins. This is not included in the reported breaks nor in the calculation
#' of density.  \cr
#' \cr
#' Breaks Algorithms: \cr
#' \cr
#' The default breaks algorithm is "dhist", which implements the varying binwidth algorithm of Lorraine Denby. 
#' This algorithm has some notable
#' advantages from a statistical point of view. Regions of high density have not only taller bins (as is usual) 
#' but more narrow bins as well. Regions of lower denisty have not only shorter, but wider bins. This makes the 
#' probability density much more 
#' immediately obvious, and captures interesting features of heavy tails and skew with greater efficacy. Algorithms
#' that yield The former oversmooths in regions of high density, and is poor at identifying sharp peaks and multimodality. 
#' By contrast, the latter variety oversmooths in regions of low density and can mask outliers and the 
#' heavy tails of more leptokurtotic distributions. For more information, see Denby & Mallows (2009).  
#' Other options include "scott", "sturges" and "fd" / "Freedman-Diaconis". Case is ignored and partial matching is used.  
#' Alternatively, a function can be supplied which will compute the
#' intended number of breaks or the actual breakpoints as a function of x. 
#' \cr
#' An example of output: \cr
#' \cr
#' \if{html}{\figure{breaks.png}{}}
#' \if{latex}{\figure{breaks.png}{}}
#' \cr
#' @references Denby, L., & Mallows, C. (2009). Variations on the Histogram. Journal of Computational and Graphical Statistics, 18(1), 21–31. doi:10.1198/jcgs.2009.0002 
#' 

hist <- function(x, breaks = "dhist", freq = FALSE, probability = !freq, 
                 include.lowest = TRUE, right = TRUE, density = NULL, angle = 45, 
                 col = "#71D795CC", border = NULL, main = NULL, 
                 xlim = range(breaks), ylim = NULL, xlab = xname, ylab, axes = TRUE, 
                 plot = TRUE, labels = FALSE, nclass = NULL, warn.unused = TRUE, rug = TRUE, 
                 rug.col = "#445edd",
                 ...) { UseMethod("hist") }

#' Modified Version of hist() with additional functionality
#'
#' @param x a vector of values for which the histogram is desired.
#' @param breaks 	one of: a vector giving the breakpoints between histogram cells, a function to compute the vector of breakpoints,
#' a single number giving the number of cells for the histogram, or a character string naming an algorithm to compute the number of cells (see ‘Details’), a function to compute the number of cells. Defaults to "dhist".
#' @param freq 	logical; if TRUE, the histogram graphic is a representation of frequencies, the counts component of the result; if FALSE, probability densities, component density, are plotted (so that the histogram has a total area of one). Defaults to FALSE.
#' @param probability an alias for !freq, for S compatibility.
#' @param include.lowest logical; if TRUE, an x[i] equal to the breaks value will be included in the first (or last, for right = FALSE) bar. This will be ignored (with a warning) unless breaks is a vector.
#' @param right logical; if TRUE, the histogram cells are right-closed (left open) intervals.
#' @param density the density of shading lines, in lines per inch. The default value of NULL means that no shading lines are drawn. Non-positive values of density also inhibit the drawing of shading lines.
#' @param angle the slope of shading lines, given as an angle in degrees (counter-clockwise).
#' @param col a colour to be used to fill the bars. The default is "#71D795CC".
#' @param border the color of the border around the bars. The default is to use the standard foreground color.
#' @param main 	these arguments to title have useful defaults here.
#' @param xlim the range of x and y values with sensible defaults. Note that xlim is not used to define the histogram (breaks), but only for plotting (when plot = TRUE).
#' @param ylim the range of x and y values with sensible defaults. Note that xlim is not used to define the histogram (breaks), but only for plotting (when plot = TRUE).
#' @param xlab these arguments to title have useful defaults here.
#' @param ylab these arguments to title have useful defaults here.
#' @param axes logical. If TRUE (default), axes are draw if the plot is drawn.
#' @param rug Should a rug be plotted under the histogram? Defaults to TRUE. 
#' @param rug.col The rug color. 
#' @param plot 	logical. If TRUE (default), a histogram is plotted, otherwise a list of breaks and counts is returned. In the latter case, a warning is used if (typically graphical) arguments are specified that only apply to the plot = TRUE case.
#' @param labels logical or character string. Additionally draw labels on top of bars, if not FALSE; see plot.histogram.
#' @param nclass 	numeric (integer). For S(-PLUS) compatibility only, nclass is equivalent to breaks for a scalar or character argument.
#' @param warn.unused 	logical. If plot = FALSE and warn.unused = TRUE, a warning will be issued when graphical parameters are passed to hist.default().
#' @param ... further arguments and graphical parameters passed to plot.histogram and thence to title and axis (if plot = TRUE).
#'
#' @return a histogram
#'
#' @examples
#' x <- rgamma(1000, 2, .25)
#' par(mfrow = c(2, 2))
#' hist(x, breaks = "dhist")
#' hist(x, breaks = "fd")
#' par(mfrow = c(2, 2))
#' hist(x, breaks = "dhist", main = "'dhist' breaks method")
#' hist(x, breaks = "fd", main = "'Freedman-Diaconis' breaks method")
#' hist(x, breaks = "scott", main = "'Scott' breaks method")
#' hist(x, breaks = "sturges", main = "'Sturges' breaks method")
#' 
#' @details 
#' The definition of histogram differs by source (with country-specific biases). 
#' R's default with equi-spaced breaks (also the default) is to plot the counts in the cells defined by breaks. 
#' Thus the height of a rectangle is proportional to the number of points falling into the cell, as is the area 
#' provided the breaks are equally-spaced. The modified hist() function included in the Bayezilla package however
#' utilizes non-equal sized breaks. More details are given below.
#' If right = TRUE (default), the histogram cells are intervals of the form (a, b], i.e., 
#' they include their right-hand endpoint,
#' but not their left one, with the exception of the first cell when include.lowest is TRUE.
#' For right = FALSE, the intervals are of the form [a, b), and include.lowest means ‘include highest’.
#' A numerical tolerance of 1e-7 times the median bin size (for more than four bins, otherwise the median is substituted) 
#' is applied when counting entries on the edges of bins. This is not included in the reported breaks nor in the calculation
#' of density.  \cr
#' \cr
#' Breaks Algorithms: \cr
#' \cr
#' The default breaks algorithm is "dhist", which implements the varying binwidth algorithm of Lorraine Denby. 
#' This algorithm has some notable
#' advantages from a statistical point of view. Regions of high density have not only taller bins (as is usual) 
#' but more narrow bins as well. Regions of lower denisty have not only shorter, but wider bins. This makes the 
#' probability density much more 
#' immediately obvious, and captures interesting features of heavy tails and skew with greater efficacy. Algorithms
#' that yield The former oversmooths in regions of high density, and is poor at identifying sharp peaks and multimodality. 
#' By contrast, the latter variety oversmooths in regions of low density and can mask outliers and the 
#' heavy tails of more leptokurtotic distributions. For more information, see Denby & Mallows (2009).  
#' Other options include "scott", "sturges" and "fd" / "Freedman-Diaconis". Case is ignored and partial matching is used.  
#' Alternatively, a function can be supplied which will compute the
#' intended number of breaks or the actual breakpoints as a function of x. 
#' \cr
#' An example of output: \cr
#' \cr
#' \if{html}{\figure{breaks.png}{}}
#' \if{latex}{\figure{breaks.png}{}}
#' \cr
#' \cr
#' @references Denby, L., & Mallows, C. (2009). Variations on the Histogram. Journal of Computational and Graphical Statistics, 18(1), 21–31. doi:10.1198/jcgs.2009.0002 
#' @export 
#' 
hist.default <- function (x, breaks = "dhist", freq = FALSE, probability = !freq, 
                          include.lowest = TRUE, right = TRUE, density = NULL, angle = 45, 
                          col = "#71D795CC", border = NULL, main = NULL, 
                          xlim = range(breaks), ylim = NULL, xlab = xname, ylab, axes = TRUE, 
                          plot = TRUE, labels = FALSE, nclass = NULL, warn.unused = TRUE, rug = TRUE, 
                          rug.col = "#445edd",
                          ...) 
{
  if (!is.numeric(x)) 
    stop("'x' must be numeric")
  xname <- paste(deparse(substitute(x), 500), collapse = "\n")
  n <- length(x <- x[is.finite(x)])
  n <- as.integer(n)
  if (is.na(n)) 
    stop("invalid length(x)")
  use.br <- !missing(breaks)
  if (use.br) {
    if (!missing(nclass)) 
      warning("'nclass' not used when 'breaks' is specified")
  }
  else if (!is.null(nclass) && length(nclass) == 1L) 
    breaks <- nclass
  use.br <- use.br && (nB <- length(breaks)) > 1L
  if (use.br) 
    breaks <- sort(breaks)
  else {
    if (!include.lowest) {
      include.lowest <- TRUE
      warning("'include.lowest' ignored as 'breaks' is not a vector")
    }
    if (is.character(breaks)) {
      breaks <- match.arg(tolower(breaks), c("sturges", "dhist",
                                             "fd", "freedman-diaconis", "scott"))
      breaks <- switch(breaks, sturges = nclass.Sturges(x), 
                       `freedman-diaconis` = , dhist = nclass.dhist(x), fd = nclass.FD(x), scott = nclass.scott(x), 
                       stop("unknown 'breaks' algorithm"))
    }
    else if (is.function(breaks)) {
      breaks <- breaks(x)
    }
    if (length(breaks) == 1) {
      if (!is.numeric(breaks) || !is.finite(breaks) || 
          breaks < 1L) 
        stop("invalid number of 'breaks'")
      if (breaks > 1e+06) {
        warning(gettextf("'breaks = %g' is too large and set to 1e6", 
                         breaks), domain = NA)
        breaks <- 1000000L
      }
      breaks <- pretty(range(x), n = breaks, min.n = 1)
      nB <- length(breaks)
      if (nB <= 1) 
        stop(gettextf("hist.default: pretty() error, breaks=%s", 
                      format(breaks)), domain = NA)
    }
    else {
      if (!is.numeric(breaks) || length(breaks) <= 1) 
        stop(gettextf("Invalid breakpoints produced by 'breaks(x)': %s", 
                      format(breaks)), domain = NA)
      breaks <- sort(breaks)
      nB <- length(breaks)
      use.br <- TRUE
    }
  }
  nB <- as.integer(nB)
  if (is.na(nB)) 
    stop("invalid length(breaks)")
  h <- diff(breaks)
  equidist <- !use.br || diff(range(h)) < 1e-07 * mean(h)
  if (!use.br && any(h <= 0)) 
    stop("'breaks' are not strictly increasing")
  freq1 <- freq
  if (is.null(freq)) {
    freq1 <- if (!missing(probability)) 
      !as.logical(probability)
    else equidist
  }
  else if (!missing(probability) && any(probability == freq)) 
    stop("'probability' is an alias for '!freq', however they differ.")
  diddle <- 1e-07 * if (nB > 5) 
    stats::median(h)
  else if (nB <= 3) 
    diff(range(x))
  else min(h[h > 0])
  fuzz <- if (right) 
    c(if (include.lowest) -diddle else diddle, rep.int(diddle, 
                                                       nB - 1L))
  else c(rep.int(-diddle, nB - 1L), if (include.lowest) diddle else -diddle)
  fuzzybreaks <- breaks + fuzz
  counts <- .Call(graphics:::C_BinCount, x, fuzzybreaks, right, include.lowest)
  if (any(counts < 0L)) 
    stop("negative 'counts'. Internal Error.", domain = NA)
  if (sum(counts) < n) 
    stop("some 'x' not counted; maybe 'breaks' do not span range of 'x'")
  dens <- counts/(n * h)
  mids <- 0.5 * (breaks[-1L] + breaks[-nB])
  r <- structure(list(breaks = breaks, counts = counts, density = dens, 
                      mids = mids, xname = xname, equidist = equidist), class = "histogram")
  if (plot) {
    plot(r, freq = freq1, col = col, border = border, angle = angle, 
         density = density, main = main, xlim = xlim, ylim = ylim, 
         xlab = xlab, ylab = ylab, axes = axes, labels = labels, 
         ...)

    invisible(r)
    if (isTRUE(rug)){
      rug(x, ticksize = 0.01, side = 1, lwd = 0.5, col = rug.col,
          quiet = getOption("warn") < 0, ...)
    }
  }
  else {
    if (warn.unused) {
      nf <- names(formals())
      nf <- nf[is.na(match(nf, c("x", "breaks", "nclass", 
                                 "plot", "include.lowest", "right")))]
      missE <- lapply(nf, function(n) substitute(missing(.), 
                                                 list(. = as.name(n))))
      not.miss <- !sapply(missE, eval, envir = environment())
      if (any(not.miss)) 
        warning(sprintf(ngettext(sum(not.miss), "argument %s is not made use of", 
                                 "arguments %s are not made use of"), paste(sQuote(nf[not.miss]), 
                                                                            collapse = ", ")), domain = NA)
    }
    r
  }
}

#' Density-based adaptive width histogram bins
#'
#' @description This algorithm has some notable
#' advantages from a statistical point of view. Regions of high density have not only taller bins (as is usual) 
#' but more narrow bins as well. Regions of lower denisty have not only shorter, but wider bins. This makes the 
#' probability density much more 
#' immediately obvious, and captures interesting features of heavy tails and skew with greater efficacy. Algorithms
#' that yield The former oversmooths in regions of high density, and is poor at identifying sharp peaks and multimodality. 
#' By contrast, the latter variety oversmooths in regions of low density and can mask outliers and the 
#' heavy tails of more leptokurtotic distributions. Alternatively, a function can be supplied which will compute the
#' intended number of breaks or the actual breakpoints as a function of x. For more information, see Denby & Mallows (2009). 
#'
#' @author Lorraine Denby
#' @seealso
#'   \url{http://pubs.research.avayalabs.com/pdfs/ALR-2007-003-paper.pdf}
#' @references L. Denby and C. Mallows. Variations on the histogram. Journal
#'   of Computational and Graphical Statistics, 18 (1):21-31, 2009.
#'   URL \url{http://pubs.amstat.org/doi/abs/10.1198/jcgs.2009.0002.}
#' @param b slope. See paper for details. Defaults to 1.5.
#' @param nbins number of bins. See Denby & Mallows (2009).
#' @param min.bins the minimum number of bins. 
#' @param rx range of data, if not taken from data.
#' @return A function that takes a single parameter, a numeric x specifying
#'   the data for which breaks are needed, and returns a vector of breaks.
#' @export
#' @examples 
#' nclass.dhist()
#' 
nclass.dhist = function(x, b = 1.5, rx = range(x), min.bins = 9) {
  
  n.bins <-function(x, min.bins = 9){
    x<-x[!is.na(x)]
    n<-length(x)
    Q<-quantile(x, c(0.025 , 0.98))
    X<-range(x)
    result<- ceiling(n^(1/3)*(X[2]-X[1])/(2*(Q[2]-Q[1])))
    result <- max(c(min.bins, result))
    names(result)<-NULL
    result
  }
  
  a = b*diff(quantile(x, c(0.025 , 0.98)))
  
  nbins <- n.bins(x, min.bins = min.bins)
  
  x <- sort(x)
  if(a == 0)
    a <- diff(range(x))/100000000
  
  if(a != Inf) {
    n <- length(x)
    h <- (rx[2] + a - rx[1])/nbins
    ybr <- rx[1] + h * (0:nbins)
    yupper <- x + (a * (1:n))/n
    
    # upper and low Fer corners in the ecdf
    ylower <- yupper - a/n
    
    cmtx <- cbind(
      cut(yupper, breaks = ybr), 
      cut(yupper, breaks = ybr, left.include = T), 
      cut(ylower, breaks = ybr),
      cut(ylower, breaks = ybr, left.include = T)
    )
    cmtx[1, 3] <- cmtx[1, 4] <- 1
    # to replace NAs when default r is used
    cmtx[n, 1] <- cmtx[n, 2] <- nbins
    
    checksum <- rowSums(cmtx) %% 4
    # will be 2 for obs. that straddle two bins
    straddlers <- (1:n)[checksum == 2]
    
    # to allow for zero counts
    if(length(straddlers) > 0) {
      counts <- table(c(1:nbins, cmtx[- straddlers, 1])) 
    } else {
      counts <- table(c(1:nbins, cmtx[, 1]))
    }
    counts <- counts - 1
    
    if(length(straddlers) > 0) {
      for(i in straddlers) {
        binno <- cmtx[i, 1]
        theta <- ((yupper[i] - ybr[binno]) * n)/a
        counts[binno - 1] <- counts[binno - 1] + (1 - theta)
        counts[binno] <- counts[binno] + theta
      }
    }
    xbr <- ybr
    xbr[-1] <- ybr[-1] - (a * cumsum(counts))/n
  } else {
    bin.size <- length(x)/nbins
    cut.pt <- c(
      min(x) - abs(min(x))/1000, 
      approx(seq(length(x)),
             x, 
             (1:(nbins - 1)) * bin.size, rule = 2)$y, 
      max(x)
    )
    aa <- hist(x, breaks = cut.pt, plot = F, probability = T)
    xbr <- aa$breaks
  }
  
  xbr
}
