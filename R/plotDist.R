#' plot a density function over a histogram
#'
#' @param x the data
#' @param dfunc The name of a density function (ie, dnorm, dt, etc)
#' @param args list of the additional argument of the function, ie, list(mu = 0, sigma = 1)
#' @param bins defaults to 20
#' @param n the number of x values to be evaluated in the superimposed density curve. Defaults to 1500.
#' @param type defaults to "irregular" using the binning method described in plotPost. Otherwise, "equal".
#' @param xname character; argument of the function in func, the name of the x axis
#' @param xlab label of the x-axis
#' @param ylab label of the y-axis
#' @param xlim restrict the range of the function
#' @param color the color of the line. defaults to "red"
#' @param size the size of the line. defaults to 2
#' @param ... Additional arguments to stat_function()
#' @examples
#' plotDist(dnorm, from = -8 to = 8)
#' @export
#' 
#' 
plotDist = function (x = NULL, dfunc, args = list(), bins = 20, n = 1500, type = "irregular", xname = "\u03C7", 
            xlab = xname, ylab = NULL, xlim = NULL,  
            color = "#f45342CC", 
            size = 2, 
            ...) 
  {
    
    sexpr <- substitute(dfunc)
    if (is.function(dfunc)) {
      funky_function <- dfunc
    }
    
    else {
      stop("please provide a function")
    }
    
    if (is.null(ylab)) 
      ylab <- paste0("ƒ ", "( ", "\u03C7", " )\n")
    
    if (is.list(args) == FALSE){
      args = as.list(args)
    }
    
    dhist <- function(x, a=5*diff(quantile(x, c(0.05,.95))), nbins = bins, rx = range(x)) {
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
    
    
    if (type == "equal") {
    tibble(x = x) %>% 
      ggplot(aes(x = x)) + 
      geom_histogram(bins = bins, aes(y = stat(density))) + 
      stat_function(
        fun = funky_function, 
                    n = n, 
                    xlim = xlim, 
                    args = args, 
                    color = color, 
                    size = size, 
                    ...) + 
      labs(x = xlab, y = ylab)
    } else{
      brks = dhist(x, nbins = bins)
      tibble(x = x) %>% 
        ggplot(aes(x = x)) + 
        geom_histogram(breaks = brks, aes(y = stat(density))) + 
        stat_function(
          fun = funky_function, 
          n = n, 
          xlim = xlim, 
          args = args, 
          color = color, 
          size = size, 
          ...) + 
        labs(x = xlab, y = ylab)
    }
      
      
  }