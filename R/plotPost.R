#' Histogram of Posterior Distribution
#'
#' @description This is a heavy modification of John Kruchke's function of the same name in the BEST
#' package. This however makes many aesthetic changes and has some different features from the
#' original function.
#' \cr
#' \cr
#' A unique feature of this function is that by default irregular sized bins are used. Why? Although
#' this may appear unsightly to you at first, it actually provides more information. There are two general
#' classes of histogram binning: equal-width and equal-area histograms. The former oversmooths in
#' regions of high density, and is poor at identifying sharp peaks and multimodality. 
#' By contrast, the latter variety
#' oversmooths in regions of low density and can mask outliers and the heavy tails of more
#' leptokurtotic distributions (such as the Student-t and Laplacian distributions). The irregular
#' binned histogram on the other hand lacks these faults, and is better at showing structure of 
#' the plotted distribution. For plotting posterior distributions I find that this can aid interpretation
#' of the posterior. Furthermore, regions of high density have not only taller bins (as is usual) but 
#' more narrow bins as well. Regions of lower denisty have not only shorter, but wider bins. This
#' makes the probability density much more immediately obvious, and captures interesting features of
#' heavy tails and skew with greater efficacy. While you can turn off this feature with type = "equal",
#' I urge you to give the irregular binning method a chance before dismissing it as unaesthetically
#' appealing. For more information, see Denby & Mallows (2009). \cr
#' \cr
#' See the details section at the end for an example of the output of this function.
#' 
#' @details 
#' An example of output: \cr
#' \cr
#' \if{html}{\figure{plotPost.png}{}}
#' \if{latex}{\figure{plotPost.png}{}}
#'
#' 
#' @references Denby, L., & Mallows, C. (2009). Variations on the Histogram. Journal of Computational and Graphical Statistics, 18(1), 21–31. doi:10.1198/jcgs.2009.0002 \cr
#' 
#'
#' @param paramSampleVec a vector containing the posterior distribution
#' @param fit a stanfit or runjags object. This can be used as an alternative to paramSampleVec,
#' but you must specify which parameter you would like to plot. If using this argument be sure to type
#' fit = "yourmodel" so that the function knows it is not intended to be a vector.
#' @param param the name of the parameter in the stanfit or runjags object you want to plot.
#' @param col the color scheme. One of "blue", "green" "red", or "purple".
#' @param xlab if desired, a custom x-axis label. If supplying a fit object the label will be
#' the parameter name, but you can override that. If paramSampleVec is supplied without a custom
#' label here, the x-axis will display the Greek theta symbol by default.
#' @param cred.level The credibility level. Defaults to 90\% (.90).
#' @param method Quantile Intervals "QI" (the default) or highest density intervals "HDI"
#' @param CItextPlace Adjust the position of the credible interval text if needed.
#' @param bins Adjust the number of bins in the histogram. Default is 20.
#' @param type Either "irregular" (default) which uses uneven sized bins proportional to probability density or "equal" for even sized bins. 
#' @param showMedian Should the median be used instead of the mean?
#' @param ROPE If you would like to display a ROPE, enter vector of three numbers, the first being the
#' comparison value, the second being the lower limit of the ROPE and the third being the upper limit of the ROPE.
#' @param ... other arguments to hist()
#' @return
#' a histogram
#' @export
#' @examples
#' plotPost()
plotPost <- function(paramSampleVec, fit = NULL, param = NULL, xlab = NULL, col = "blue", ROPE = NULL, cred.level = 0.90, method = "QI", CItextPlace = 0.7, type = "irregular", bins = 20, showMedian = FALSE, ...) {
  
  # dhist.
  # An another algorithm for computing histogram breaks.  Produces irregular bins.
  # Provided by Lorraine Denby
  
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
  
  old.par <- par(no.readonly = TRUE) # save default, for resetting... 
  on.exit(par(old.par))     #and when we quit the function, restore to original values
  
  
  if (is.null(fit) != TRUE){
    if (is.null(param) == TRUE)  {
      stop("please choose a single parameter to plot with the 'param' argument")
      } 
    stan <- inherits(fit, "stanfit")
    if (stan == TRUE) {
      paramSampleVec <- as.matrix(fit)
      paramSampleVec = as.vector(paramSampleVec[,which(colnames(paramSampleVec) == param)])
    } 
    else if (class(fit) == "runjags") {
      paramSampleVec = as.matrix(combine.mcmc(fit, collapse.chains = TRUE, vars = param))
    }
    
    param.label <- noquote(param)
  }
  else {
    if (is.null(param)){
      param.label <- expression(theta)
    }
    else if (!is.null(param)){
      param.label <- noquote(param)
    }
  }
  
  if (is.null(xlab) != TRUE){
    param.label <-  noquote(xlab)
  }
  
  #light, dark, dark, light
  if (col == "blues" || col == "blue"){
    ColorScheme = c("#6495edAD", "#0d2f6dCC", "#0a2556", "#046eff")
  }
  if (col == "purples" || col == "purple"){
    ColorScheme = c("#d4aad4CC" , "#400040CC", "#270027", "#a700a7")
  }
  if (col == "reds" || col == "red"){
    ColorScheme = c("#e87a7aCC", "#8b1a1aCC", "#910202", "#e80039")
  }
  if (col == "greens" || col == "green"){
    ColorScheme= c("#66cc66CC", "#194d33CC", "#133a26", "#39a636")
  }
  
  # Get point estimate information for the title
  
  if ( showMedian==FALSE ) {
    Param <- round(mean( paramSampleVec ), 2)
    estimate.text = paste("Mean", ":", Param)
  } else {
    Param <- round(median( paramSampleVec ), 2)
    estimate.text = paste("Median", ":", Param)
  }
  
  # Get HDI information for the title
  if (method == "QI" || method == "ETI") {
    HDI <- round(cred_interval( paramSampleVec, cred.level = cred.level, method = "QI"), 2)
    interval.text = paste(100*cred.level, "% QI", ":", "[", HDI[1], "," ,  HDI[2], "]")
  }
  else if(method == "HDI" || method == "HDP") {
    HDI <- round(cred_interval( paramSampleVec, cred.level = cred.level, method = "HDI"), 2)
    interval.text = paste(100*cred.level, "% HDI", ":", "[", HDI[1], "," , HDI[2], "]")
  }
  
  main.title = paste(estimate.text, "\n", interval.text)
  
  # Deal with ... argument:
  dots <- list(...)
  if(length(dots) == 1 && class(dots[[1]]) == "list")
    dots <- dots[[1]]
  defaultArgs <- list(xlab= param.label,
                      yaxt="n", ylab="", main= main.title, cex.lab=1.5,
                      cex=1.4, col=ColorScheme[1], border=ColorScheme[2], bty="n", xaxt = "n",
                      family = 'serif',
                      lwd=6, freq=FALSE,
                      xlim = cred_interval(paramSampleVec, cred.level=0.999999, method = method))
  useArgs <- modifyList(defaultArgs, dots)
  
  # Get breaks argument
  if (type == "irregular"){
    breaks <- dhist
  } else {
    breaks <- bins
  }
  
  histinfo <- hist(paramSampleVec, breaks=breaks, plot=FALSE)
  histinfo$xname <- useArgs$xlab
  
  
  oldpar <- par(xpd=TRUE) ; on.exit(par(oldpar))
  
  
  plot.histogram.args.names <- c("freq", "density", "angle", "border",
                                 "main", "sub", "xlab", "ylab", "xlim", "ylim", "axes", "labels",
                                 "add") # plot.histogram not exported, so need to cheat!
  selPlot <- names(useArgs) %in%
    c(plot.histogram.args.names, names(par(no.readonly=TRUE)))
  plotArgs <- useArgs[selPlot]
  plotArgs$lwd <- 1
  plotArgs$x <- histinfo
  do.call(plot, plotArgs, quote=TRUE)
  
  # Display the HDI.
  
  if (method == "QI" || method == "ETI") {
    HDI <- cred_interval( paramSampleVec, cred.level = cred.level, method = "QI")
    SE <-  cred_interval( paramSampleVec, cred.level = .682, method = "QI")
    lines(HDI, c(0,0), lwd=4.5, lend='butt', col = ColorScheme[3])
    lines(SE, c(0,0), lwd=5.45, lend='butt', col = ColorScheme[4])
  }  else if(method == "HDI" || method == "HDP") {
    HDI <- cred_interval( paramSampleVec, cred.level = cred.level, method = "HDI")
    SE <-  cred_interval( paramSampleVec, cred.level = .682, method = "HDI")
    lines(HDI, c(0,0), lwd=4.5, lend='butt', col = ColorScheme[3])
    lines(SE, c(0,0), lwd=5.45, lend='butt', col = ColorScheme[4])
  }
  
  
  if (is.null(ROPE) == FALSE) {
    compVal = ROPE[1]
    ROPE = ROPE[-1]
    ROPEtextHt = 0.55*max(histinfo$density)
    # Display the ROPE.
    lines( c(compVal,compVal) , c(0.96*ROPEtextHt,0) , lty="dashed" , lwd= 2,
           col= "#9bacc9")
    ropeCol = "#19386b"
    pcInROPE = ( sum( paramSampleVec > ROPE[1] & paramSampleVec < ROPE[2] )
                 / length( paramSampleVec ) )
    lines( c(ROPE[1],ROPE[1]) , c(0.96*ROPEtextHt,0) , lty="dotted" , lwd=3 ,
           col=ropeCol )
    lines( c(ROPE[2],ROPE[2]) , c(0.96*ROPEtextHt,0) , lty="dotted" , lwd=3 ,
           col=ropeCol)
    text( mean(ROPE) , ROPEtextHt ,
          bquote( .(round(100*pcInROPE))*"% in ROPE" ) ,
          adj=c(.5,0) , cex=1 , col=ropeCol )
  }
  
  points(x = Param, y = 0, pch = 18, cex = 1.75)
  axis(1, col = NA, tck = 0, family = 'serif')
  return(invisible(histinfo))
}