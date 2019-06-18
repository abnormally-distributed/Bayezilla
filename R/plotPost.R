#' Histogram of Posterior Distribution
#'
#' @description This is a heavy modification of John Kruchke's function of the same name in the BEST
#' package. This however makes many aesthetic changes and has some different features from the
#' original function.
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
#' @param breaks Adjust the number of bins in the histogram
#' @param showMedian Should the median be used instead of the mean?
#' @param ROPE If you would like to display a ROPE, enter vector of three numbers, the first being the
#' comparison value, the second being the lower limit of the ROPE and the third being the upper limit of the ROPE.
#' @param ... other arguments to hist()
#' @return
#' a histogram
#' @export
#' @examples
#' plotPost()
plotPost <- function(paramSampleVec, fit = NULL, param = NULL, xlab = NULL, col = "blue", ROPE = NULL, cred.level = 0.90, method = "QI", CItextPlace = 0.7, breaks = 40, showMedian = FALSE, ...) {

  
  old.par <- par(no.readonly = TRUE) # save default, for resetting... 
  on.exit(par(old.par))     #and when we quit the function, restore to original values
  if (is.null(fit) != TRUE){
    if (is.null(param) == TRUE)
      stop("please choose a single parameter to plot with the 'param' argument")
    paramSampleVec = as.matrix(combine.mcmc(fit,collapse.chains = TRUE, vars = param))
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
    breaks <- breaks
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
