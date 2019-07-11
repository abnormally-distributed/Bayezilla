#' Kernel Density Plot
#'
#' @param x the numeric vector of interest
#' @param rug should a rug of the data be plotted underneath the density? Defaults to TRUE.
#' @param col The fill color of the density curve. Defaults to "#45a3ffCC"
#' @param border The border color of the density curve. Defaults to "#000f1e"
#' @param rug.col The color of the data rug. Defaults to "#62a3e2".
#' @param main.title Defaults to "Density Plot".
#' @param ... Arguments to pass to density()
#'
#' @return a plot
#' @export 
#'
#' @examples densPlot(rnorm(1000))
densPlot <- function(x, rug = TRUE, col = "#45a3ffCC", border = "#000f1e", rug.col = "#62a3e2", main.title = "Density Plot", ...){
  d <- density(x, ...)
  plot(d, yaxt="n", ylab="", main= main.title, cex.lab=1.3, cex=1.8, bty="n", family = 'serif', lwd=1)
  polygon(d, col=col, border = border)
  if (isTRUE(rug)){
    rug(x, ticksize = 0.01, side = 1, lwd = 0.5, col = rug.col,
        quiet = getOption("warn") < 0)
  }
}
