#' Plot a large number of variable values as vertical lines
#'
#' @param paramSamples the parameters as a vector of point values
#' @param main.title defaults to "parameters"
#' @param lwd line width. defaults to 2.5
#' @param cex.axis axis text size. defaults to .85.
#' @return
#' a plot
#' @export
#'
#' @examples
#' spikePlot()
spikePlot <- function(paramSamples, main.title = "parameters", lwd = 2.5, cex.axis = .85){

    color.assign = function(X){
      colors = rep(0, length(X))
      for (i in 1:length(X)){
        if (sign(X[i]) == 1){
          colors[i] <- "#0073e6BF"
        } else if (sign(X[i]) == -1){
          colors[i] <- "#e60000BF"
        } else{
          colors[i] <- "#000000BF"
        }
      }
      return(colors)
    }
    
    cols = color.assign(paramSamples)
    
    par(mar = c(3, 3, 3, 1))
    plot(paramSamples, type = "h", 
         lwd = lwd, 
         col = cols, 
         yaxt="n", 
         xlab="Index",
         main= main.title, 
         cex.lab=1.12,
         bty="n", 
         xaxt = "n", 
         family = 'serif')
    points(x = which(paramSamples==0), 
           y = rep(0, 
           length(which(paramSamples==0))), 
           pch = 20, 
           cex = .5, 
           col = "#000000A1")
    axis(1, col = NA, tck = 0, family = 'serif',
         cex.axis = cex.axis)
    axis(2, col = NA, tck = 0, family = 'serif',
         cex.axis = cex.axis)
}
