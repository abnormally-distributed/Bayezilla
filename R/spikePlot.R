#' Plot a large number of variable values as vertical lines
#'
#' @param paramSamples the parameters as a vector of point values
#' @param main.title defaults to blank ""
#' @param lwd line width. defaults to 2.5
#'
#' @return
#' a plot
#' @export
#'
#' @examples
#' spikePlot()
spikePlot <- function(paramSamples, main.title = "", lwd = 2.5){

  
color.assign = function(X){
  colors = rep(0, length(X))
  for (i in 1:length(X)){
    if (sign(X[i]) == 1){
      colors[i] <- "#0073e6CC"
    } else if (sign(X[i]) == -1){
      colors[i] <- "#e60000CC"
    }
    else{
      colors[i] <- "black"
    }
  }
  return(colors)
}

cols = color.assign(paramSamples)
par(mar = c(2.5, 2.5, 3, 1))
plot(paramSamples, type = "h", 
     lwd = lwd, 
     col = cols, 
     yaxt="n", 
     xlab="Variable Index",
     main= main.title, 
     cex.lab=1.12,
     bty="n", 
     xaxt = "n", 
     family = 'serif')
axis(1, col = NA, tck = 0, family = 'serif',
     cex.axis = .85)
axis(2, col = NA, tck = 0, family = 'serif',
     cex.axis = .85)
}