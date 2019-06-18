#' Plot a large number of variable values as vertical lines
#'
#' @param x the model
#' @param estimate "mean" (the default) or "median"
#' @param keeppars The list of specific variables to keep if passing an runjags object. Defaults to "beta"
#' @param droppars list of parameters to exclude
#' @param main.title defaults to "parameters"
#' @param lwd line width. defaults to 2.5
#' @param cex.axis axis text size. defaults to .85.
#' @return
#' a plot
#' @export
#'
#' @examples
#' spikePlot()
spikePlot <- function(x, estimate = "mean", keeppars = "beta", droppars = c("ySim", "log_lik", "lp__"), main.title = "parameters", lwd = 2.5, cex.axis = .85){

  stan <- inherits(x, "stanfit")
  if (stan == TRUE) {
    paramSamples <- as.matrix(x)
    wch = unique(unlist(sapply(droppars, function(z) which(regexpr(z, colnames(paramSamples)) == 1))))
    if (length(wch) != 0){
      paramSamples <- paramSamples[,-wch]
    }
    if (!is.null(keeppars)) {
      wch = unique(unlist(sapply(keeppars, function(z) which(regexpr(z, colnames(paramSamples)) == 1))))
      if (length(wch) != 0){
        paramSamples <- paramSamples[,wch]
      }
    }
  }
  else if (class(x) == "runjags"){
    paramSamples <- runjags::combine.mcmc(x, collapse.chains = TRUE)
    paramSamples <- as.matrix(paramSamples)
    wch = unique(unlist(sapply(droppars, function(z) which(regexpr(z, colnames(paramSamples)) == 1))))
    if (length(wch) != 0){
      paramSamples <- paramSamples[,-wch]
    }
    if (!is.null(keeppars)) {
      wch = unique(unlist(sapply(keeppars, function(z) which(regexpr(z, colnames(paramSamples)) == 1))))
      if (length(wch) != 0){
        paramSamples <- paramSamples[,wch]
      }
    }
  }
  else {
    paramSamples <- as.matrix(x)
    wch = unique(unlist(sapply(droppars, function(z) which(regexpr(z, colnames(paramSamples)) == 1))))
    if (length(wch) != 0){
      paramSamples <- paramSamples[,-wch]
    }
    if (!is.null(keeppars)) {
      wch = unique(unlist(sapply(keeppars, function(z) which(regexpr(z, colnames(paramSamples)) == 1))))
      if (length(wch) != 0){
        paramSamples <- paramSamples[,wch]
      }
    }
  }
  
  if (estimate == "mean"){
    paramSamples <- colMeans(paramSamples)
  } 
  else if (estimate == "median"){
    paramSamples <- apply(paramSamples, 2, median)
  }
  
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
