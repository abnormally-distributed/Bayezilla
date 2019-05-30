
spikePlot <- function(b, main.title = "parameters"){
  

  
color.assign = function(X){
  colors = rep(0, length(X))
  for (i in 1:length(X)){
    if (sign(X[i]) == 1){
      colors[i] <- "#0073e6CC"
    } else if (sign(X[i]) <= 0){
      colors[i] <- "#e60000CC"
    }
  }
  return(colors)
}

cols = color.assign(b)

par(mar = c(3, 3, 3, 1))
plot(b, type = "h", 
     lwd = 2.25, 
     col = cols, 
     yaxt="n", 
     xlab="Index",
     main= main.title, 
     cex.lab=1.12,
     bty="n", 
     xaxt = "n", 
     family = 'serif')
axis(1, col = NA, tck = 0, family = 'serif',
     cex.axis = .75)
axis(2, col = NA, tck = 0, family = 'serif',
     cex.axis = .75)
}