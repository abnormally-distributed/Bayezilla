#' plot a density function over a histogram
#'
#
#' @param dfunc The name of a density function (ie, dnorm, dt, etc)
#' @param args list of the additional argument of the function, ie, list(mu = 0, sigma = 1)
#' @param x the data
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
plotDist = function (dfunc, args = list(), x = NULL, bins = 20, n = 1500, type = "irregular", xname = "\u03C7", 
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