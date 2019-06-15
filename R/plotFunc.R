#' plot a function
#'
#
#' @param func The name of a function
#' @param from the lower limit of the x-axis
#' @param to the upper limit of the x-axis
#' @param n the number of x values to be evaluated. Defaults to 1500.
#' @param xname character; argument of the function in func, the name of the x axis
#' @param xlab label of the x-axis
#' @param ylab label of the y-axis
#' @param xlim restrict the range of the function
#' @param args list of the additional argument of the function
#' @param color the color of the line. defaults to "red"
#' @param size the size of the line. defaults to .75
#' @param ... Additional arguments to stat_function()
#' @examples
#' ggfunc(dnorm, from = -8 to = 8)
#' @export
#' 
#' 
plotFunc = 
function (func, from = NULL, to = NULL, n = 1500, xname = "\u03C7", 
          xlab = xname, ylab = NULL, xlim = NULL, args = list(), 
          color = "#005d8bCC", 
          size = .75, 
          ...) 
{
  
  sexpr <- substitute(func)
  if (is.function(func)) {
    funky_function <- func
  }
  
  else {
    stop("please provide a function")
  }
  
  if (is.null(ylab)) 
    ylab <- paste0("ƒ ", "( ", "\u03C7", " )")

  tibble(x = c(from, to)) %>% 
      ggplot(aes(x = x)) + 
      stat_function(fun = funky_function, 
                    n = n, 
                    xlim = xlim, 
                    args = args, 
                    color = color, 
                    size = size, 
                    ...) + 
                    labs(x = xlab, y = ylab)
}