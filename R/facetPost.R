#' Plot Multiple Posterior Distributions 
#' 
#' 
#' @param fit a stanfit or runjags object. This can be used as an alternative to paramSampleVec,
#' but you must specify which parameter you would like to plot. If using this argument be sure to type
#' fit = "yourmodel" so that the function knows it is not intended to be a vector.
#' @param pars variables to keep
#' @param droppars variables to drop. Defaults to c("log_lik", "ySim", "yTest", "deviance")
#' @param ncol number of columns in the layout
#' @param nrow number of rows in the layout
#' @param col the color scheme. One of "blue", "green" "red", or "purple".
#' @param cred.level The credibility level. Defaults to 90\% (.90).
#' @param method Quantile Intervals "QI" (the default) or highest density intervals "HDI"
#' @param showMedian Should the median be used instead of the mean?
#' comparison value, the second being the lower limit of the ROPE and the third being the upper limit of the ROPE.
#' 
#' @return A plot
#' @export
#'
#' @examples
#' facetPost()
#'
facetPost = function(fit, pars = NULL, droppars = c("log_lik", "ySim", "yTest", "deviance"), col = "blue",
                     nrow = 4, ncol = 2, method = "QI", showMedian = FALSE, cred.level = .90){
 
  old.par <- par(no.readonly = TRUE) # save default, for resetting...
  on.exit(par(old.par))     #and when we quit the function, restore to original values
  
  codaObject <- extractPost(fit, pars = pars, droppars = droppars)
  parNames = colnames(codaObject)

  par(mar= c(3.25, 3, 3.5, 2) , oma= c(.25,.25,.25,.25) , mgp=c(1.5, 0.125, 0))
  total = nrow*ncol
  layout(matrix(1:total, nrow=nrow))
  
  for (i in 1:length(parNames)){
      plotPost(paramSampleVec = codaObject[,i], param = paste0(parNames[i]), col = col, 
               showMedian = showMedian,
               method = method, 
               cred.level = cred.level)
  }
}