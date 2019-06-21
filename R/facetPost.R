#' Plot Multiple Posterior Distributions 
#'  
#' @description This lets you plot several different posterior histograms at once in a faceted plot. 
#' The number of columns, rows, and several other features can be customized. It is also easy to drop
#' or keep certain parameters so you can view only the ones of interest. \cr 
#' \cr
#' \cr
#' A unique feature of this function is that by default irregular sized bins are used. Why? Although
#' this may appear unsightly to you at first, it actually provides more information. There are two general
#' classes of histogram binning: equal-width and equal-area histograms. The former oversmooths in
#' regions of high density, and is poor at identifying sharp peaks and multimodality. By contrast, the latter variety
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
#' @references Denby, L., & Mallows, C. (2009). Variations on the Histogram. Journal of Computational and Graphical Statistics, 18(1), 21–31. doi:10.1198/jcgs.2009.0002 \cr
#' 
#' @param fit a stanfit or runjags object. This can be used as an alternative to paramSampleVec,
#' but you must specify which parameter you would like to plot. If using this argument be sure to type
#' fit = "yourmodel" so that the function knows it is not intended to be a vector.
#' @param type  Either "irregular" (default) which uses uneven sized bins proportional to probability density or "equal" for even sized bins. 
#' @param bins Adjust the number of bins in the histogram. 
#' @param keeppars variables to keep
#' @param droppars variables to drop. Defaults to c("ySim", "log_lik", "lp__")
#' @param ncol number of columns in the layout
#' @param nrow number of rows in the layout
#' @param col the color scheme. One of "blue", "green" "red", or "purple".
#' @param cred.level The credibility level. Defaults to 90\% (.90).
#' @param method Quantile Intervals "QI" (the default) or highest density intervals "HDI"
#' @param showMedian Should the median be used instead of the mean?
#' @param comparison value, the second being the lower limit of the ROPE and the third being the upper limit of the ROPE.
#' 
#' @return A plot
#' @export
#'
#' @examples
#' facetPost()
#' 
#' @details 
#' 
#' An example of output: \cr
#' \cr
#' \if{html}{\figure{facetPost.png}{}}
#' \if{latex}{\figure{facetPost.png}{}}
#'
facetPost = function(fit, keeppars = NULL, droppars = c("ySim", "log_lik", "lp__"), col = "blue",  nrow = 4, ncol = 2, bins = 20, type = "irregular", method = "QI", showMedian = FALSE, cred.level = .90){
  
  old.par <- par(no.readonly = TRUE) # save default, for resetting...
  on.exit(par(old.par))     #and when we quit the function, restore to original values
  
  stan <- inherits(fit, "stanfit")
  if (stan == TRUE) {
    codaObject <- as.matrix(fit)
    wch = unique(unlist(sapply(droppars, function(z) which(regexpr(z, colnames(codaObject)) == 1))))
    if (length(wch) != 0){
      codaObject <- codaObject[,-wch]
    }
    if (!is.null(keeppars)) {
      wch = unique(unlist(sapply(keeppars, function(z) which(regexpr(z, colnames(codaObject)) == 1))))
      if (length(wch) != 0){
        codaObject <- codaObject[,wch]
      }
    }
  }
  else if (class(fit) == "runjags"){
    codaObject <- runjags::combine.mcmc(fit, collapse.chains = TRUE)
    codaObject <- as.matrix(codaObject)
    wch = unique(unlist(sapply(droppars, function(z) which(regexpr(z, colnames(codaObject)) == 1))))
    if (length(wch) != 0){
      codaObject <- codaObject[,-wch]
    }
    if (!is.null(keeppars)) {
      wch = unique(unlist(sapply(keeppars, function(z) which(regexpr(z, colnames(codaObject)) == 1))))
      if (length(wch) != 0){
        codaObject <- codaObject[,wch]
      }
    }
  }
  else {
    codaObject <- as.matrix(fit)
    wch = unique(unlist(sapply(droppars, function(z) which(regexpr(z, colnames(codaObject)) == 1))))
    if (length(wch) != 0){
      codaObject <- codaObject[,-wch]
    }
    if (!is.null(keeppars)) {
      wch = unique(unlist(sapply(keeppars, function(z) which(regexpr(z, colnames(codaObject)) == 1))))
      if (length(wch) != 0){
        codaObject <- codaObject[,wch]
      }
    }
  }
  
  parNames = colnames(codaObject)
  
  par(mar= c(3.25, 3, 3.5, 2) , oma= c(.25,.25,.25,.25) , mgp=c(1.5, 0.125, 0))
  total = nrow*ncol
  layout(matrix(1:total, nrow=nrow))
  
  for (i in 1:length(parNames)){
    plotPost(paramSampleVec = codaObject[,i], param = paste0(parNames[i]), col = col, type = type, bins = bins,  
             showMedian = showMedian,
             method = method, 
             cred.level = cred.level)
  }
}
