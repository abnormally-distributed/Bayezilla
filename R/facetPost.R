#' Plot Multiple Posterior Distributions 
#'  
#' @description This lets you plot several different posterior histograms at once in a faceted plot. 
#' The number of columns, rows, and several other features can be customized. It is also easy to drop
#' or keep certain parameters so you can view only the ones of interest. \cr 
#' 
#' @param fit a stanfit or runjags object. This can be used as an alternative to paramSampleVec,
#' but you must specify which parameter you would like to plot. If using this argument be sure to type
#' fit = "yourmodel" so that the function knows it is not intended to be a vector.
#' @param breaks a number, function name, or function.
#' @param keeppars variables to keep
#' @param droppars variables to drop. Defaults to c("ySim", "log_lik", "lp__")
#' @param ncol number of columns in the layout
#' @param nrow number of rows in the layout
#' @param col the color scheme. One of "blue", "green" "red", or "purple".
#' @param cred.level The credibility level. Defaults to 90\% (.90).
#' @param method Quantile intervals "QI" (the default), or highest density intervals "HDI"
#' @param est One of "mean" (the default), "median", or "mode". The mode gives the value corresponding to the highest density of
#' a kernel density estimate of the marginal posterior. 
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
facetPost = function(fit, keeppars = NULL, droppars = c("ySim", "log_lik", "lp__"), col = "blue", nrow = 4, ncol = 2, breaks = "dhist", method = "QI", est = "mean", cred.level = .90){
  
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
    plotPost(fit = codaObject[,i], param = paste0(parNames[i]), col = col, breaks = breaks,  
             est = est,
             method = method, 
             cred.level = cred.level)
  }
}
