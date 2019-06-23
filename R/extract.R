#' Extract a posterior to a data frame 
#' 
#' This lets you extract the contents of a posterior to a data frame. You can
#' select parameters to keep and drop.
#'
#' @param fit the model
#' @param keeppars parameters to keep. Defaults to NULL (implying anything is kept that's not in droppars)
#' @param droppars parameters to drop. Defaults to  c("log_lik", "ySim", "yTest", "deviance", "omega", "lambda")
#'
#' @return
#' a data frame
#' @export
#'
#' @examples
#' extract(out)
extractPost =  function (fit, keeppars = NULL, droppars = c("log_lik", "ySim", "yTest", "deviance", "omega", "lambda")) 
{

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
  
  as.data.frame(codaObject)
}