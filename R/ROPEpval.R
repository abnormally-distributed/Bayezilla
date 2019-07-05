#' Get ROPE p-values
#'
#' @description Calculates the percentage of the posterior distribution that lies within the defined ROPE interval. Odds in favor
#' of the parameter lying within the ROPE are also provided. This utilizes the
#' full posterior distribution, rather than the percentage of the credible interval that overlaps with the ROPE as in the ROPEtest function.
#' This is designed to rule out effects that are smaller than a specified magnitude. 
#' 
#' 
#' Some guidelines for interpreting the values: \cr
#' \cr
#' 
#' \if{html}{\figure{ROPEpval.png}{}}
#' \if{latex}{\figure{ROPEpval.png}{}}
#'
#' @param x a stanfit or runjags object. 
#' @param keeppars The name of the parameters. Can use "beta" to match up with all betas, ie, "beta[1]", "beta[2]", etc. Defaults to c("Intercept", "beta").
#' @param ROPE A numeric vector of two numbers. Defaults to c(-.1, .1), assuming standardized regression coefficients or Cohen's d... Adjust to what is reasonable!
#' @export
#' @examples
#' ROPEpval()
#'
ROPEpval = function (x, 
                 ROPE = c(-.1, .1), 
                 keeppars = c("Intercept", "beta"))
{
  
  droppars = NULL
  stan <- inherits(x, "stanfit")
  if (stan == TRUE) {
    ss <- as.matrix(x)
    wch = unique(unlist(sapply(droppars, function(z) which(regexpr(z, colnames(ss)) == 1))))
    if (length(wch) != 0){
      ss <- ss[,-wch]
    }
    if (!is.null(keeppars)) {
      wch = unique(unlist(sapply(keeppars, function(z) which(regexpr(z, colnames(ss)) == 1))))
      if (length(wch) != 0){
        ss <- ss[,wch]
      }
    }
  }
  else if (class(x) == "runjags"){
    ss <- runjags::combine.mcmc(x, collapse.chains = TRUE)
    ss <- as.matrix(ss)
    wch = unique(unlist(sapply(droppars, function(z) which(regexpr(z, colnames(ss)) == 1))))
    if (length(wch) != 0){
      ss <- ss[,-wch]
    }
    if (!is.null(keeppars)) {
      wch = unique(unlist(sapply(keeppars, function(z) which(regexpr(z, colnames(ss)) == 1))))
      if (length(wch) != 0){
        ss <- ss[,wch]
      }
    }
  }
  else {
    ss <- as.matrix(x)
    wch = unique(unlist(sapply(droppars, function(z) which(regexpr(z, colnames(ss)) == 1))))
    if (length(wch) != 0){
      ss <- ss[,-wch]
    }
    if (!is.null(keeppars)) {
      wch = unique(unlist(sapply(keeppars, function(z) which(regexpr(z, colnames(ss)) == 1))))
      if (length(wch) != 0){
        ss <- ss[,wch]
      }
    }
  }

  names = colnames(ss)
  
  rope = function(x) {
  
  interval <- c(min(x), max(x))
  names(interval) <- c("lower", "upper")

  if (length(ROPE) != 2)
    stop("Argument ROPE needs to be a vector of length two.",
         call. = F)
  if (ROPE[1] > ROPE[2]) {
    tmp <- ROPE[2]
    ROPE[2] <- ROPE[1]
    ROPE[1] <- tmp
  }
  
  original.x = x
  x <- sort(x)
  x.rope <- dplyr::between(x, interval[1], interval[2])
  x <- x[which(x.rope == TRUE)]
  r <- dplyr::between(x, ROPE[1], ROPE[2])
  rope.pct = sum(r)/length(x)
  rope.pct
  }
  
  pvals = sapply(1:ncol(ss), function(x) rope(ss[,x]))
  data.frame("term" = names, "ROPEpvalue" = pvals)
  
}  
  