#' Get the posterior odds between a point estimator and a null hypothesis
#'
#' @description Evaluates the density at the null hypothesis vs the density at the estimate, and calculates a
#' posterior p-value as the percentage of the density at the null hypothesis out of the value of density at the  
#' point estimate. If the value of the density at the estimate and null are equivalent, then the probability
#' of the null hypothesis is 1. 
#' 
#' The posterior odds ratio in favor of the null hypothesis is calculated as p / (1 - p).
#' 
#' While similar to the Savage-Dickey density ratio method of calculating Bayes Factors, the posterior odds
#' ratio evaluates P(H0 | Data) / P(H | Data ), while the Bayes Factor evaluates P(H0 | Data) / P(H0).
#' That is, the Bayes Factor evaluates the increase in probability relative to the prior, while the posterior odds
#' is the relative strength of evidence in favor of the null hypothesis. \cr
#' \cr
#' Guidelines for interpretation are given in a figure below, however, these are intended as rough 
#' guidelines for interpreting a continuous measure of evidence. The levels of evience given 
#' below are not intended to be significance thresholds, but rather, "mile-markers". 
#' 
#' \if{html}{\figure{postodds.png}{}}
#' \if{latex}{\figure{postodds.png}{}}
#' \cr
#'
#' @param fit a runjags or stanfit object. Alternatively, a numeric vector can be provided.
#' @param keeppars The name of the parameters. Can use "beta" to match up with all betas, ie, "beta[1]", "beta[2]", etc. Defaults to c("Intercept", "beta").
#' @param H0 a single value or a vector of values for the null hypothesis. Defaults to zero but this is not appropriate
#' for a binomial test. Be sure to pick a reasonable null hypothesis. If only one values is provided for multiple parameters this
#' value will be used for all tests.
#' @param method whether the mean (default), "median", or a kernel density estimate of the "mode" should be used for the hypothesis test
#' @export
#' @examples
#' postOdds()
#' 
postOdds = function(fit, keeppars = c("Intercept", "beta"), H0 = 0, method = "mean"){
  
  posterior.odds = function(posterior, H0, method){
    
    if (method == "mean"){
      estimate = mean(posterior)
    }
    else if (method == "median"){
      estimate = median(posterior)
    }
    else if (method == "mode"){
      
      contMode = function(x) {
        d <- density(x, from = min(x), to = max(x), n = length(x), kernel = "t")
        d$x[which.max(d$y)]
      }
      
      estimate = contMode(posterior)
      
    }
    
    dfun = function(t){
      
      
      fun = approxfun(x = density(posterior, from=min(posterior), to=max(posterior), n=length(posterior), kernel="o")$x, 
                      y = density(posterior, from=min(posterior), to=max(posterior), n=length(posterior), kernel="o")$y)
      d = fun(t)
      if (is.na(d)){
        d <- 0
      }
      return(d)
    }
    
    dens.est  = dfun(estimate)
    dens.null = dfun(H0)
    p = (((dens.null) / (dens.est))  / (1 + ((dens.null) / (dens.est))))
    c("p" = round(p, 3), "odds" = round(p / (1 - p), 3))
  }
  
  if (is.null(fit)){
    stop("Please provide a vector, a runjags, or stanfit model object.")
  }
  
  droppars = NULL
  
  if (is.vector(fit)) {
    paramSampleVec = fit
    posterior.odds(paramSampleVec, H0, method)
  }
  
  else {
    stan <- inherits(fit, "stanfit")
    if (stan == TRUE) {
      ss <- as.matrix(fit)
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
    else if (class(fit) == "runjags"){
      ss <- runjags::combine.mcmc(fit, collapse.chains = TRUE)
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
      ss <- as.matrix(fit)
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
    if (length(H0)==1){
      H0 = rep(H0, ncol(ss))
    }
    pvals = sapply(1:ncol(ss), function(x) round(posterior.odds(ss[,x], H0[x], method = method), 3))
    pvals = as.data.frame(t(pvals))
    rownames(pvals) = names
    pvals
  }
  
}


