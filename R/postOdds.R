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
#' @param param the name of the parameter to be tested
#' @param H0 a single value or a vector of values for the null hypothesis. Defaults to zero but this is not appropriate
#' for a binomial test. Be sure to pick a reasonable null hypothesis.
#' @param method whether the mean (default), "median", or a kernel density estimate of the "mode" should be used for the hypothesis test
#' @export
#' @examples
#' postOdds()
#' 
postOdds = function(fit, param = NULL, H0 = 0, method = "mean"){
  
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
    p = (dens.null) / (dens.est)
    c("p" = p, "odds" = p / (1 - p))
  }
  
    if (is.null(fit)){
      stop("Please provide a runjags or stanfit model object.")
    }
    if (is.null(param) & !is.vector(fit))  {
      stop("Please choose a single parameter to test with the 'param' argument.")
    } 
  
    stan <- inherits(fit, "stanfit")
    if (is.vector(fit)) {
      paramSampleVec = fit
    }
    else if (stan == TRUE) {
      paramSampleVec <- as.matrix(fit)
      paramSampleVec = as.vector(paramSampleVec[,which(colnames(paramSampleVec) == param)])
    } 
    else if (class(fit) == "runjags") {
      paramSampleVec = as.vector(combine.mcmc(fit, collapse.chains = TRUE, vars = param))
    }
     
  
    posterior.odds(paramSampleVec, H0, method)
  
}


