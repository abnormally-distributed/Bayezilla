#' Get the posterior odds between a point estimator and a null hypothesis
#'
#' @description Evaluates the density at the null hypothesis vs the density at the estimate, and calculates the
#' posterior odds as P(H0 | Data) / P(theta_hat | Data). This function is not an estimate of the Bayes
#' Factor. The Savage-Dickey Density Ratio Bayes Factors evaluate P(H0 | Data) / P(H0) 
#' and characterize the change in probability between the null hypothesis in the prior and posterior. 
#' Rather, this characterizes the relative posterior probabilities between a null hypothesis and a chosen
#' point estimate and is not reliant on the prior distribution. \cr
#' \cr
#' Smaller values indicate greater evidence against the null. The maximum odds in favor of the null is 1, because
#' if the estimate is the same value as the null hypothesis, the density ratio will consist of density / density = 1.
#' This function also returns a "posterior p-value" which is the density ratio (odds) converted
#' to a probability through the formula p = odds / (1 + odds), which gives the probability
#' of an error if one rejects the point estimate. Hence, if the estimate is at zero and the
#' null hypothesis is also zero, the density ratio = 1, so the p-value here
#' would be 1 / 2 = 0.50. 
#' 
#' Guidelines for interpretation are given in a figure below, however, these are intended as rough 
#' guidelines for interpreting a continuous measure of evidence. The levels of evience given 
#' below are not intended to be significance thresholds, but rather, "mile-markers". 
#' 
#' \if{html}{\figure{postodds.png}{}}
#' \if{latex}{\figure{postodds.png}{}}
#' \cr
#' \cr
#' \if{html}{\figure{postodds2.png}{}}
#' \if{latex}{\figure{postodds2.png}{}}
#'
#'
#' @param fit a runjags or stanfit object
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
        lower = cred_interval(round(x, 5), cred.level = 0.15, method = "HDI")[1]
        upper = cred_interval(round(x, 5), cred.level = 0.15, method = "HDI")[2]
        x = x[-c(which(x < lower), which(x > upper))]
        dens = density(round(x, 5), kernel = "t", n = 1024 * 100)
        dens$x[which.max(dens$y)] 
      }
      estimate = contMode(posterior)
      
    }
    dens.est  = EnvStats::demp(estimate, posterior, discrete = FALSE, density.arg.list = list(kernel = "t", n = 1024 * 100, adjust = 1))
    dens.null = EnvStats::demp(H0, posterior, discrete = FALSE, density.arg.list = list(kernel = "t", n = 1024 * 100, adjust = 1))
    odds = (dens.null + 1) / (dens.est + 1) 
    if (odds > 1){
      odds = 1
    }
    c("odds" = odds, "p" = odds / (1+odds) )
  }
  
    if (is.null(fit)){
      stop("Please provide a runjags or stanfit model object.")
    }
    if (is.null(param) == TRUE)  {
      stop("Please choose a single parameter to test with the 'param' argument.")
    } 
  
    stan <- inherits(fit, "stanfit")
    if (stan == TRUE) {
      paramSampleVec <- as.matrix(fit)
      paramSampleVec = as.vector(paramSampleVec[,which(colnames(paramSampleVec) == param)])
    } 
    else if (class(fit) == "runjags") {
      paramSampleVec = as.vector(combine.mcmc(fit, collapse.chains = TRUE, vars = param))
    }
  
    posterior.odds(paramSampleVec, H0, method)
  
}


