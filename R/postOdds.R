#' Get the posterior odds between a point estimator and a null hypothesis
#'
#' Evaluates the density at the null hypothesis vs the density at the estimate, and calculates the
#' posterior odds as P(H0 | Data) / P(theta_hat | Data). This function is not an estimate of the Bayes
#' Factor. The Savage-Dickey Density Ratio Bayes Factors evaluate P(H0 | Data) / P(H0) 
#' and characterize the change in probability between the null hypothesis in the prior and posterior. 
#' Rather, this characterizes the relative posterior probabilities between a null hypothesis and a chosen
#' point estimate and is not reliant on the prior distribution. \cr
#' \cr
#' Smaller values indicate greater evidence against the null. A suggestion: Values equal to or larger than 0.33 
#' can be regarded as no evidence against the null, values smaller than 0.33 but larger than 0.2 as weak evidence
#' against the null, 0.2 as moderate evidence against the null, 0.1 as strong evidence against the null, 
#' and smaller than 0.1 as very strong evidence against the null. However these are NOT arbitrary cutoffs
#' as in null-hypothesis-significance testing with p-values, but rather, rough guidelines indexing
#' a continuous measure of evidential strength. 
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
        density(x, kernel = "optcosine", n = 10000)$x[which.max(density(x, kernel = "optcosine", n = 10000)$y)]
      }
      estimate = contMode(posterior)
      
    }
    dens.est  = EnvStats::demp(estimate, posterior, discrete = FALSE)
    dens.null = EnvStats::demp(H0, posterior, discrete = FALSE)
    (dens.null + 2e-16) / (dens.est + 2e-16) 
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


