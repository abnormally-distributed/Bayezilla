#' Generate predicted values from posterior samples for new observations
#'
#' @description Note that this just generates the posterior for the predicted means / expected values.
#' for a set of new observations. This is not a posterior predictive distribution generator,
#' which draws samples from the outcome distribution. Hence, credible intervals for 
#' predictions derived in this function represent the credible intervals for the 
#' mean regression line, not the posterior predictive distribution which would have
#' considerably wider intervals.
#'
#' @param samples a matrix of MCMC samples containing the intercept and coefficients
#' @param formula the formula used for the original model
#' @param newdata a data frame of new data
#' @param link the link function used for the original model. For normal, student's t, laplace
#' and similar use "identity". For poisson glms use "log". For binomial models use either
#' "logitProb" if you want the predicted probabilities or "logitBin" if you want the classifications
#' split into the binary outcome, such that predicted probabilities < .50 are 0 and prob. > .50
#' are 1. 
#'
#' @return
#' a matrix
#' @export
#'
#' @examples
#' linpred(posterior, Sepal.Length ~ ., iris[testset,], "identity" )
#' 
linpred = function (samples, formula, newdata, link = "identity") 
{
  xmat = model.matrix(formula, newdata)
  
  if (link == "identity") {
    
    predfun = function(samples, x){
      samples %*% t(xmat)
    }
    
    out = t(apply(samples, 1, function(z) as.vector(predfun(z, xmat))))
    return(out)
  }
  
  if (link == "log") {
    
    predfun = function(samples, x){
      exp(samples %*% t(xmat))
    }
    
    out =t(apply(samples, 1, function(z) as.vector(predfun(z, xmat))))
    return(out)
  }
  
  if (link == "logitProb") {
    
    predfun = function(samples, x){
      mu = exp(samples %*% t(xmat))
      mu / (1 + mu)
    }
    
    out = t(apply(samples, 1, function(z) as.vector(predfun(z, xmat))))
    return(out)
  }
  
  if (link == "logitBin") {
    
    predfun = function(samples, x){
      mu = exp(samples %*% t(xmat))
      mu / (1 + mu)
    }
    
    out = t(apply(samples, 1, function(z) as.vector(predfun(z, xmat))))
    out[which(out < .5)] <- 0
    out[which(out > .5)] <- 1
    return(out)
  }
  
}