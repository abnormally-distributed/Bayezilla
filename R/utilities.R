#' Raise a number to a power
#' 
#' @description Ever use the command "pow(x, 3)" only to find that R does not have a pow function? I do this all the time
#' out of habit because both JAGS and Stan use this, and matlab has the similar "power" function (but power does something else in R). 
#' This function exists to save frustration if you use the pow command to raise a number to a power by avoiding 
#' annoying error messages.
#'
#' @param x a number
#' @param lambda a power
#'
#' @return
#' a number
#' @export
#'
#' @examples
#' pow(3, 4)
pow <- function(x, lambda=NULL){
  if (is.null(lambda)){
    cat(crayon::blue(crayon::bold("Please enter a power to which 'x' should be raised as an argument to 'lambda'.")))
  }
  return(x^lambda)
}

#' Square a number
#' 
#' @description Ever use the command "square(x)" only to find that R does not have a "square" function, despite the fact
#' that it has a "sqrt" function? I do this all the time out of habit because I frequently use that when coding in Stan. 
#' This function exists to save me (or you) frustration by not incurring annoying error messages if using square(x) 
#' out of habit.
#'
#' @param x a number
#'
#' @return
#' a number
#' @export
#'
#' @examples
#' square(2)
square <- function(x){
  return(x^2)
}


#' Scale a data frame or matrix
#'
#' @details an improvement of the base R scale function. Unlike R's standard scale function, this
#' allows for factor columns or character columns to be present in the data, and will simply
#' leave those untouched without throwing an error.
#' @param data a data frame or vector
#' @param scale.type one of "sd" (default) for mean centered standard deviation scaled, "medmad" for
#' median centered median absolute deviation scaled, "medsd" for median centered standard
#' deviation scaled, "mean" for mean centered but unscaled, or "med" for median centered
#' but unscaled. Also accepts 1, 2, 3, 4, and 5 respectively as input types.
#' @export
#' @return A data frame
#' @examples
#' scale(data)
#'

scale = function (data, scale.type = 1)
{
  
  if (isTRUE(is.vector(data))) {
    Vector = "YES"
    data = cbind.data.frame(x = data)
  } else {
    Vector = "NO"
    data = as.data.frame(data)
  }
  
  if (scale.type == "sd" || scale.type == 1) {
    Scale1 = function(x) {
      Mean = mean(x)
      Sd = sd(x)
      x = x - Mean
      x = x/Sd
      return(x)
    }
    ind <- sapply(data, is.numeric)
    scaled.data = data
    scaled.data[ind] <- lapply(scaled.data[ind], Scale1)
    
    if (Vector=="YES") {
      return(as.vector(as.matrix(scaled.data)))
    }
    else {
      return(as.data.frame(scaled.data))
    }
  }
  if (scale.type == "medmad" || scale.type == 2) {
    Scale2 = function(x) {
      Median = median(x)
      Mad = mad(x)
      x = x - Median
      x = x/Mad
      return(x)
    }
    ind <- sapply(data, is.numeric)
    scaled.data = data
    scaled.data[ind] <- lapply(scaled.data[ind], Scale2)
    if (Vector=="YES") {
      return(as.vector(as.matrix(scaled.data)))
    }
    else {
      return(as.data.frame(scaled.data))
    }
  }
  if (scale.type == "medsd" || scale.type == 3) {
    Scale3 = function(x) {
      Median = median(x)
      Mad = sd(x)
      x = x - Median
      x = x/Mad
      return(x)
    }
    ind <- sapply(data, is.numeric)
    scaled.data = data
    scaled.data[ind] <- lapply(scaled.data[ind], Scale3)
    if (Vector=="YES") {
      return(as.vector(as.matrix(scaled.data)))
    }
    else {
      return(as.data.frame(scaled.data))
    }
  }
  if (scale.type == "mean" || scale.type == 4) {
    Scale4 = function(x) {
      Mean = mean(x)
      x = x - Mean
      return(x)
    }
    ind <- sapply(data, is.numeric)
    scaled.data = data
    scaled.data[ind] <- lapply(scaled.data[ind], Scale4)
    if (Vector=="YES") {
      return(as.vector(as.matrix(scaled.data)))
    }
    else {
      return(as.data.frame(scaled.data))
    }
  }
  if (scale.type == "med" || scale.type == 5) {
    Scale5 = function(x) {
      Median = median(x)
      x = x - Median
      return(x)
    }
    ind <- sapply(data, is.numeric)
    scaled.data = data
    scaled.data[ind] <- lapply(scaled.data[ind], Scale5)
    if (Vector=="YES") {
      return(as.vector(as.matrix(scaled.data)))
    }
    else {
      return(as.data.frame(scaled.data))
    }
  }
}


#' Compute the precision of a vector
#'
#' An R function to compute the precision of a vector, 1 / var(x)
#'
#' @param x a vector
#'
#' @return
#' the precision
#' @export
#'
#' @examples
#' prec(x)
prec <- function(x) { 1 / var(x)}

#' Compute the posterior model weights from a vector of Information Criteria 
#' 
#' Supply a vector containing AIC, AICc, WAIC, LOO-IC, BIC, DIC, etc. 
#'
#' First a delta score is calculated, which is the difference between each model's
#' score and the minimum model score of the vector. The best model has a delta of zero.
#' 
#' The model weights are calculated from the delta scores as in the equation below:
#'
#' \deqn{w_{m}=\frac{\exp \left(-\frac{1}{2} \Delta_{m}\right)}{\sum_{j=1}^{M} \exp \left(-\frac{1}{2} \Delta_{j}\right)}}:
#' 
#' See: 
#' Wagenmakers, EJ. & Farrell, S. AIC model selection using Akaike weights
#' Psychonomic Bulletin & Review (2004) 11: 192. https://doi.org/10.3758/BF03206482
#' 
#' Burnham, Kenneth P., Anderson, David R. (2002) Model Selection and Multimodel Inference: A Practical Information-Theoretic Approach
#' 2nd Ed. Springer. https://doi.org/10.1007/b97636
#'
#' @param ICs 
#'
#' @return
#' a data frame of delta scores and model weights
#' @export
#'
#' @examples
#' ICs(c(-923.13, -232.45, -896.12))
modelWeights = function (ICs) 
{
  delta = ICs - min(ICs)
  probs = round(exp(-0.5 * delta), 4) / sum(round(exp(-0.5 * delta), 4))
  names = paste0("model", 1:length(probs))
  cbind.data.frame("model" = names, "delta" = delta, "prob[best.model]" = probs)
}


#' Parameterize a Gamma Distribution by Mean and SD
#'
#' @param mean mean
#' @param sd standard deviation
#' @param plot defaults to TRUE
#'
#' @return
#' a plot and/or parameters
#' @export
#'
#' @examples
#' gammaPars(3, 3)
gammaPars = function (mean, sd, plot = TRUE) 
{
  ra = (mean)/(sd^2)
  sh = (mean^2)/(sd^2)
  params = c(sh, ra)
  names(params) = c("shape", "rate")
  if (plot == TRUE) {
    rGamma = function(n, shape = 1, rate = 1){
      return(qgamma(seq(1/n, 1 - 1/n, length.out = n), shape, rate))
    }
    
    hist(rGamma(20000, sh, ra), col = "#20ed80", bty = "n", border = "#07562c", axes = FALSE,
         family = "serif", main = paste0("Gamma( Shape = ", sh, ",  Rate = ", ra, " )"), 
         xlab = NULL, ylab = NULL, breaks = 50)
    axis(1, col = NA, tck = 0, family = "serif")
  }
  return(params)
}