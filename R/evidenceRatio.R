#' Obtain evidence ratios (pseudo-Bayes Factors) comparing two models
#'
#' @description This function is used to compare two models. The output can be interpreted
#  just as one would interpret a Bayes Factor as representing the relative likelihood of one
#' model compared to the other. Instead of using marginal likelihoods this function utilizes
#' either the Watanabe-Akaike Information Criterion or the Leave-One-Out Cross-Validation Information
#' Criterion. The default is to use LOO-IC. See the IC() function for more information on how
#' these are calculated. Models with a larger evidence ratio are favored over the alternative.
#' Posterior model probabilities are also provided for each model to aid interpretation.
#'
#'
#' @param modelA the first model
#' @param modelB the second model
#' @param method one of "LOO" (the default) or "WAIC"
#'
#' @return
#' a data frame
#' @export
#'
#' @examples
#' evidenceRatio()
evidenceRatio = function (modelA, modelB, method = "LOO")
{
  A = IC(modelA)
  B = IC(modelB)

  WAIC.A = log(exp(-0.5 * A[1]))
  WAIC.B = log(exp(-0.5 * B[1]))
  delta = c(A[1], B[1]) - min(c(A[1], B[1]))
  prob.waic = round((exp(-0.5 * delta)) / (sum(exp(-0.5 * delta)) ), digits = 3)

  LOO.A = log(exp(-0.5 * A[2]))
  LOO.B = log(exp(-0.5 * B[2]))
  delta = c(A[2], B[2]) - min(c(A[2], B[2]))
  prob.loo = round((exp(-0.5 * delta)) / (sum(exp(-0.5 * delta)) ), digits = 3)

  WAIC.AB = signif(exp(WAIC.A - WAIC.B), 3)
  WAIC.BA = signif(exp(WAIC.B - WAIC.A), 3)

  LOO.AB = signif(exp(LOO.A - LOO.B), 3)
  LOO.BA = signif(exp(LOO.B - LOO.A), 3)

  if (method == "LOO"){
    mat = matrix(c(LOO.AB, LOO.BA, prob.loo[1], prob.loo[2]), nrow = 4, ncol = 1)
    rownames(mat) = c("Evidence.Ratio_AB", "Evidence.Ratio_BA", "P(A)", "P(B)")
    colnames(mat) = "LOO-IC Statistics"
    return(as.data.frame(mat))
  }
  if (method == "WAIC"){
    mat = matrix(c(WAIC.AB, WAIC.BA, prob.waic[1], prob.waic[2]), nrow = 4, ncol = 1)
    rownames(mat) = c("Evidence.Ratio_AB", "Evidence.Ratio_BA", "P(A)", "P(B)")
    colnames(mat) = "WAIC Statistics"
    return(as.data.frame(mat))
  }
}
