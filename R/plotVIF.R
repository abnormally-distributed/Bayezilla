#' Variance Inflation Factor Plot
#'
#' @description Calculate the VIF for a least squares or generalized linear model. This is used to
#' diagnose the ill effects of multicollinearity and collinearity in a regression model. This can 
#' be helpful in deciding to keep or drop a variable, or more preferably in some cases, use a regularized model.
#' If the VIFs are all rather low, then using \code{\link[Bayezilla]{glmBayes}} is safe. If some are higher, but the
#' matrix is still full rank, you might wish to use \code{\link[Bayezilla]{apcGlm}}. Otherwise, if the model is not 
#' full rank, not positive definite, or has a very high conditioning number, you may wish to use a ridge regression 
#' estimatior such as \code{\link[Bayezilla]{ridge}}. To obtain information about rank,
#' positive definiteness, and the condition number, use the \code{\link[Bayezilla]{vitals}} function. \cr
#' \cr
#' To use this simply input the formula, data,
#' and family exactly as you would do with the glm() function. A horizontal dash is marked at 5, indicating
#' a common point where many argue the variance inflation is problematic. Some have lower conservative (2) 
#' thresholds and some have higher liberal (10) thresholds, but 5 is one of the more common figures, i.e., 
#' Sheather (2009).
#' \cr \cr
#' What this means is that if a variable is inflating the variance of the estimation by a factor of 
#' 5, the standard error of the corresponding coefficient is 2.236068 higher than it would 
#' be if it were not correlated with other variables. 
#' \cr
#' An example of output: \cr
#' \cr
#' \if{html}{\figure{VIF.png}{}}
#' \if{latex}{\figure{VIF.png}{}}
#' \cr
#' 
#' @references Sheather, Simon (2009). A modern approach to regression with R. New York, NY: Springer. ISBN 978-0-387-09607-0. 
#'
#' @param formula the formula
#' @param data the data
#' @param family the glm family. Defaults to "gaussian"
#' @export
#' @examples
#' plotVIF()
#'
plotVIF = function(formula, data, family="gaussian"){
  
  calculate.vifs =  function (mod, ...) 
  {
    if (any(is.na(coef(mod)))) 
      stop("there are perfectly collinear predictors in the model. cannot calculate VIFs as the variance
           inflation is infinite.")
    v <- vcov(mod)
    assign <- attr(model.matrix(mod), "assign")
    if (names(coefficients(mod)[1]) == "(Intercept)") {
      v <- v[-1, -1]
      assign <- assign[-1]
    }
    else warning("No intercept: vifs may not be appropriate..")
    terms <- labels(terms(mod))
    n.terms <- length(terms)
    if (n.terms < 2) 
      stop("model contains fewer than 2 terms")
    R <- cov2cor(v)
    detR <- det(R)
    result <- matrix(0, n.terms, 3)
    rownames(result) <- terms
    colnames(result) <- c("GVIF", "Df", "GVIF^(1/(2*Df))")
    for (term in 1:n.terms) {
      subs <- which(assign == term)
      result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, 
                                                                         -subs]))/detR
      result[term, 2] <- length(subs)
    }
    if (all(result[, 2] == 1)) 
      result <- result[, 1]
    else result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
    result
  }
  
  model = glm(formula=formula,data=data,family=family)
  vifs = calculate.vifs(model)
  if (isTRUE(is.vector(vifs))) {
    vifs = as.data.frame(vifs)
    colnames(vifs) = "VIF"
  }
  
  vifs = as.data.frame(vifs)
  colnames(vifs) = "VIF"
  vifs$Variable = rownames(vifs)

  ggplot(data = vifs, aes(x=Variable, y=VIF, fill=VIF)) +
  geom_col(fill="steelblue") +
  geom_hline(yintercept=5, linetype="dashed") +
  guides(fill=FALSE)
  
}
