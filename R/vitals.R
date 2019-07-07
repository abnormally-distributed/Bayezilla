#' Quickly calculate the rank, condition, positive definiteness, and variance inflation factors for a model
#'
#' @param formula a model formula
#' @param data a data frame
#' @param matrix if not using formula, you can provide a model matrix. The intercept column will be removed
#' if detected.
#' @param y the vector of y values if providing the model matrix
#' @param family The glm family being used (required for the VIFs calculation). One of "gaussian" (the default), "binomial", "Gamma", "inverse.gaussian", "poisson", "quasi", "quasibinomial", or "quasipoisson"
#' 
#' @description This lets you quickly and easily calculate a few important things to know about a model matrix. \cr
#' \cr
#' The first is the column rank of a matrix. If this returns a value other than the number of columns full.rank will display FALSE.\cr
#' \cr
#' The second is the ratio the maximum to minimum singular value of the model matrix, which gives the conditioning number, C. 
#' When a  matrix has a bad condition number that indicates the model is ill-conditioned. This means 
#' a linear regression will be extremely sensitive to even the smallest errors and give untrustworthy results.
#' In the worst case, if C is infinite the model matrix is singular and shrinkage methods must be used to 
#' obtain a solution at all. Alternatively, the conditioning number can be defined as the ratio of the minimum to maximum singular 
#' value of the model matrix, Cinv. Here, larger is better.  \cr
#' \cr
#' A rough estimate of how many digits the estimated y values will have can be manually calculated via the formula mean(D) - log(C),
#' where mean(D) is the average number of decimal digits in the entries of the vector y. 
#' \cr
#' The third 'vital sign' is a check for positive definiteness of the covariance matrix, cov(Matrix).
#' If the covariance matrix has zero or negative eigenvalues, it will fail the positive definiteness check.
#' \cr
#' \cr
#' Also caculated are the variance inflation factors. These indicate the factor by which the standard error is inflated for
#' a variable due to correlation with other variables. A common threshold for a VIF being bad is 5. Other commonly used
#' thresholds are 3, 6, 8, and, 10. However, these should be interpreted with some degree of caution. If the model has
#' a sufficiently low conditioning number and passes all of the other checks a decent fit may still be obtained (O’Brien, 2007). 
#' \cr
#' \cr
#' A sufficiently bad conditioning number, a lack of positive-definiteness, or
#' lack of full rank can, in the worst case, mean the model matrix may be non-invertible and OLS or GLMs will not work. 
#' Shrinkage methods will have to be used to fit the model.  Fortunately, this package provides a great number of regularized regression models. 
#' Otherwise, sufficiently bad values can result in 
#' the model being "ill-posed", which means one of the three following conditions of a well-posed mathematical problem
#' are violated: \cr
#' \cr
#' 1. A solution exists \cr
#' 2. A unique solution exists \cr
#' 3. The output of a function changes continuously with the input(s) \cr
#' \cr
#' \cr
#' @return
#' a list
#' @export
#' @examples
#' vitals(matrix = model.matrix(Sepal.Width ~ ., iris)[,-1], y = iris$Sepal.Width)
#' vitals(Sepal.Width ~ ., iris)
#' 
#' @references O’brien, R.M.(2007) A Caution Regarding Rules of Thumb for Variance Inflation Factors Qual Quant 41: 673. https://doi.org/10.1007/s11135-006-9018-6
#' 
vitals = function(formula = NULL, data = NULL, matrix = NULL, y = NULL, family = "gaussian"){
  
  sigdig = function(x) 
  { 
    str<-as.character(x) 
    if(is.na(strsplit(str,"\\.")[[1]][2])) return(0) 
    else return(nchar(strsplit(str,"\\.")[[1]][2]))   
  } 
  
  calculate.vifs =  function (mod, ...) 
  {
    if (any(is.na(coef(mod)))) 
      stop("There are perfectly collinear predictors in the model. Cannot calculate VIFs as the variance
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
    else result <- result[, 1]
    result
  }
  
  if (!is.null(matrix)){
    matrix = as.matrix(matrix)
    if(length(unique(matrix[,1])) == 1){
      cat(crayon::magenta("The first column is a constant. The first column has been removed since it must be the intercept column."))
      matrix = matrix[,-1] 
    }
    d = svd(matrix)$d
    rank = suppressWarnings(suppressMessages(Matrix::rankMatrix(matrix, method = "qrLINPACK")))
    full.rank = isTRUE(ncol(matrix) - rank == 0)
    
    if (isTRUE(full.rank)){
      if (!is.null(y)){
        model = glm(formula= y ~ . , data= cbind.data.frame(y = y, matrix),family=family)
        vifs = calculate.vifs(model)
      }
      else{
        stop("Please provide a vector of y values with the model matrix.")
      }
    } else {
      vifs = "Cannot calculate VIFs, model is not full rank."
    }
    
    list(
      data.frame(
        rank =  rank,
        full.rank = full.rank,
        C = max(d) / min(d),
        Cinv = min(d) / max(d), 
        positive.definite = fBasics::isPositiveDefinite(cov(matrix))
      ),
      VIFs = vifs
    )
  }
  
  if (!is.null(formula)){
    if (is.null(data)){
      stop("You must provide a data frame to the 'data' argument when using a formula")
    }
    matrix = as.matrix(model.matrix(formula, data))[,-1]
    y = model.frame(formula, data)[,1]
    d = svd(matrix)$d
    rank = suppressWarnings(suppressMessages(Matrix::rankMatrix(matrix, method = "qrLINPACK")))
    full.rank = isTRUE(ncol(matrix) - rank == 0)

    if (isTRUE(full.rank)){
        model = model = glm(formula=formula,data=data,family=family)
        vifs = calculate.vifs(model)
    } else {
      vifs = "Cannot calculate VIFs, model is not full rank."
    }

    list(
      data.frame(
        rank =  rank,
        full.rank = full.rank,
        C = max(d) / min(d),
        Cinv = min(d) / max(d), 
        positive.definite = fBasics::isPositiveDefinite(cov(matrix))
      ),
      VIFs = vifs
    )
  }
}
