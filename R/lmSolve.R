#' Least Squares with Moore-Penrose Inverse
#'
#' A fast solver for simple least squares models. Returns coefficients only.
#' Solves with moore-penrose pseudoinversion which is very fast and tolerant
#' to non-positive definite model matrices.
#'
#' @param formula formula
#' @param data data frame
#' @param tolerance tolerance
#'
#' @return
#' a vector
#' @export
#'
#' @examples
#' lmSolve()
#'
lmSolve = function(formula, data, tolerance = 1e-3){
  X = as.matrix(model.matrix(formula, data))
  Y = model.frame(formula, data)[,1]
  betas = as.vector(XtXinv(X, tol = tolerance) %*% t(X) %*% Y)
  names(betas) = colnames(X)
  return(betas)
}

#' pseudoinverse
#'
#' Solve for the inverse of a matrix
#'
#' @param X a matrix
#' @param tol tolerance
#'
#' @return
#' a matrix
#' @export
#'
#' @examples
#' pseudoinverse()
#'
pseudoinverse  =
  function (X, tol = 1e-3)
  {
    X = as.matrix(X)

    positive.svd <- function(m, tol) {
      s <- svd(m)
      if (missing(tol)) {
        tol <- max(dim(m)) * max(s$d) * .Machine$double.eps
      }
      Positive <- s$d > tol
      return(list(d = s$d[Positive], u = s$u[, Positive, drop = FALSE],
                  v = s$v[, Positive, drop = FALSE]))
    }
    fast.svd <- function(m, tol) {
      n <- dim(m)[1]
      p <- dim(m)[2]
      EDGE.RATIO <- 2
      if (n > EDGE.RATIO * p) {
        return(psmall.svd(m, tol))
      }
      else if (EDGE.RATIO * n < p) {
        return(nsmall.svd(m, tol))
      }
      else {
        return(positive.svd(m, tol))
      }
    }
    pseudoinverse <- function(m, tol) {
      msvd <- fast.svd(m, tol)
      if (length(msvd$d) == 0) {
        return(array(0, dim(m)[2:1]))
      }
      else {
        return(msvd$v %*% (1/msvd$d * t(msvd$u)))
      }
    }
    return(as.matrix(pseudoinverse(X, tol = tol)))
  }


#' Obtain the inverse crossproduct of a matrix
#'
#' Returns the pseudoinverse of the crossproduct. Use instead of solve(crossprod(X))
#'
#' @param X the matrix
#' @param tol tolerance
#'
#' @return
#' a matrix
#' @export
#'
#' @examples
#' XtXinv()
#'
XtXinv = function (X, tol = 1e-3)
{
  X = as.matrix(X)
  return(pseudoinverse(crossprod(X), tol = tol))
}
