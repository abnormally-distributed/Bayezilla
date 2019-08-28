#' Minimum Norm Regresion with Moore-Penrose Inverse
#'
#' A fast solver for simple least squares models. Returns coefficients only.
#' Solves with moore-penrose pseudoinversion which is very fast and tolerant
#' to non-positive definite model matrices.
#'
#' @param formula formula
#' @param data data frame
#' @param w optional weights
#' @param tolerance tolerance
#'
#' @return
#' a vector
#' @export
#'
#' @examples
#' lmSolve()
#'
lmSolve = function(formula, data, w = NULL, tolerance = 1e-3){
  X = as.matrix(model.matrix(formula, data))
  Y = model.frame(formula, data)[,1]
  if (is.null(w)){
    betas = as.vector(XtXinv(X, tol = tolerance) %*% t(X) %*% Y)
    names(betas) = colnames(X)
    return(betas)
  } else{
    betas = as.vector(pseudoinverse(t(X) %*% diag(w) %*% X, tol = tolerance) %*% t(X) %*% diag(w) %*% Y)
    names(betas) = colnames(X)
    return(betas)
  }

}

#' Bootstrap minimum norm regression
#'
#' Bootstrap variant of lmSolve. Returns samples of coefficients only.
#' Solves with moore-penrose pseudoinversion which is very fast and tolerant
#' to non-positive definite model matrices.
#'
#' @param formula formula
#' @param data data frame
#' @param iter number of iterations. defaults to 4000.
#' @param tolerance tolerance
#'
#' @return
#' a vector
#' @export
#'
#' @examples
#' bootlmSolve()
#'
bootlmSolve = function(formula, data, iter = 4000, tolerance = 1e-3){
  
  dirichlet_weights <- matrix( rexp(NROW(data) * iter, 1) , ncol = NROW(data), byrow = TRUE)
  dirichlet_weights <- dirichlet_weights / rowSums(dirichlet_weights)
  
  linear_model <- function(formula, data, w, tolerance){
    lmSolve(formula, data, w, tolerance)
  }
  
  stat_call <- quote(linear_model(formula, data, w))
  
  boot_sample <- apply(dirichlet_weights, 1, function(w) {
    eval(stat_call)
  })
  
  return(as.mcmc(t(boot_sample)))
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

    nsmall.svd =
    function (m, tol) 
    {
      B = m %*% t(m)
      s = svd(B, nv = 0)
      if (missing(tol)) 
        tol = dim(B)[1] * max(s$d) * .Machine$double.eps
      Positive = s$d > tol
      d = sqrt(s$d[Positive])
      u = s$u[, Positive, drop = FALSE]
      v = crossprod(m, u) %*% diag(1/d, nrow = length(d))
      return(list(d = d, u = u, v = v))
    }
    
    psmall.svd = 
    function (m, tol) 
    {
      B = crossprod(m)
      s = svd(B, nu = 0)
      if (missing(tol)) 
        tol = dim(B)[1] * max(s$d) * .Machine$double.eps
      Positive = s$d > tol
      d = sqrt(s$d[Positive])
      v = s$v[, Positive, drop = FALSE]
      u = m %*% v %*% diag(1/d, nrow = length(d))
      return(list(d = d, u = u, v = v))
    }
    
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
XtXinv = function (X, robust = TRUE, tol = 1e-3)
{
  X = as.matrix(X)
  return(pseudoinverse(crossprod(X), tol = tol))
}
