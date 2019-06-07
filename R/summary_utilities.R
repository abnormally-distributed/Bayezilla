#' Get a tidy summary of a bayesian model
#'
#' @description Get summary including QI or HDI credible intervals
#'
#' @param x the set of posterior samples to be summarized
#' @param digits number of digits to round estimates to. Defaults to 3.
#' @param estimate.method one of "mean" or "median", or "both". If "mean" or "median" the
#' estimate will be in column title "estimate". If "both" (the default) the mean will be the "estimate" column and the median will be in the "median" column.
#' @param cred.method one of "HDI" (highest density intervals) or "QI" (equal tailed \ quantile intervals).
#' @param cred.level the credibility level. defaults to 0.90
#' @param keep.vars the list of specific variables to keep if passing an runjags object.
#' @param droppars list of parameters to exclude
#' @param ess set to TRUE (default) to include ess
#' @param knitr if you want to return the summary with knitr's kable function set to TRUE. Default is FALSE.
#' @param type output type for kable. Defaults to "markdown"
#' @param ... other arguments to pass to knitr::kable
#' @return a tibble or knitr::kable output
#' @export
#' @examples
#' post_summary()
#'
post_summary = function (x, digits = 3, estimate.method = "both", ess = TRUE, cred.level = 0.90, cred.method = "QI", keep.vars = NA, droppars = NA, knitr = FALSE, type = "markdown", ...)
  {
    stan <- inherits(x, "stanfit")
    if (stan == TRUE) {
      ss <- as.matrix(x)
    }
    else if (class(x) == "runjags"){
      ss <- runjags::combine.mcmc(x, collapse.chains = TRUE, vars = keep.vars)
      ss <- as.matrix(ss)
    }
    else {
      ss <- as.matrix(x)
      ss <- ss[, !colnames(ss) %in% droppars, drop = FALSE]
    }

    if (estimate.method == "mean") {
      estimate.method <- match.arg(estimate.method, c("mean",
                                                      "median"))
      m <- switch(estimate.method, mean = colMeans(ss), median = apply(ss, 2, stats::median))
      ret <- data.frame(estimate = m, std.error = apply(ss, 2, sd))
    }
    else if (estimate.method == "median") {
      estimate.method <- match.arg(estimate.method, c("mean",
                                                      "median"))
      m <- switch(estimate.method, mean = colMeans(ss), median = apply(ss, 2, stats::median))
      ret <- data.frame(estimate = m, std.error = apply(ss, 2, sd))
    }
    else if (estimate.method == "both"){
      exp.val <- colMeans(ss)
      q50 <- apply(ss, 2, stats::median)
      ret <- data.frame(estimate = exp.val, median = q50,
                        std.error = apply(ss, 2, sd))
    }

    digits = digits
    cred.method <- match.arg(cred.method, c("QI", "HDI"))
    SS = as.data.frame(ss)
    ci <- switch(cred.method, QI = t(apply(SS, 2, function(s)
      Bayezilla:::cred_interval(s, method = "QI", cred.level = cred.level))),
      HDI = t(apply(SS,2, function(s) cred_interval(s, method = "HDI", cred.level = cred.level))))
    ci = as.data.frame(ci)
    colnames(ci) = c("cred.low", "cred.high")
    ret <- data.frame(ret, ci)

    ret = Bayezilla:::post_tibble(ret)
    ret$term = factor(ret$term, levels = ret$term)

    if (estimate.method == "both"){
      ret = ret %>% mutate(estimate = round(estimate, digits = digits),
                           median = round(median, digits = digits),
                           std.error =  round(std.error, digits = digits),
                           cred.low =  round(cred.low, digits = digits),
                           cred.high =  round(cred.high, digits = digits))
    } else {
      ret = ret %>% mutate(estimate = round(estimate, digits = digits),
                           std.error =  round(std.error, digits = digits),
                           cred.low =  round(cred.low, digits = digits),
                           cred.high =  round(cred.high, digits = digits))
    }
    if (ess == "TRUE"){
      ess = apply(ss, 2, coda::effectiveSize)
      ret$ess <- ess
    }

    if (knitr == TRUE){
      knitr::kable(ret, format = type, digits = digits, ...)
    }
    else {
      return(ret)
    }
  }

#' sign_probs
#'
#' Summarize the sign probabilities of the coefficients in a runjags object
#'
#' @param model the model
#' @param labs variable labels
#'
#' @return
#' a data frame
#' @export
#' @examples
#' sign_probs()
sign_probs = function(model, labs = NULL){
  x = as.matrix(combine.mcmc(model,collapse.chains = TRUE))
  x = x[, grep(colnames(x), pattern =  "beta")]
  x = sign(x)

  neg = x
  neg[which(sign(neg) == 1)] <- 0
  neg = colSums(abs(neg)) / nrow(neg)

  pos = x
  pos[which(sign(pos) == -1)] <- 0
  pos = colSums(abs(pos)) / nrow(pos)

  zeros = 1 - (colSums(abs(x)) / nrow(x))

  out = data.frame("p(b < 0)" = neg, "p(b > 0)" = pos, "p(b = 0)" = zeros)
  if (is.null(labs) == FALSE){
    rownames(out) <- labs
  }
  out = rownames_to_column(out)
  out = as.data.frame(out)
  colnames(out) = c("variable", "p(b < 0)", "p(b > 0)", "p(b = 0)")
  z = rep(" ", nrow(out))
  p = which(out$`p(b = 0)` < .50)
  z[p] <- "*"
  p = which(out$`p(b = 0)` < .20)
  z[p] <- "**"
  out = cbind.data.frame(out, z)
  colnames(out) = c("variable", "p(b < 0)", "p(b > 0)", "p(b = 0)", " ")
  return(out)
}


#' Convert posterior summary to tibble
#' @keywords internal
#' @examples
#' post_tibble()
post_tibble = function (x, newnames = NULL, newcol = "term")
{
  unrowname = function (x)
  {
    rownames(x) <- NULL
    x
  }

  if (!is.null(newnames) && length(newnames) != ncol(x)) {
    stop("newnames must be NULL or have length equal to number of columns")
  }
  if (all(rownames(x) == seq_len(nrow(x)))) {
    ret <- data.frame(x, stringsAsFactors = FALSE)
    if (!is.null(newnames)) {
      colnames(ret) <- newnames
    }
  }
  else {
    ret <- data.frame(...new.col... = rownames(x), unrowname(x),
                      stringsAsFactors = FALSE)
    colnames(ret)[1] <- newcol
    if (!is.null(newnames)) {
      colnames(ret)[-1] <- newnames
    }
  }
  dplyr::as_tibble(ret)
}

#' Calculate Equal Tailed or Highest Density Interval
#'
#'
#' @param object vector of posterior samples
#' @param cred.level the confidence level. defaults to  .90
#' @param method type of credible interval. One of "QI" or default "HDI".
#' @export
#' @examples
#' cred_interval()

cred_interval = function (object, cred.level = .90, method="QI")
{
  if (method=="HDI"){
    object = as.vector(as.matrix(object))
    result <- c(NA_real_, NA_real_)
    if (is.numeric(object)) {
      attributes(object) <- NULL
      x <- sort.int(object, method = "quick")
      n <- length(x)
      if (n > 0) {
        exclude <- n - floor(n * cred.level)
        low.poss <- x[1:exclude]
        upp.poss <- x[(n - exclude + 1):n]
        best <- which.min(upp.poss - low.poss)
        result <- c(low.poss[best], upp.poss[best])
        names(result) <- c("hdi.low", "hdi.high")
      }
    }
  } else if (method=="QI"){
    object = as.vector(as.matrix(object))
    alpha = 1-cred.level
    upper = cred.level + (alpha/2)
    lower = alpha/2
    result <- quantile(object,probs=c(lower,upper), type=2)
    names(result)  <- c("cred.low", "cred.high")
  }
  return(result)
}



#' Extract Point Estimate of Correlation Matrix 
#'
#' @param fit the stanfit or runjags object. Must contain only variables with the name "Rho[,]".
#'
#' @return
#' a matrix
#' @export
#'
#' @examples
#' extractCormat(out)
extractCormat = function (fit) 
{
  stan <- inherits(fit, "stanfit")
  if (stan == TRUE) {
    post <- as.matrix(fit)
  }
  else if (class(fit) == "runjags"){
    post <- combine.mcmc(fit, collapse.chains = TRUE)
    post <- as.matrix(post)
  }
  P = sqrt(ncol(post))
  estimates = colMeans(post)
  post = matrix(estimates, P, P)
  return(post)
}