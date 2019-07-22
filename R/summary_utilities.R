#' Get a tidy summary of a bayesian model
#'
#' @description Get a detailed summary of your posterior distribution
#' in the tidy data frame format.
#'
#' @param x the set of posterior samples to be summarized
#' @param digits number of digits to round estimates to. Defaults to 2.
#' @param estimate.method one of "mean" or "median", or "both". If "mean" or "median" the
#' estimate will be in column title "estimate". If "both" (the default) the mean will be the "estimate" column and the median will be in the "median" column.
#' @param cred.method one of "HDI" (highest density intervals) or "QI" (equal tailed \ quantile intervals).
#' @param cred.level the credibility level. defaults to 0.90
#' @param keeppars the list of specific variables to keep
#' @param droppars list of parameters to exclude
#' @param ess set to TRUE to include ess
#' @param knitr if you want to return the summary with knitr's kable function set to TRUE. Default is TRUE.
#' @param type output type for kable. Defaults to "markdown"
#' @param ... other arguments to pass to knitr::kable
#' @return a tibble or knitr::kable output
#' @export
#' @examples
#' post_summary()
#'
post_summary = function (x, digits = 2, estimate.method = "both", ess = FALSE, cred.level = 0.90, cred.method = "QI", keeppars = NULL, droppars = c("ySim", "log_lik", "lp__"), knitr = TRUE, type = "markdown", ...)
  {
    stan <- inherits(x, "stanfit")
    if (stan == TRUE) {
      ss <- as.matrix(x)
      wch = unique(unlist(sapply(droppars, function(z) which(regexpr(z, colnames(ss)) == 1))))
      if (length(wch) != 0){
        ss <- ss[,-wch]
      }
      if (!is.null(keeppars)) {
        wch = unique(unlist(sapply(keeppars, function(z) which(regexpr(z, colnames(ss)) == 1))))
        if (length(wch) != 0){
          ss <- ss[,wch]
        }
      }
    }
    else if (class(x) == "runjags"){
      ss <- runjags::combine.mcmc(x, collapse.chains = TRUE)
      ss <- as.matrix(ss)
      wch = unique(unlist(sapply(droppars, function(z) which(regexpr(z, colnames(ss)) == 1))))
      if (length(wch) != 0){
        ss <- ss[,-wch]
      }
      if (!is.null(keeppars)) {
        wch = unique(unlist(sapply(keeppars, function(z) which(regexpr(z, colnames(ss)) == 1))))
        if (length(wch) != 0){
          ss <- ss[,wch]
        }
      }
    }
    else {
      ss <- as.matrix(x)
      wch = unique(unlist(sapply(droppars, function(z) which(regexpr(z, colnames(ss)) == 1))))
      if (length(wch) != 0){
        ss <- ss[,-wch]
      }
      if (!is.null(keeppars)) {
        wch = unique(unlist(sapply(keeppars, function(z) which(regexpr(z, colnames(ss)) == 1))))
        if (length(wch) != 0){
          ss <- ss[,wch]
        }
      }
    }

    if (estimate.method == "mean") {
      estimate.method <- match.arg(estimate.method, c("mean",
                                                      "median"))
      m <- switch(estimate.method, mean = colMeans(ss), median = apply(ss, 2, function(x) stats::median(zapsmall(x, digits = 3))))
      ret <- data.frame(estimate = m, std.error = apply(ss, 2, sd))
    }
    else if (estimate.method == "median") {
      estimate.method <- match.arg(estimate.method, c("mean",
                                                      "median"))
      m <- switch(estimate.method, mean = colMeans(ss), median = apply(ss, 2, function(x) stats::median(zapsmall(x, digits = 3))))
      ret <- data.frame(estimate = m, std.error = apply(ss, 2, sd))
    }
    else if (estimate.method == "both"){
      exp.val <- colMeans(ss)
      q50 <- apply(ss, 2, function(x) stats::median(zapsmall(x, digits = 3)))
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

#' Calculate Equal Tailed Quantile Intervals or Highest Density Intervals
#' 
#' @param object vector of posterior samples
#' @param cred.level the confidence level. defaults to  .90
#' @param method type of credible interval. One of  "HDI", or default "QI".
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
#' @param fit the stanfit or runjags object.
#' @param type either "Rho" or "parRho"  which extracts columns labeled "Rho" or "parRho" respectively. 
#' @return
#' a matrix
#' @export
#'
#' @examples
#' extractCormat(out)
extractCormat = function (fit, type = "Rho") 
{
  stan <- inherits(fit, "stanfit")
  if (stan == TRUE) {
    post <- extractPost(fit, keeppars = type)
  }
  else if (class(fit) == "runjags"){
    post <- extractPost(fit, keeppars = type)
  }
  P = sqrt(ncol(post))
  estimates = colMeans(post)
  post = matrix(estimates, P, P)
  diag(post) <- 1
  return(post)
}



#' Get the marginal modes of a posterior distribution
#'
#' @param x set of posterior samples to be summarized in a runjags or stanfit object
#' @param keeppars list of specific variables to keep
#' @param droppars list of parameters to exclude
#' @return a vector of modes
#' @export
#'
#' @examples
#' marginalModes(fit)
marginalModes = function(x, keeppars = NULL, droppars = NULL){
  
  stan <- inherits(x, "stanfit")
  if (stan == TRUE) {
    ss <- as.matrix(x)
    wch = unique(unlist(sapply(droppars, function(z) which(regexpr(z, 
                                                                   colnames(ss)) == 1))))
    if (length(wch) != 0) {
      ss <- ss[, -wch]
    }
    if (!is.null(keeppars)) {
      wch = unique(unlist(sapply(keeppars, function(z) which(regexpr(z, 
                                                                     colnames(ss)) == 1))))
      if (length(wch) != 0) {
        ss <- ss[, wch]
      }
    }
  }
  else if (class(x) == "runjags") {
    ss <- runjags::combine.mcmc(x, collapse.chains = TRUE)
    ss <- as.matrix(ss)
    wch = unique(unlist(sapply(droppars, function(z) which(regexpr(z, 
                                                                   colnames(ss)) == 1))))
    if (length(wch) != 0) {
      ss <- ss[, -wch]
    }
    if (!is.null(keeppars)) {
      wch = unique(unlist(sapply(keeppars, function(z) which(regexpr(z, 
                                                                     colnames(ss)) == 1))))
      if (length(wch) != 0) {
        ss <- ss[, wch]
      }
    }
  }
  else {
    ss <- as.matrix(x)
    wch = unique(unlist(sapply(droppars, function(z) which(regexpr(z, 
                                                                   colnames(ss)) == 1))))
    if (length(wch) != 0) {
      ss <- ss[, -wch]
    }
    if (!is.null(keeppars)) {
      wch = unique(unlist(sapply(keeppars, function(z) which(regexpr(z, 
                                                                     colnames(ss)) == 1))))
      if (length(wch) != 0) {
        ss <- ss[, wch]
      }
    }
  }
  
  contMode = function(x) {
    d <- density(x, n = length(x), kernel = "triangular")
    d$x <- density(x, from = min(x), to = max(x), n = length(x), kernel = "triangular")$x
    round(d$x[which.max(d$y)], 5)
  }
  
  modes = as.vector(apply(ss, 2, contMode))
  names(modes) = colnames(ss)
  return(modes)
}


#' Approximate the maximum joint posterior density estimate 
#' 
#' @description This function calculates the density estimate for each column of a 
#' joint posterior distribution, then takes the logarithm of the estimate densities. The log-
#' densities are then summed to give the row-wise log-posterior probability function. The 
#' row with the highest log-posterior density is returned which corresponds to the monte carlo
#' sample with the highest joint probability, or MAP estimate. Note, however, this is an approximation
#' to the MAP in the sense that MAP estimation proper is done through optimization.
#'
#' @param x the set of posterior samples to be summarized in an runjags or stanfit object.
#' @param keeppars the list of specific variables to keep. Defaults to NULL. 
#' @param droppars list of parameters to exclude from the calculation of the density. The default is c("ySim", "log_lik", "lp__", "Deviance", "BIC", "delta"),
#' but this should be adjusted to make sure that derived and generated quantities such as predictions, deviance, indicator variables, and so forth --
#' anything not directly part of the posterior distribution of the model parameters -- are removed. Note, however, that the function will still
#' return all quantities in the model. The droppars argument only applies to the density estimation itself. 
#' @return a vector of modes
#' @export
#'
#' @examples
#' jointMode(fit)
jointMode = function(x, keeppars = NULL, droppars = c("ySim", "log_lik", "lp__", "Deviance", "delta", "BIC")){
  
  stan <- inherits(x, "stanfit")
  
  if (stan == TRUE) {
    ss <- as.matrix(x)
    ss1 <- ss
    wch = unique(unlist(sapply(droppars, function(z) which(regexpr(z,  colnames(ss)) == 1))))
    if (length(wch) != 0) {
      ss <- ss[, -wch]
    }
    if (!is.null(keeppars)) {
      wch = unique(unlist(sapply(keeppars, function(z) which(regexpr(z,  colnames(ss)) == 1))))
      if (length(wch) != 0) {
        ss <- ss[, wch]
      }
    }
  }
  else if (class(x) == "runjags") {
    ss <- runjags::combine.mcmc(x, collapse.chains = TRUE)
    ss <- as.matrix(ss)
    ss1 <- ss
    wch = unique(unlist(sapply(droppars, function(z) which(regexpr(z, colnames(ss)) == 1))))
    if (length(wch) != 0) {
      ss <- ss[, -wch]
    }
    if (!is.null(keeppars)) {
      wch = unique(unlist(sapply(keeppars, function(z) which(regexpr(z, colnames(ss)) == 1))))
      if (length(wch) != 0) {
        ss <- ss[, wch]
      }
    }
  }
  else {
    ss <- as.matrix(x)
    ss1 <- ss
    wch = unique(unlist(sapply(droppars, function(z) which(regexpr(z, 
                                                                   colnames(ss)) == 1))))
    if (length(wch) != 0) {
      ss <- ss[, -wch]
    }
    if (!is.null(keeppars)) {
      wch = unique(unlist(sapply(keeppars, function(z) which(regexpr(z, 
                                                                     colnames(ss)) == 1))))
      if (length(wch) != 0) {
        ss <- ss[, wch]
      }
    }
  }
  
  logdensity = function(x) {
    d <- density(zapsmall(round(x, 3), 2), from = min(x), to = max(x), n = length(x), kernel = "triangular")
    lp = log(d$y) 
    #dx = d$x
    #cbind.data.frame(lp = lp, x = dx)
    round(lp, 4)
  }
  
  logpost = apply(ss, 2, function(x) logdensity(x))
  round(ss1[which.max(rowSums(logpost)) , ], 5)
}



#' Approximate the maximum joint posterior density estimate 
#' 
#' @description This function calculates the kernel density estimate for each column of a 
#' joint posterior distribution, then takes the logarithm of the estimate densities. The log-
#' densities are then summed to give the row-wise log-posterior probability function. The 
#' row with the highest log-posterior density is returned which corresponds to the monte carlo
#' sample with the highest joint probability, or MAP estimate. Note, however, this is an approximation
#' to the MAP in the sense that MAP estimation proper is done through optimization. \cr
#' \cr
#' Frankly this function is included largely for didactic purposes. MAP estimation, whether by approximating via
#' finding the MCMC iteration with the highest joint log-probability or directly optimizing an objective function,
#' will be often poor with higher dimensional problems. This is due to the fact that the mode of a multidimensional
#' function may be no where near the expected values of the marginal distributions. This stems from the curse of
#' dimensionality. \cr
#'
#' @param x the set of posterior samples to be summarized in an runjags or stanfit object.
#' @param keeppars the list of specific variables to keep. Defaults to NULL. 
#' @param droppars list of parameters to exclude from the calculation of the density. The default is c("ySim", "log_lik", "Deviance", "BIC", "delta"),
#' but this should be adjusted to make sure that derived and generated quantities such as predictions, deviance, indicator variables, and so forth --
#' anything not directly part of the posterior distribution of the model parameters -- are removed. Note, however, that the function will still
#' return all quantities in the model. The droppars argument only applies to the density estimation itself. 
#' @return a vector of modes
#' @export
#'
#' @examples
#' jointMode(fit)
jointMode = function(x, keeppars = NULL, droppars = c("ySim", "log_lik", "Deviance", "delta", "BIC")){
  
  stan <- inherits(x, "stanfit")
  
  if (stan == TRUE){
    ss <- as.matrix(x)
    if(isTRUE(suppressWarnings(any(colnames(ss), "lp__")))){
      ssp = as.data.frame(ss)
      stop(return(ss[which.max(ssp$`lp__`), ]))
    }
    else {
    ss1 <- ss
    wch = unique(unlist(sapply(droppars, function(z) which(regexpr(z,  colnames(ss)) == 1))))
    if (length(wch) != 0) {
      ss <- ss[, -wch]
    }
    if (!is.null(keeppars)) {
      wch = unique(unlist(sapply(keeppars, function(z) which(regexpr(z,  colnames(ss)) == 1))))
      if (length(wch) != 0) {
        ss <- ss[, wch]
        }
      }
    }
  }
  else if (class(x) == "runjags") {
    ss <- runjags::combine.mcmc(x, collapse.chains = TRUE)
    ss <- as.matrix(ss)
    ss1 <- ss
    wch = unique(unlist(sapply(droppars, function(z) which(regexpr(z, colnames(ss)) == 1))))
    if (length(wch) != 0) {
      ss <- ss[, -wch]
    }
    if (!is.null(keeppars)) {
      wch = unique(unlist(sapply(keeppars, function(z) which(regexpr(z, colnames(ss)) == 1))))
      if (length(wch) != 0) {
        ss <- ss[, wch]
      }
    }
  }
  else {
    ss <- as.matrix(x)
    ss1 <- ss
    wch = unique(unlist(sapply(droppars, function(z) which(regexpr(z, 
                                                                   colnames(ss)) == 1))))
    if (length(wch) != 0) {
      ss <- ss[, -wch]
    }
    if (!is.null(keeppars)) {
      wch = unique(unlist(sapply(keeppars, function(z) which(regexpr(z, 
                                                                     colnames(ss)) == 1))))
      if (length(wch) != 0) {
        ss <- ss[, wch]
      }
    }
  }
  
  logdensityY = function(x) {
    d <- density(x, from = min(x), to = max(x), n = length(x), kernel = "g")
    lp = log(d$y / sum(d$y))
    lp
  }
  
  logdensityX = function(x) {
    d <- density(x, from = min(x), to = max(x), n = length(x), kernel = "g")
    zapsmall(round(d$x, 3) , 2)
  }
  
  logpost = apply(ss, 2, function(x) logdensityY(x))
  dens = apply(ss1, 2, function(x) logdensityX(x))
  round(dens[which.max(rowSums(logpost)) , ], 4)
}

