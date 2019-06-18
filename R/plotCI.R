#' Plot Credible Intervals of a Bayesian Model 
#'
#' @param fit the model fit
#' @param keeppars variables to keep
#' @param droppars  variables to drop. defaults to c("ySim", "log_lik", "lp__")
#' @param cred.level Confidence level. Defaults to .90.
#' @param method credible interval methods. Defaults to QI (quantile intervals).
#' @param H0 the null hypothesis. defaults to 0.
#' @param SEbars Defaults to FALSE. If TRUE, superimposes 68.2 percent intervals over wider ones.
#' @param ROPE If you would like, ROPE limits added to the plot.
#' @param point.shape default to 18 (diamond)
#' @param rope.color rope color
#' @param SEbar.color standard error bar color
#' @param point.color point color
#' @param bar.color bar color
#' @param point.size point size
#' @param bar.size bar size
#' @param ... other arguments to pass to post_summary, which is called internally
#' @export
#' @examples
#' plotCI()
#'
plotCI <- function(fit, keeppars = NULL, droppars = c("ySim", "log_lik", "lp__"), 
                            estimate = "mean",
                            cred.level = 0.90, 
                            method = "QI",
                            ROPE = NULL, 
                            rope.color = "black", 
                            SEbars = FALSE,
                            SEbar.color = "darkslategray3", 
                            point.color = "firebrick2", 
                            point.size = 2.75, 
                            point.shape = 18,
                            bar.color = "skyblue4", 
                            bar.size = 0.8, 
                            H0 = 0, 
                            ...) {
  
  stan <- inherits(fit, "stanfit")
  if (stan == TRUE) {
    codaObject <- as.matrix(fit)
    wch = unique(unlist(sapply(droppars, function(z) which(regexpr(z, colnames(codaObject)) == 1))))
    if (length(wch) != 0){
      codaObject <- codaObject[,-wch]
    }
    if (!is.null(keeppars)) {
      wch = unique(unlist(sapply(keeppars, function(z) which(regexpr(z, colnames(codaObject)) == 1))))
      if (length(wch) != 0){
        codaObject <- codaObject[,wch]
      }
    }
  }
  else if (class(fit) == "runjags"){
    codaObject <- runjags::combine.mcmc(fit, collapse.chains = TRUE)
    codaObject <- as.matrix(codaObject)
    wch = unique(unlist(sapply(droppars, function(z) which(regexpr(z, colnames(codaObject)) == 1))))
    if (length(wch) != 0){
      codaObject <- codaObject[,-wch]
    }
    if (!is.null(keeppars)) {
      wch = unique(unlist(sapply(keeppars, function(z) which(regexpr(z, colnames(codaObject)) == 1))))
      if (length(wch) != 0){
        codaObject <- codaObject[,wch]
      }
    }
  }
  else {
    codaObject <- as.matrix(fit)
    wch = unique(unlist(sapply(droppars, function(z) which(regexpr(z, colnames(codaObject)) == 1))))
    if (length(wch) != 0){
      codaObject <- codaObject[,-wch]
    }
    if (!is.null(keeppars)) {
      wch = unique(unlist(sapply(keeppars, function(z) which(regexpr(z, colnames(codaObject)) == 1))))
      if (length(wch) != 0){
        codaObject <- codaObject[,wch]
      }
    }
  }
  
  
  summary <- apply(codaObject, 2, function(x) cred_interval(object = x, cred.level = cred.level, method = method))
  rownames(summary) = c("cred.low", "cred.high")
  summary = as.data.frame(t(summary))

  if (estimate == "mean"){
    estimate <- colMeans(codaObject)
  }
  else {
    estimate <- apply(codaObject, 2, median)
  }
  
  sum <- cbind.data.frame(
    terms = factor(colnames(codaObject), levels = colnames(codaObject)),
    estimates = estimate, lower = summary$cred.low,
    upper = summary$cred.high)
  
  if (is.null(ROPE) == FALSE) {
    plot <- sum %>%
      ggplot(aes(y = forcats::fct_rev(terms), x = estimates, xmin = lower, xmax = upper)) + 
      labs(y = "term", x = "\u03B8") +
      geom_vline(xintercept = ROPE[1], linetype = "dotted", size = 0.7, color = rope.color) + 
      geom_vline(xintercept = ROPE[2], linetype = "dotted", size = 0.7, color = rope.color) 
  } else if (is.null(ROPE) == TRUE){
    plot <- sum %>%
      ggplot(aes(
        y = forcats::fct_rev(terms), x = estimates, xmin = lower, xmax = upper)) + 
      labs(y = "term", x = "\u03B8") +
      geom_vline(xintercept = H0, linetype = "dotted", size = 0.7, color = rope.color)
  }
  
  if (SEbars == TRUE) {
    
    SE <- apply(codaObject, 2, function(x) cred_interval(object = x, cred.level = 0.682, method = "QI"))
    rownames(SE) = c("SElower", "SEupper")
    SE = as.data.frame(t(SE))
    
    sum$SElower = SE$SElower
    sum$SEupper = SE$SEupper
    
    plot <- plot + 
      ggplot2::geom_errorbarh(color = bar.color,
      height = 0.29, size = bar.size, alpha = 1) +
      ggstance::geom_linerangeh(xmin = sum$SElower, xmax = sum$SEupper, color = SEbar.color, size = bar.size * 1.45, alpha = .80, position = position_nudge(y = 0, x = 0)) + 
      geom_point(size = point.size, color = point.color, shape = point.shape)
  }
  
  else {
    plot <- plot + 
      ggplot2::geom_errorbarh(color = bar.color, height = 0.29, size = bar.size, alpha = 1) + 
      geom_point(size = point.size, color = point.color, shape = point.shape)
  }
  
  plot(plot)
}
