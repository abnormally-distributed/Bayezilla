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

#' Compute the kullback leibler divergence
#'
#' @param a the first distribution
#' @param b the second distribution
#' @param base exponent base, defaults to exp(1)
#' @export
#' @examples
#' kld()

kld = function (a, b, base = exp(1))
{
  if (!is.vector(a))
    a <- as.vector(a)
  if (!is.vector(b))
    b <- as.vector(b)
  n1 <- length(a)
  n2 <- length(b)
  if (!identical(n1, n2))
    stop("a and b must have the same length.")
  if (any(!is.finite(a)) || any(!is.finite(b)))
    stop("a and b must have finite values.")
  if (any(a <= 0))
    a <- exp(a)
  if (any(b <= 0))
    b <- exp(b)
  a[which(a < .Machine$double.xmin)] <- .Machine$double.xmin
  b[which(b < .Machine$double.xmin)] <- .Machine$double.xmin
  a <- a/sum(a)
  b <- b/sum(b)
  kld.a.b <- a * (log(a, base = base) - log(b, base = base))
  kld.b.a <- b * (log(b, base = base) - log(a, base = base))
  sum.kld.a.b <- sum(kld.a.b)
  sum.kld.b.a <- sum(kld.b.a)
  mean.kld <- (kld.a.b + kld.b.a)/2
  mean.sum.kld <- (sum.kld.a.b + sum.kld.b.a)/2
  out <- list(kld.a.b = kld.a.b, kld.b.a = kld.b.a,
              mean.kld = mean.kld, sum.kld.a.b = sum.kld.a.b, sum.kld.b.a = sum.kld.b.a,
              mean.sum.kld = mean.sum.kld, intrinsic.discrepancy = min(sum.kld.a.b,
                                                                       sum.kld.b.a))
  return(out)
}

#' Compute Jeffrey's Distance between two distributions
#'
#' Returns the Jeffrey's Distance
#'
#' @param a the first distribution
#' @param b the second distribution
#' @param base exponent base, defaults to exp(1)
#' @export
#' @examples
#' jeffreys_dist()
jeffreys_dist = function (a, b, base = exp(1)) {
  if (!is.vector(a))
    a <- as.vector(a)
  if (!is.vector(b))
    b <- as.vector(b)
  n1 <- length(a)
  n2 <- length(b)
  if (!identical(n1, n2))
    stop("a and b must have the same length.")
  if (any(!is.finite(a)) || any(!is.finite(b)))
    stop("a and b must have finite values.")
  if (any(a <= 0))
    a <- exp(a)
  if (any(b <= 0))
    b <- exp(b)
  a[which(a < .Machine$double.xmin)] <- .Machine$double.xmin
  b[which(b < .Machine$double.xmin)] <- .Machine$double.xmin
  a <- a/sum(a)
  b <- b/sum(b)
  kld.a.b <- a * (.5 * (log(a, base = base) - log(b, base = base)))
  kld.b.a <-  (b * (log(b, base = base) - log(a, base = base)))
  sum.kld.a.b <- sum(kld.a.b)
  sum.kld.b.a <- .5 * sum(kld.b.a)
  sum.kld.a.b + sum.kld.b.a
}

#' Compute Jaynes' Psi between two distributions
#'
#' Returns the Jaynes' Psi divergence
#'
#' @param a the first distribution
#' @param b the second distribution
#' @export
#' @examples
#' jaynes_psi()
jaynes_psi = function (a, b) {
  if (!is.vector(a))
    a <- as.vector(a)
  if (!is.vector(b))
    b <- as.vector(b)
  n1 <- length(a)
  n2 <- length(b)
  if (!identical(n1, n2))
    stop("a and b must have the same length.")
  if (any(!is.finite(a)) || any(!is.finite(b)))
    stop("a and b must have finite values.")
  if (any(a <= 0))
    a <- exp(a)
  if (any(b <= 0))
    b <- exp(b)
  a[which(a < .Machine$double.xmin)] <- .Machine$double.xmin
  b[which(b < .Machine$double.xmin)] <- .Machine$double.xmin
  a <- a/sum(a)
  b <- b/sum(b)
  dist = a  * (log10(a) - log10(b))
  dist = sum(dist)
  N = 10*n1
  N * dist
}

#' Binomial Confusion matrix for classification problems
#'
#' An alternative to the caret package's confusion matrix function.
#'
#' @param predictions the predicted classes.
#' @param actual the real outcome variable.
#' @export
#' @examples
#' confusion.matrix()
#'

confusion.matrix = function(predictions, actual) {
  
  table = table(predictions, actual)
  table = as.data.frame.matrix(table)
  
  A = table[1,1]
  B = table[1,2]
  C = table[2,1]
  D = table[2,2]
  
  MCC = ((A*D) - (C* B)) / sqrt((A+C)*(A+B)*(D+C)*(D+B))
  
  
  Accuracy = (A+D)/(A+B+C+D)
  NIR = max((A+C)/(A+B+C+D) , 1  - (A+C)/(A+B+C+D))
  Prevalence = (A+C)/(A+B+C+D)
  
  a.prior = 0.5
  b.prior = 0.5
  z = (A+D)
  N = (A+B+C+D)
  estimate.accuracy =  (z + a.prior) / (N + a.prior + b.prior)
  alpha1 = z + a.prior
  beta1 = N - z + b.prior
  lower = qbeta((1  - .90  )/ 2, alpha1, beta1)
  upper = qbeta(.90 +((1  - .90  )/ 2), alpha1, beta1)
  p = pbeta(NIR, alpha1, beta1, lower.tail = FALSE)
  prior.dens = dbeta(NIR, .5, .5)
  posterior.dens = dbeta(NIR, alpha1, beta1)
  BF = prior.dens/posterior.dens
  
  Sensitivity = A/(A+C)
  Specificity = D/(B+D)
  
  PPV = (Sensitivity * Prevalence)/((Sensitivity*Prevalence) + ((1-Specificity)*(1-Prevalence)))
  NPV = (Specificity * (1-Prevalence))/(((1-Sensitivity)*Prevalence) + ((Specificity)*(1-Prevalence)))
  FDR = 1-PPV
  FNDR = 1-NPV
  
  matrix0 = matrix(0, nrow=13, ncol=3)
  matrix0[1,] <- c("Classification", paste0("True", levels(actual)[1]), paste0("True", levels(actual)[2]))
  matrix0[c(2,3),1] <- c(levels(actual)[1], levels(actual)[2])
  matrix0[2,2] = A
  matrix0[2,3] = B
  matrix0[3,2] = C
  matrix0[3,3] = D
  matrix0[4,] = c("", "", "")
  matrix0[5,] = c("Accuracy", round(Accuracy, 4), "")
  matrix0[6,] = c("90% Credible Interval", round(lower,4), round(upper,4))
  matrix0[7,] = c("No Information Rate", round(NIR,3), " ")
  matrix0[8,] = c("Prob(Acc > NIR | Data)", round(p,4), paste0("(BF=", round(BF,4), ")"))
  matrix0[9,] = c("Sensitivity", round(Sensitivity,4), " ")
  matrix0[10,] = c("Specificity", round(Specificity,4), " ")
  matrix0[11,] = c("False Discovery Rate", round(FDR,4), " ")
  matrix0[12,] = c("False Nondiscovery Rate", round(FNDR,4), " ")
  matrix0[13,] = c(paste0("Matthew's Correlation"," ", "(" ,noquote("\u3d5"), ")"), round(MCC, 4), paste0(noquote(paste0("\u3d5", "\u00B2")), "=", round(MCC^2, 4)))
  knitr::kable(as.tibble(matrix0), col.names = c("","",""))
  
}

#' Conduct a ROPE test
#'
#' Conduct a test of practical equivalence by defining a region of practical equivalence (ROPE) around the null hypothesis H0.
#' This is most appropriate when you wish to test whether an effect is larger in magnitude than a cutoff that defines an effect
#' size that is practically unimportant. Hence, this is a test for practical significance and requires a well defined idea of what
#' practical is. 
#'
#' @param fit a stanfit or runjags object. 
#' @param param the name of the parameter in the stanfit or runjags object you want to plot.
#' @param ROPE a numeric vector of two numbers. Defaults to c(-.1, .1), assuming a Cohen's d... Adjust to what is reasonable!
#' @param interval.method One of "QI" (the default), "HDI"
#' @param cred.level the confidence level for interval.method. Defaults to  .90
#' @param plot plot the output. defaults to TRUE.
#' @param line.color color marking H0. "black" is default.
#' @param interval.color defaults to purple/violet "#8c00ffCC"
#' @param density.color defaults to blue "#2CBAFFCC"
#' @param rope.color defaults to red "#FF0000CC"
#' @export
#' @examples
#' ROPE()
#'
ROPE = function (fit, 
                      ROPE = c(-.1, .1), 
                      H0 = 0, 
                      plot = TRUE, 
                      param = NULL,
                      line.color = "black", 
                      interval.color = "#8c00ffCC", 
                      density.color = "#2CBAFFCC",
                      rope.color = "#FF0000CC",
                      interval.method = "QI",
                      cred.level = .90)
{
  
  if (is.null(param)) {stop("You must specify a single parameter. Either you have not provided a parameter name to the param argument, or have provided more than one.")}
  
  stan <- inherits(fit, "stanfit")
  if (stan == TRUE) {
    paramSampleVec <- as.matrix(fit)
    paramSampleVec = as.vector(paramSampleVec[,which(colnames(paramSampleVec) == param)])
  } 
  else if (class(fit) == "runjags") {
    paramSampleVec = as.matrix(combine.mcmc(fit, collapse.chains = TRUE, vars = param))
  }
  
  if (is.null(param)){
    param.label <- expression(theta)
  }
  else if (!is.null(param)){
    param.label <- noquote(param)
  }
  
  x = as.vector(as.matrix(paramSampleVec))
  
  if (interval.method == "QI") {
    interval <- Bayezilla:::cred_interval(as.vector(as.matrix(x)), 
                                          cred.level = cred.level, method = "QI")
    names(interval) <- c("lower", "upper")
  }
  else if (interval.method == "HDI") {
    interval <- Bayezilla:::cred_interval(as.vector(as.matrix(x)),
                                          cred.level = cred.level, method = "HDI")
    names(interval) <- c("lower", "upper")
  }
  
  if (length(ROPE) != 2)
    stop("Argument ROPE needs to be a vector of length two.",
         call. = F)
  if (ROPE[1] > ROPE[2]) {
    tmp <- ROPE[2]
    ROPE[2] <- ROPE[1]
    ROPE[1] <- tmp
  }
  
  original.x = x
  x <- sort(x)
  x.rope <- dplyr::between(x, interval[1], interval[2])
  x <- x[which(x.rope == TRUE)]
  r <- dplyr::between(x, ROPE[1], ROPE[2])
  rope.pct = sum(r)/length(x)
  rope.pct = data.frame(rope = rope.pct)
  interval = as.data.frame(t(interval))
  rope.pct$outside.rope = 1 - rope.pct$rope
  colnames(rope.pct) <- c("inside", "outside")
  .neff <- length(x)
  result <- dplyr::case_when(interval$lower > ROPE[2] ~ "reject null",
                             interval$upper < ROPE[1] ~ "reject null", interval$lower >=
                               ROPE[1] & interval$upper <= ROPE[2] ~ "accept null",
                             TRUE ~ "undecided")
  
  result = cbind.data.frame(rope.pct, result, stringsAsFactors = FALSE)
  names(result) = c("Inside ROPE", "Outside ROPE", "Decision")
  
  if (plot == TRUE) {
    xdf = data.frame(posterior.samples = as.vector(original.x))
    ROPE = data.frame(t(ROPE))
    plot <- xdf %>% ggplot2::ggplot(aes(x = posterior.samples,
                                        y = ..density..)) + 
      stat_density(fill = density.color, alpha = 1, adjust = 3)
    
    if (result$`Inside ROPE` == 0) {
      plot <- plot + theme(legend.title = element_blank()) +
        ggplot2::geom_vline(xintercept = H0, color = line.color,
                            linetype = "dotted", size = 1) + 
        geom_segment(aes(y = 0, yend = 0, xend = upper, x = lower), 
                     data = interval, size = 2, alpha = .90, color = interval.color) + 
        labs(x = param.label)
    }
    else {
      d <- ggplot2::ggplot_build(plot)$data[[1]]
      plot <- plot + 
        geom_area(data = subset(d, x > ROPE$X1 & x < ROPE$X2), fill = rope.color, aes(x = x, y = y), alpha = 0.85) + 
        theme(legend.title = element_blank()) +
        ggplot2::geom_vline(xintercept = H0, color = line.color, linetype = "dotted", size = 1) + 
        geom_segment(aes(y = 0, yend = 0, xend = upper, x = lower), data = interval,size = 2, alpha = .90, color = interval.color) +
        labs(x = param.label)
    }
    print(result)
    plot
  }
  else {
    return(result)
  }
}
