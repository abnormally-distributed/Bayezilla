#' Scatterplots with Marginal Histograms
#'
#' @param x the x-axis variable or a formula
#' @param y the y-axis variable or a data frame if x is a formula
#' @param xlab the x axis label
#' @param ylab the y axis label
#' @param col color scheme. one of "blue" (the default), "red", "green", or "purple"
#' @param x.breaks breaks in x axis histogram. defaults to 15.
#' @param y.breaks breaks in y axis histogram. defaults to 15. 
#' @param smooth should a smoother line be added. Defaults to FALSE. 
#' @param regline should a regression line be added? Defaults to TRUE, with both robust and ordinary least squares lines being plotted.
#'
#' @return
#' a plot
#' @export
#'
#' @examples
#' scatterPlot(iris$Sepal.Width, iris$Petal.Length, xlab = "Sepal Width", ylab = "Petal Length", col = "purple")
scatPlot = function(x, y, xlab="x", ylab="y", col = "blues", hist = TRUE, x.breaks = "dhist", y.breaks = "dhist", smooth = FALSE, regline = TRUE){
  
  old.par <- par(no.readonly = TRUE) # save default, for resetting... 
  on.exit(par(old.par))     #and when we quit the function, restore to original values
  
  if (class(x) == "formula"){
    
    if (class(y) != "data.frame"){
      stop("When providing a formula to x, y must be a data frame")
    } else if (class(y) == "data.frame"){
      mf = model.frame(x, y)
      y = as.vector(mf[,2])
      x = as.vector(mf[,1])
      data = data.frame(y = y, x = x)
    }
    
  } else {
    
    data = data.frame(y = y, x = x)
    
  }
  
  #light, dark, dark, light, smoothline 
  if (col == "blues" || col == "blue"){
    ColorScheme = c("#6495edCC", "#0d2f6dCC", "#0a2556CC", "#6dabff5E", "#1f7dffA1", "#022e7dCC")
  }
  if (col == "purples" || col == "purple"){
    ColorScheme = c("#d4aad4CC" , "#400040CC", "#270027CC", "#e035e05E", "#8347b3A1", "#33253dCC")
  }
  if (col == "reds" || col == "red"){
    ColorScheme = c("#e87a7aCC", "#8b1a1aCC", "#910202CC", "#f44b755E", "#ec0518A1", "#4f0d0dCC")
  }
  if (col == "greens" || col == "green"){
    ColorScheme= c("#66cc66BF", "#194d33CC", "#133a26CC", "#61db5e5E", "#42d43eA1", "#163324CC")
  }
  
  tukey.wts = function (r, t = 4.685){
    return((1 - pmin(1, abs(r/t))^2)^2)
  }
  
  xrange <- range(x)
  yrange <- range(y)
  
  if (hist == TRUE){
    zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
    layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
    xhist = hist(x, plot=FALSE, breaks = x.breaks)
    yhist = hist(y, plot=FALSE, breaks = y.breaks)
    top = max(c(xhist$counts, yhist$counts))
    par(mar=c(4, 4,.5,.5), oma = c(1.5, 1.5, 0, 0))
    plot(x, y, xlab = xlab, ylab = ylab, xlim=xrange, ylim=yrange, xaxt = "n", yaxt = "n", lwd = 1.225, cex = 1.25, lty = 1, pch = 21, bty="n", family = "serif", col = ColorScheme[2], bg = ColorScheme[1], cex.lab=1.25)
    if(smooth == TRUE){
      lines(smooth.spline(x, y), col= ColorScheme[5], lwd = 3)
    } 
    if(regline == TRUE){
      resids = model.frame(y ~ x, data = data)[,1] - as.vector(lmSolve(y ~ x, data = data) %*% t(model.matrix(y ~ x, data = data)))
      w = tukey.wts((resids - mean(resids)) / sd(resids))
      robust = round(coef(MASS::rlm(y ~ x, data = data, scale.est = "Huber", init = "lts", method = "MM", psi = MASS::psi.hampel, w = w, acc = 1e-3, maxit = 500)), 3)
      ordinary = coef(lm(y ~ x))
      abline(robust[1], robust[2], col = ColorScheme[6], lwd = 3)
      abline(ordinary[1], ordinary[2], col = "#1b1e24A1", lwd = 2.5, lty = 2)
      legend("topright", legend=c("Robust", "OLS"),
             col=c(ColorScheme[6], "#1b1e24A1"), lty=c(1,3), lwd = 2, cex = 1, inset=.05)
    }
    axis(1, col = NA, tck = 0, family = "serif", cex = 1.5)
    axis(2, col = NA, tck = 0, family = "serif", cex = 1.5)
    par(mar=c(0,4,3,2))
    barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0, border = ColorScheme[3], col = ColorScheme[4])
    par(mar=c(4,0,2,3))
    barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE, border = ColorScheme[3], col = ColorScheme[4])
    par(oma=c(4,4,0,0))
    
  } 
  
  if (hist == FALSE){
    plot(x, y, xlab = xlab, ylab = ylab, xlim=xrange, ylim=yrange, xaxt = "n", yaxt = "n", lwd = 1.225, cex = 1.25, lty = 1, pch = 21, bty="n", family = "serif", col = ColorScheme[2], bg = ColorScheme[1], cex.lab=1.25)
    if(smooth == TRUE){
      lines(smooth.spline(x, y), col= ColorScheme[5], lwd = 3)
    } 
    if(regline == TRUE){
      resids = model.frame(y ~ x, data = data)[,1] - as.vector(lmSolve(y ~ x, data = data) %*% t(model.matrix(y ~ x, data = data)))
      w = tukey.wts((resids - mean(resids)) / sd(resids))
      robust = round(coef(MASS::rlm(y ~ x, data = data, scale.est = "Huber", init = "lts", method = "MM", psi = MASS::psi.hampel, w = w, acc = 1e-3, maxit = 500)), 3)
      
      ordinary = coef(lm(y ~ x))
      abline(robust[1], robust[2], col = ColorScheme[6], lwd = 3)
      abline(ordinary[1], ordinary[2], col = "#1b1e24A1", lwd = 2, lty = 2)
      legend("topright", legend=c("Robust", "OLS"),
             col=c(ColorScheme[6], "#1b1e24A1"), lty=c(1,3), lwd = 2, cex = 1, inset=.05)
    }
    
    axis(1, col = NA, tck = 0, family = "serif", cex = 1.5)
    axis(2, col = NA, tck = 0, family = "serif", cex = 1.5)
  }
}

#' Visualize your data with a scatterplot matrix
#'
#' @description This is a function to help quickly visualize the relationships among your data.
#' Like the other pre-themed plotting functions in this package, this utilizes a
#' minimalist Tufte inspired aesthetic.
#'
#' @param x The data frame
#' @param smooth Should the smoother line be enabled?
#' @param digits Number of significant digits. defaults to 2.
#' @param method "pearson", "spearman", or "kendall"
#' @param pch The symbol used for plotting 
#' @param lm TRUE or FALSE (Defaults to FALSE)
#' @param cor TRUE or FALSE (Defaults to TRUE)
#' @param jitter TRUE or FALSE (Defaults to FALSE)
#' @param amount Amount of jittering
#' @param hist.col The histogram color. Defaults to "#0054f9A6".
#' @param points.col The color of the data points. Defaults to "#35a8f59E".
#' @param smooth.col The color of the smooth line. Defaults to "#663399BF"
#' @param breaks Breaks in histogram. Defaults to "FD" method.
#' @param cex.cor Size of points
#' @param cex.num The size of the numbers in the upper right triangle of the scatterplot. Defaults to 1.5
#' @param smoother If TRUE, then smooth.scatter the data points – slow but pretty with lots of subjects
#' @param ... other arguments
#'
#' @return
#' a plot
#' @export
#'
#' @examples
#' scatMat(iris)
scatMat <-
  function(x,
           smooth = TRUE,
           digits = 2,
           method = "pearson",
           pch = 19,
           lm = FALSE,
           cor = TRUE,
           jitter = FALSE,
           amount = 2,
           hist.col = "#0054f9A6",
           points.col = "#35a8f59E",
           smooth.col = "#663399BF",
           breaks = "FD",
           cex.cor = 1,
           cex.num = 1.5,
           smoother = FALSE,
           ...) {
    density <- FALSE
    scale <- TRUE
    wt <- NULL
    show.points <- TRUE
    rug <- FALSE
    
    "panel.hist.density" <-
      function(x, ...) {
        usr <- par("usr")
        on.exit(par(usr))
        par(usr = c(usr[1], usr[2], 0, 1.5))
        tax <- table(x)
        if (length(tax) < 11) {
          breaks <- as.numeric(names(tax))
          y <- tax / max(tax)
          interbreak <- min(diff(breaks)) * (length(tax) - 1) / 41
          rect(breaks - interbreak, 0, breaks + interbreak, y, col = hist.col)
        } else {
          h <- hist(x, breaks = breaks, plot = FALSE)
          breaks <- h$breaks
          nB <- length(breaks)
          y <- h$counts
          y <- y / max(y)
          rect(breaks[-nB], 0, breaks[-1], y, col = hist.col)
        }
        if (density) {
          tryd <-
            try(d <- density(x,
                             na.rm = TRUE,
                             bw = "nrd",
                             adjust = 1.5
            ), silent = TRUE)
          if (class(tryd) != "try-error") {
            d$y <- d$y / max(d$y)
            lines(d)
          }
        }
        if (rug) {
          rug(x)
        }
      }
    
    
    "panel.cor" <-
      function(x, y, prefix = "", ...) {
        usr <- par("usr")
        on.exit(par(usr))
        par(usr = c(0, 1, 0, 1))
        if (is.null(wt)) {
          r <- cor(x, y, use = "pairwise", method = method)
        }
        else {
          r <- cor.wt(data.frame(x, y), w = wt[, c(1:2)])$r[1, 2]
        }
        txt <-
          format(c(round(r, digits), 0.123456789), digits = digits)[1]
        txt <- paste(prefix, txt, sep = "")
        cex <- cex.num
        text(0.5, 0.5, txt, cex = cex)
      }
    
    "panel.smoother" <-
      function(x,
               y,
               pch = par("pch"),
               col.smooth = smooth.col,
               span = 2 / 3,
               iter = 3,
               lwd = 3,
               ...) {
        xm <- mean(x, na.rm = TRUE)
        ym <- mean(y, na.rm = TRUE)
        xs <- sd(x, na.rm = TRUE)
        ys <- sd(y, na.rm = TRUE)
        r <- cor(x, y, use = "pairwise", method = method)
        
        if (jitter) {
          x <- jitter(x, factor = amount)
          y <- jitter(y, factor = amount)
        }
        
        if (smoother) {
          if (show.points) {
            points(x, y, pch = pch, col = points.col, ...)
          }
          axis(1, col = NA, tck = 0)
          axis(2, col = NA, tck = 0)
        } else {
          if (show.points) {
            points(x, y, pch = pch, col = points.col, ...)
          }
        }
        ok <- is.finite(x) & is.finite(y)
        if (any(ok)) {
          if (smooth) {
            
            if (length(unique(x)) < 10 || length(unique(y)) < 10){
              smoothed = stats::lowess(x[ok], y[ok], f = span, iter = iter)
            } else {
              smoothed = smooth.spline(x[ok], y[ok])
            }
            
            lines(smoothed ,
                  col = smooth.col,
                  lwd = 3
            )
          }
        }
      }
    
    "panel.lm" <-
      function(x,
               y,
               pch = par("pch"),
               col.lm = "red",
               ...) {
        ymin <- min(y)
        ymax <- max(y)
        xmin <- min(x)
        xmax <- max(x)
        ylim <- c(min(ymin, xmin), max(ymax, xmax))
        xlim <- ylim
        if (jitter) {
          x <- jitter(x, factor = amount)
          y <- jitter(y, factor = amount)
        }
        if (smoother) {
          smoothScatter(x, y, add = TRUE, nrpoints = 0)
        } else {
          if (show.points) {
            points(
              x,
              y,
              pch = pch,
              ylim = ylim,
              xlim = xlim,
              col = points.col,
              ...
            )
          }
        }
        ok <- is.finite(x) & is.finite(y)
        if (any(ok)) {
          lml <- lm(y ~ x)
        }
      }
    
    
    #######
    
    old.par <-
      par(no.readonly = TRUE) # save default, for resetting...
    on.exit(par(old.par)) # and when we quit the function, restore to original values
    
    par(
      xaxt = "n",
      yaxt = "n",
      bty = "null",
      family = "serif"
    )
    
    if (missing(cex.cor)) {
      cex.cor <-
        1
    } # this allows us to scale the points separately from the correlations
    
    for (i in 1:ncol(x)) {
      # treat character data as numeric
      if (is.character(x[[i]])) {
        x[[i]] <- as.numeric(as.factor(x[[i]]))
        colnames(x)[i] <- paste(colnames(x)[i], "*", sep = "")
      }
    }
    n.obs <- nrow(x)
    
    if (!lm) {
      # the basic default is here
      if (cor) {
        pairs(
          x,
          diag.panel = panel.hist.density,
          upper.panel = panel.cor,
          lower.panel = panel.smoother,
          pch = pch,
          cex = .90,
          bg = points.col,
          ...
        )
        axis(1, col = NA, tck = 0)
        axis(2, col = NA, tck = 0)
      }
      else {
        pairs(
          x,
          diag.panel = panel.hist.density,
          upper.panel = panel.smoother,
          lower.panel = panel.smoother,
          pch = pch,
          bg = points.col,
          cex = .90,
          ...
        )
        axis(1, col = NA, tck = 0)
        axis(2, col = NA, tck = 0)
      }
    } else {
      # lm is TRUE
      if (!cor) {
        # this case does not show the correlations, but rather shows the regression lines above and below the diagonal
        
        pairs(
          x,
          diag.panel = panel.hist.density,
          upper.panel = panel.lm,
          lower.panel = panel.lm,
          pch = pch,
          bg = points.col,
          cex = .90,
          ...
        )
        axis(1, col = NA, tck = 0)
        axis(2, col = NA, tck = 0)
      } else {
        # the normal case is to show the regressions below and the rs above
        
        pairs(
          x,
          diag.panel = panel.hist.density,
          upper.panel = panel.cor,
          lower.panel = panel.lm,
          pch = pch,
          bg = points.col,
          cex = .90,
          ...
        )
        axis(1, col = NA, tck = 0)
        axis(2, col = NA, tck = 0)
      }
    }
  }