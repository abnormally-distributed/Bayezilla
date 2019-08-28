#' Bivariate kernel density plot
#'
#' Calculates the joint density of two variables and then plots them using persp.
#'
#' @param x the data from which the estimate is to be computed. For the default method a numeric vector: long vectors are not supported.
#' @param y the data from which the estimate is to be computed. For the default method a numeric vector: long vectors are not supported.
#' @param plot defaults to TRUE. Otherwise, returns a list of x, y, and z.
#' @param xlim x-, y- and z-limits. These should be chosen to cover the range of values of the surface
#' @param ylim x-, y- and z-limits. These should be chosen to cover the range of values of the surface
#' @param zlim x-, y- and z-limits. These should be chosen to cover the range of values of the surface
#' @param xlab titles for the axes. N.B. These must be character strings; expressions are not accepted. Numbers will be coerced to character strings.
#' @param ylab titles for the axes. N.B. These must be character strings; expressions are not accepted. Numbers will be coerced to character strings.
#' @param zlab titles for the axes. N.B. These must be character strings; expressions are not accepted. Numbers will be coerced to character strings.
#' @param main main and sub title, as for title.
#' @param sub main and sub title, as for title.
#' @param theta angles defining the viewing direction. theta gives the azimuthal direction and phi the colatitude.
#' @param phi angles defining the viewing direction. theta gives the azimuthal direction and phi the colatitude.
#' @param r the distance of the eyepoint from the centre of the plotting box.
#' @param d a value which can be used to vary the strength of the perspective transformation. Values of d greater than 1 will lessen the perspective effect and values less and 1 will exaggerate it.
#' @param scale before viewing the x, y and z coordinates of the points defining the surface are transformed to the interval [0,1]. If scale is TRUE the x, y and z coordinates are transformed separately. If scale is FALSE the coordinates are scaled so that aspect ratios are retained. This is useful for rendering things like DEM information.
#' @param expand a expansion factor applied to the z coordinates. Often used with 0 < expand < 1 to shrink the plotting box in the z direction.
#' @param col the color(s) of the surface facets. Transparent colours are ignored. This is recycled to the (nx-1)(ny-1) facets.
#' @param viridis.opt A character string indicating the viridis palette variant to use. Four options are available: "magma" (or "A"), "inferno" (or "B"), "plasma" (or "C"), "viridis" (or "D", the default option) and "cividis" (or "E").
#' @param border the color of the line drawn around the surface facets. The default of NA will disable the drawing of borders: this is sometimes useful when the surface is shaded. NULL corresponds to par("fg").
#' @param ltheta if finite values are specified for ltheta and lphi, the surface is shaded as though it was being illuminated from the direction specified by azimuth ltheta and colatitude lphi.
#' @param lphi if finite values are specified for ltheta and lphi, the surface is shaded as though it was being illuminated from the direction specified by azimuth ltheta and colatitude lphi.
#' @param shade the shade at a surface facet is computed as ((1+d)/2)^shade, where d is the dot product of a unit vector normal to the facet and a unit vector in the direction of a light source. Values of shade close to one yield shading similar to a point light source model and values close to zero produce no shading. Values in the range 0.5 to 0.75 provide an approximation to daylight illumination.
#' @param box 	should the bounding box for the surface be displayed. The default is TRUE.
#' @param axes should ticks and labels be added to the box. The default is TRUE. If box is FALSE then no ticks or labels are drawn.
#' @param nticks character: "simple" draws just an arrow parallel to the axis to indicate direction of increase; "detailed" draws normal ticks as per 2D plots.
#' @param ticktype the (approximate) number of tick marks to draw on the axes. Has no effect if ticktype is "simple".
#' @param plotly whether or not to use plotly for the plotting. If TRUE, most of the arguments except
#' col, xlab, ylab, and zlab are ignored.
#' 
#' @return either a list or a plot
#' @export
#'
#' @examples
#' epdf2d(wines$Alcohol, wines$Malic.acid)
epdf2d = function(x, y, xlim = range(x), ylim = range(y),
                  zlim = range(z, na.rm = TRUE),
                  xlab = "x", ylab = "y", zlab = "Joint PDF",
                  main = NULL, sub = NULL,
                  theta = 40, phi = 10, r = exp(1), d = 1,
                  scale = TRUE, expand = 1,
                  col = NULL, viridis.opt = "D" , border = NA, ltheta = -135, lphi = 15,
                  shade = NA, box = TRUE, axes = TRUE, nticks = 4,
                  ticktype = "detailed", plot = TRUE, plotly = FALSE){
  

  dx = density(x, n = 43, kernel = "o", bw = "nrd0")$y
  dy = density(y, n = 43, kernel = "o", bw = "nrd0")$y
  
  wch = c(max(dx), max(dy))[which.max(c(max(dx), max(dy)))]
  dx = dx * solve(wch, 100)
  dy = dy * solve(wch, 100)
  if (any(dx < 0 )){
    dx = dx + abs(min(dx))
  }
  if (any(dy < 0 )){
    dy = dy + abs(min(dy))
  }
  dx = ((dx-min(dx))/(max(dx)-min(dx))) * max(dnorm(x, mean(x), sd(x)))
  dy = ((dy-min(dy))/(max(dy)-min(dy))) * max(dnorm(y, mean(y), sd(y)))
  
  z = tcrossprod(dx, dy)
  x = seq(min(x), max(x), len = 43)
  y = seq(min(y), max(y), len = 43)
  X = list(x = sort(x), y = sort(y), z = z)
  
  
  if (plot){
    
    if (plotly){
      
      if (is.null(col)){
        colors = substr(viridis::viridis(1849, alpha = 1, begin = .084, end = 1, option = viridis.opt), 1, 7)
      } else {
        if (is.function(col)) {
          color_function <- col
          colors = color_function(1849)
        } else {
          stop("Please provide a color palette function")
        }
      }
      
      require(plotly)
      
      xlab <- list(
        title = xlab ,
        titlefont = par()$family
      )
      ylab <- list(
        title = ylab,
        titlefont = par()$family
      )
      zlab <- list(
        title = zlab,
        titlefont = par()$family
      )
      
      plot_ly(z = ~X$z, x = ~ X$x, y = ~ X$y, colors = colors, showscale=FALSE) %>% 
        add_surface() %>%
        layout(scene = list(xaxis = xlab, yaxis = ylab, zaxis = zlab), showlegend = FALSE)
      
    } else{
      
      if (is.null(col)){
        colors = substr(viridis::viridis(1849, alpha = 1, begin = 0, end = 1, option = viridis.opt), 1, 7)
        z.facet.center <- (z[-1, -1] + z[-1, -43] + z[-43, -1] + z[-43, -43])/4
        z.facet.range = cut(z.facet.center, 1849)
        col = colors[z.facet.range]
      } else {
        if (is.function(col)) {
          color_function <- col
          colors = color_function(1849)
          z.facet.center <- (z[-1, -1] + z[-1, -43] + z[-43, -1] + z[-43, -43])/4
          z.facet.range = cut(z.facet.center, 1849)
          col = colors[z.facet.range]
        } else {
          stop("Please provide a color palette function")
        }
      }
      
      
      persp(X, xlim = xlim, ylim = ylim, zlim = zlim, xlab = paste0("\n", xlab), ylab = paste0(ylab, "\n\n"),
            zlab = paste0("\n", zlab), main = main, sub = sub, theta = theta, phi = phi, r = r, 
            d= d, scale = scale, expand = expand, col = col, border = border, 
            ltheta = ltheta, lphi = lphi, shade = shade, box = box, axes = axes, 
            nticks = nticks, ticktype = ticktype)
    }
  } else if (!plot){
    return(X)
  }
}

#' Bivariate empirical cumulative density function plot 
#'
#' Calculates the joint empirical cumulative probability density function of two variables and then plots them using persp.
#'
#' @param x the data from which the estimate is to be computed. For the default method a numeric vector: long vectors are not supported.
#' @param y the data from which the estimate is to be computed. For the default method a numeric vector: long vectors are not supported.
#' @param plot defaults to TRUE. Otherwise, returns a list of x, y, and z.
#' @param xlim x-, y- and z-limits. These should be chosen to cover the range of values of the surface
#' @param ylim x-, y- and z-limits. These should be chosen to cover the range of values of the surface
#' @param zlim x-, y- and z-limits. These should be chosen to cover the range of values of the surface
#' @param xlab titles for the axes. N.B. These must be character strings; expressions are not accepted. Numbers will be coerced to character strings.
#' @param ylab titles for the axes. N.B. These must be character strings; expressions are not accepted. Numbers will be coerced to character strings.
#' @param zlab titles for the axes. N.B. These must be character strings; expressions are not accepted. Numbers will be coerced to character strings.
#' @param main main and sub title, as for title.
#' @param sub main and sub title, as for title.
#' @param theta angles defining the viewing direction. theta gives the azimuthal direction and phi the colatitude.
#' @param phi angles defining the viewing direction. theta gives the azimuthal direction and phi the colatitude.
#' @param r the distance of the eyepoint from the centre of the plotting box.
#' @param d a value which can be used to vary the strength of the perspective transformation. Values of d greater than 1 will lessen the perspective effect and values less and 1 will exaggerate it.
#' @param scale before viewing the x, y and z coordinates of the points defining the surface are transformed to the interval [0,1]. If scale is TRUE the x, y and z coordinates are transformed separately. If scale is FALSE the coordinates are scaled so that aspect ratios are retained. This is useful for rendering things like DEM information.
#' @param expand a expansion factor applied to the z coordinates. Often used with 0 < expand < 1 to shrink the plotting box in the z direction.
#' @param col the color(s) of the surface facets. Transparent colours are ignored. This is recycled to the (nx-1)(ny-1) facets.
#' @param viridis.opt A character string indicating the viridis palette variant to use. Four options are available: "magma" (or "A"), "inferno" (or "B"), "plasma" (or "C"), "viridis" (or "D", the default option) and "cividis" (or "E").
#' @param border the color of the line drawn around the surface facets. The default of NA will disable the drawing of borders: this is sometimes useful when the surface is shaded. NULL corresponds to par("fg").
#' @param ltheta if finite values are specified for ltheta and lphi, the surface is shaded as though it was being illuminated from the direction specified by azimuth ltheta and colatitude lphi.
#' @param lphi if finite values are specified for ltheta and lphi, the surface is shaded as though it was being illuminated from the direction specified by azimuth ltheta and colatitude lphi.
#' @param shade the shade at a surface facet is computed as ((1+d)/2)^shade, where d is the dot product of a unit vector normal to the facet and a unit vector in the direction of a light source. Values of shade close to one yield shading similar to a point light source model and values close to zero produce no shading. Values in the range 0.5 to 0.75 provide an approximation to daylight illumination.
#' @param box 	should the bounding box for the surface be displayed. The default is TRUE.
#' @param axes should ticks and labels be added to the box. The default is TRUE. If box is FALSE then no ticks or labels are drawn.
#' @param nticks character: "simple" draws just an arrow parallel to the axis to indicate direction of increase; "detailed" draws normal ticks as per 2D plots.
#' @param ticktype the (approximate) number of tick marks to draw on the axes. Has no effect if ticktype is "simple".
#' @param plotly whether or not to use plotly for the plotting. If TRUE, most of the arguments except
#' col, xlab, ylab, and zlab are ignored.
#'
#' @return either a list or a plot
#' @export
#'
#' @examples
#' ecdf2d(wines$Alcohol, wines$Malic.acid)
ecdf2d = function(x, y, xlim = range(x), ylim = range(y),
                  zlim = range(z, na.rm = TRUE),
                  xlab =" x", ylab = "y", zlab = "Joint CDF",
                  main = NULL, sub = NULL,
                  theta = -20, phi = 15, r = exp(1), d=25, 
                  scale = TRUE, expand = 1,
                  col = NULL, viridis.opt = "D" , border = NA, ltheta = -135, lphi = 0,
                  shade = NA, box = TRUE, axes = TRUE, nticks = 4,
                  ticktype = "detailed", plot = TRUE, plotly = FALSE){
  
  dx = ecdf(x)
  dy = ecdf(y)
  
  xfun = splinefun(jitter(x, amount = 0.0001), dx(x), method = "monoH.FC")
  yfun = splinefun(jitter(y, amount = 0.0001), dy(y), method = "monoH.FC")
  
  x = seq(min(x), max(x), len = 43)
  y = seq(min(y), max(y), len = 43)
  
  dx = xfun(x)
  dy = yfun(y)
  
  if (any(dx < 0 )){
    dx = dx + abs(min(dx))
  }
  if (any(dy < 0 )){
    dy = dy + abs(min(dy))
  }
  dx = (dx-min(dx))/(max(dx)-min(dx))
  dy = (dy-min(dy))/(max(dy)-min(dy))
  
  z = tcrossprod(dx, dy)
  X = list(x = sort(x), y = sort(y), z = z)
  

  if (plot){
    
    if (plotly){
      
      if (is.null(col)){
        colors = substr(viridis::viridis(1849, alpha = 1, begin = 0, end = 1, option = viridis.opt), 1, 7)
      } else {
        if (is.function(col)) {
          color_function <- col
          colors = color_function(1849)
          z.facet.center <- (z[-1, -1] + z[-1, -43] + z[-43, -1] + z[-43, -43])/4
          z.facet.range = cut(z.facet.center, 1849)
          colors = colors[z.facet.range]
        } else {
          stop("Please provide a color palette function")
        }
      }
      
      require(plotly)
      
      xlab <- list(
        title = xlab ,
        titlefont = par()$family
      )
      ylab <- list(
        title = ylab,
        titlefont = par()$family
      )
      zlab <- list(
        title = zlab,
        titlefont = par()$family
      )
      
      plot_ly(z = ~X$z, x = ~ X$x, y = ~ X$y, colors = colors, showscale=FALSE) %>% 
        add_surface() %>%
        layout(scene = list(xaxis = xlab, yaxis = ylab, zaxis = zlab), showlegend = FALSE)
      
    } else{
      
      if (is.null(col)){
        colors = substr(viridis::viridis(1849, alpha = 1, begin = .084, end = 1, option = viridis.opt), 1, 7)
        z.facet.center <- (z[-1, -1] + z[-1, -43] + z[-43, -1] + z[-43, -43])/4
        z.facet.range = cut(z.facet.center, 1849)
        col = colors[z.facet.range]
      } else {
        if (is.function(col)) {
          color_function <- col
          colors = color_function(1849)
          z.facet.center <- (z[-1, -1] + z[-1, -43] + z[-43, -1] + z[-43, -43])/4
          z.facet.range = cut(z.facet.center, 1849)
          col = colors[z.facet.range]
        } else {
          stop("Please provide a color palette function")
        }
      }

    persp(X, xlim = xlim, ylim = ylim, zlim = zlim, xlab = paste0("\n", xlab), ylab = paste0(ylab, "\n\n") ,
          zlab = paste0("\n", zlab), main = main, sub = sub, theta = theta, phi = phi, r = r, 
          d= d, scale = scale, expand = expand, col = col, border = border, 
          ltheta = ltheta, lphi = lphi, shade = shade, box = box, axes = axes, 
          nticks = nticks, ticktype = ticktype, ...)
    }

  } else if (!plot){
    return(X)
  }
  
  
}