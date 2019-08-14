#' Use the image function to plot the variables included at each iteration of an SSVS model
#'
#' @description This lets you view the presence or absence of a coefficient
#' in the model sampeld over N MCMC chains. Positive values and negative 
#' values are colored differently, which lets you assess the stability
#' of the sign over different models. 
#'
#' @param model the runjags object
#' @param col_labs column labels
#' @param colors colors defaults to magenta, dark grey, and blue, c("#99004C", "#1c1c1c", "#0065CC")
#' @param font font for the axis labels. Defaults to "serif"
#' @param ... other arguments to image()
#'
#'
#' @export
#'
#' @examples
#' plotIncl()
#'
#' @return
#' an image
#'
#' @details 
#' 
#' An example of output: \cr
#' \cr
#' \if{html}{\figure{imageIncl.png}{}}
#' \if{latex}{\figure{imageIncl.png}{}}
#'
plotIncl <- function(model, col_labs = NULL, font = "serif", colors=c("#99004C", "#1c1c1c", "#0065CC"), ...){
  x = as.matrix(combine.mcmc(model,collapse.chains = TRUE))
  x = x[,grep(colnames(x), pattern =  "beta")]
  x = sign(x)
  if (is.null(col_labs) == FALSE){
    colnames(x) <- col_labs
  }

  min <- min(x)
  max <- max(x)
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(x))
  }

  # Reverse Y axis
  reverse <- nrow(x) : 1
  yLabels <- yLabels[reverse]
  x <- x[reverse,]

  # Data Map
  image(1:length(xLabels), 1:length(yLabels), t(x), col=colors, xlab="",
        ylab="", axes=FALSE, zlim=c(min,max))
  if( !is.null(title) ){
    title(main=title)
  }
  par(family = font)
  axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
}




#' plot inclusion probabilities from SSVS models
#'
#' @description lollipop plot with ggplot2 of variable inclusion probabilities for SSVS models.
#'
#' @param model the model
#' @param col_labs variable names
#' @param font font for the axis labels. Defaults to "serif"
#'
#' @return
#' a plot
#'
#' @export
#'
#' @examples
#' plotPips()
#'
plotPips = function(model, col_labs = NULL, font = "serif"){
  x = as.matrix(combine.mcmc(model,collapse.chains = TRUE))
  x = x[,grep(colnames(x), pattern =  "delta")]
  if (is.null(col_labs) == FALSE){
    colnames(x) <- col_labs
  }

  data = rownames_to_column(as.data.frame(colMeans(x)))
  
  colnames(data) = c("variable", "pip")
  ggplot(data = data, aes(x = variable, y = pip, color = ifelse(pip < .50, "Fail", ifelse(pip <= .60 & pip > .50, "Borderline", "Pass")), fill = ifelse(pip < .50, "Fail", ifelse(pip <= .60 & pip > .50, "Borderline", "Pass")))) +
    geom_segment( aes(x=variable, xend=variable , y=0, yend= pip - .001), size = 1.125, alpha = .50, linetype = "solid") +
    geom_point(size = 5, alpha = .90, shape = 21) +
    scale_fill_manual(values =  c("red", "#51fa51" , "#709e61"), aesthetics = "fill") +
    scale_color_manual(values =  c("darkred", "#05a905" , "#788b74"), aesthetics = "color") +
    theme(legend.position="none") +
    coord_cartesian(ylim = c(0, 1)) +
    geom_hline(yintercept = .50, colour="#990000", linetype="longdash") + 
    theme(text = element_text(family = font))
}


#' Plot distribution of num. of included vars. over iterations of an SSVS model
#'
#'
#' @param model a runjags object with inclusion indicators labeled "delta"
#' @param binwidth defaults to .5
#'
#' @return
#' a plot
#' @export
#'
#' @examples
#' plotNumvar()
#'
plotNumvar = function(model, binwidth = .5){
  x = as.matrix(combine.mcmc(model,collapse.chains = TRUE))
  x = x[,grep(colnames(x), pattern =  "delta")]
  x = rowSums(x)
  ggplot(data = data.frame(num.variables = as.factor(x)), aes(x=num.variables)) +
    geom_bar(width = binwidth, fill = "darkblue")
}

