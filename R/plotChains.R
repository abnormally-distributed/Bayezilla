#' MCMC Diagnostic Plots
#'
#' Displays traceplots and cumulative means. Draws on Edward Tufte's minimalist aesthetic
#' ideas for the design.
#'
#' Traceplots should look very
#' fuzzy and chaotic. If they look more ordered and correlated that indicates a problem
#' with the sampler or model specification.
#'
#' Cumulative mean plots display how quickly the model converges onto the expected value.
#' The chains should all converge onto the same value or very close value very quickly.
#' If they do not, they indicates you either need to draw more samples, use a longer
#' adaptation and/or burnin/warmup for the MCMC chains, or perhaps reconsider the model
#' if those do not work.
#'
#' @param fit the stanfit or runjags object
#' @param pars the names of the parameters to be plotted.
#' @param col the color scheme. One of "blue" (the default), "darkblue", "purples", "reds",
#' "greens", or "mix"
#'
#' @return
#' a base R plot
#' @export
#'
#' @examples
#' plotChains()
plotChains = function(fit, pars = c("Intercept", "sigma"), col = "blues") {
  stan <- inherits(fit, "stanfit")
  if (stan == TRUE) {
    codaObject <- rstan::As.mcmc.list(fit, pars = pars)
  }
  else if (class(fit) == "runjags"){
    codaObject<- runjags::combine.mcmc(fit, collapse.chains = FALSE, vars = pars)
  }
  getCummean <- function (codaObject){
    cumulative.mean = function(d) {cumsum(d)/seq_along(d)}
    codaObject2 <- as.mcmc.list(lapply(codaObject, function(d) {as.mcmc(apply(d, 2, cumulative.mean))}))
    as.mcmc.list(codaObject2)
  }
  codaObject2 <- getCummean(codaObject)
  parNames = varnames(codaObject)

  if (col == "blues" || col == "blue"){
    ColorScheme = c("#6a9bc38A", "#000000F5","#4184e1DB", "#00bfffBF")
  }
  if (col == "darkblues" || col == "darkblue"){
    ColorScheme = c("#d1e1ec8A", "#6497b1F5","#03396cDB", "#011f4bBF")
  }
  if (col == "purples" || col == "purple"){
    ColorScheme = c("#e5cce58A", "#a64ca6F5","#660066DB", "#400040BF")
  }
  if (col == "reds" || col == "red"){
    ColorScheme = c("#DCBCBCDB", "#7C0000B0", "#8F27278A", "#B97C7C8F")
  }
  if (col == "greens" || col == "green"){
    ColorScheme= c("#d9f2e6B0", "#66cc99ED", "#2d8659BF", "#194d33A1")
  }
  if (col == "mix" || col == "multi"){
    ColorScheme= c("#11ff11", "#000000", "#ff0000B0", "#0061beED")
  }

  par( mar=0.5+c(3,3,.75, .10) , oma=0.1+c(0,0,1.25,.20) , mgp=c(2, 0.125, 0) ,
       cex.lab=1, cex.main = 1.35, family = 'serif')
  layout(matrix(1:4,nrow=2))
  require(coda)
  for (i in 1:length(parNames)){
    coda::traceplot(codaObject[, pars[i]], xaxt = "n", yaxt = "n", lwd = .125, lty = 1, bty="n", cex.lab = 1.1, cex.axis = .65, type = "l",  main= noquote(parNames[i]), ylab="Value" , col=ColorScheme)
    axis(1, col = NA, tck = 0)
    axis(2, col = NA, tck = 0)
    coda::traceplot(codaObject2[,pars[i]], xaxt = "n", yaxt = "n", lwd = 2.5, lty = 1, bty="n" , cex.lab = 1.1, cex.axis = .65, type = "l", main= noquote(parNames[i]), ylab="Cumulative Mean" , col=ColorScheme)
    axis(1, col = NA, tck = 0 )
    axis(2, col = NA, tck = 0)
    }

}
