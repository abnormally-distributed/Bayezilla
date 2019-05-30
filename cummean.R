getCummean <- function (mcmcout){
  cumulative.mean = function(d) {cumsum(d)/seq_along(d)}
  codaObject <- as.mcmc.list(lapply(mcmcout, function(z) {as.mcmc(apply(z, 2, cumulative.mean))}))
  return (codaObject)
}
