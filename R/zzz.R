.onLoad <- function(libname, pkgname){
  .runjags.options <<- list(modules = "dic off", inits.warning=FALSE , rng.warning=FALSE, summary.warning = FALSE, blockcombine.warning = FALSE, blockignore.warning = FALSE)
}

.onAttach<- function(libname, pkgname){
  .runjags.options <<- list(modules = "dic off", inits.warning=FALSE , rng.warning=FALSE, summary.warning = FALSE, blockcombine.warning = FALSE, blockignore.warning = FALSE)
}
