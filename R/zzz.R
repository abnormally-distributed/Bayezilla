.onLoad <- function(libname, pkgname){
  invisible(require(runjags))
  invisible(runjags::runjags.options(list(modules = "dic off", inits.warning=FALSE , rng.warning=FALSE, summary.warning = FALSE, blockignore.warning = FALSE)))
}
