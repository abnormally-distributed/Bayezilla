.onLoad <- function(libname, pkgname){
  .runjags.options <<- list(modules = "dic off", inits.warning=FALSE , rng.warning=FALSE, summary.warning = FALSE, blockcombine.warning = FALSE, blockignore.warning = FALSE)
  invisible()
}

.onAttach<- function(libname, pkgname){
  .runjags.options <<- list(modules = "dic off", inits.warning=FALSE , rng.warning=FALSE, summary.warning = FALSE, blockcombine.warning = FALSE, blockignore.warning = FALSE)
  invisible()
  }

.onLoad <- function(libname, pkgname){
  par(family = "serif", font.lab = 2, font.main = 2)
  invisible()
}

.onAttach<- function(libname, pkgname){
    par(family = "serif", font.lab = 2, font.main = 2)
    invisible()
}

.onLoad <- function(libname, pkgname) {
  ggplot2::theme_set(Bayezilla:::theme_min())
  invisible()
}

.onAttach<- function(libname, pkgname){
  suppressMessages(suppressWarnings(assignInNamespace("hist.default", Bayezilla:::hist.default, "Bayezilla")))
}