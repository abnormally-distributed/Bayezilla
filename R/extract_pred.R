#' Extract predicted y values (posterior predictions)
#'
#' @param x the stanfit or runjags object containing posterior predictions
#' @param predname default name is "ySim".
#' @export
#' @examples
#' extract_post_pred()
#'
extract_pred =function (x, predname = "ySim")
{
  stan <- inherits(x, "stanfit")
  if (stan == TRUE) {
    ss <- as.matrix(x, pars = predname)
    extrac
  }
  else if (class(x) == "runjags"){
    ss <- runjags::combine.mcmc(x, collapse.chains = TRUE, vars = predname)
    ss <- as.matrix(ss)
  }

  df = as.matrix(ss)
  rownames(df) = c(paste0("yRep[", 1:nrow(df), "]"))
  return(df)
}
