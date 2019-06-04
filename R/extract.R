#' Extract a posterior to a data frame 
#' 
#' This lets you extract the contents of a posterior to a data frame. You can
#' select parameters to keep and drop.
#'
#' @param fit the model
#' @param pars parameters to keep. Defaults to NULL (implying anything is kept that's not in droppars)
#' @param droppars parameters to drop. Defaults to  c("log_lik", "ySim", "yTest", "deviance", "omega", "lambda")
#'
#' @return
#' a data frame
#' @export
#'
#' @examples
#' extract(out)
extractPost =  function (fit, pars = NULL, droppars = c("log_lik", "ySim", "yTest", "deviance", "omega", "lambda")) 
{
  stan <- inherits(fit, "stanfit")
  if (stan == TRUE) {
    ss <- rstan::As.mcmc.list(fit)
  }
  else if (class(fit) == "runjags"){
    ss <- runjags::combine.mcmc(fit, collapse.chains = FALSE)
  }
  else {
    ss <- as.matrix(as.mcmc(fit))
  }
  
  ss <- as.matrix(ss)
  
  select_parameters <- function (explicit = character(), patterns = character(), complete = character()) 
  {
    stopifnot(is.character(explicit), is.character(patterns), 
              is.character(complete))
    if (!length(explicit) && !length(patterns)) 
      return(complete)
    if (length(explicit)) {
      if (!all(explicit %in% complete)) {
        not_found <- which(!explicit %in% complete)
        stop("Some 'pars' don't match parameter names: ", 
             paste(explicit[not_found], collapse = ", "))
      }
    }
    if (!length(patterns)) {
      return(unique(explicit))
    }
    
    else {
      regex_pars <- unlist(lapply(seq_along(patterns), function(j) {
        grep(patterns[j], complete, value = TRUE)
      }))
      
      if (!length(regex_pars)) 
        stop("No matches for 'regex_pars'.", call. = FALSE)
    }
    unique(c(explicit, regex_pars))
  }
  
  
  drop <- select_parameters(patterns = droppars, complete = colnames(ss))
  drop <- sapply(1:length(drop), function(i) which(colnames(ss) == drop[i]))
  ss <- ss[, -drop]
  as.data.frame(ss)
}