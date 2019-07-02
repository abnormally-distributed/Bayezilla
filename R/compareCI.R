#' Check the sensitivity of credible intervals for the same data fit to different models to different priors
#'
#' @param list a list of tidy summaries where the summaries are all of the same model with the same parameters. This
#' may require some preparation on your part. 
#' @param adj the amount of space between intervals, deaults to .75
#' @param model.names an optional character vector of model names.
#' @export
#' @examples
#' s1 = post_summary(model1, keeppars = c("Intercept", "beta"))
#' s2 = post_summary(model2, keeppars = c("Intercept", "beta"))
#' s3 = post_summary(model3, keeppars = c("Intercept", "beta"))
#' compareCI(list(s1, s2, s3))
#'
compareCI = function(list, model.names="default", adj=.75){
  
  if (model.names=="default") {
    Model = as.character(c(sapply(1:length(list), function(n) rep(paste0("model", n), nrow(list[[1]])))))
  } else {
    Model = as.character(c(sapply(1:length(model.names), function(n) rep(model.names[n], nrow(list[[1]])))))
  }
  
  posts = as.data.frame(dplyr::bind_rows(list))
  
  low = which(colnames(posts)=="conf.low" | colnames(posts)=="cred.low")
  high = which(colnames(posts)=="conf.high" | colnames(posts)=="cred.high")
  
  colnames(posts)[low]  = "lower"
  colnames(posts)[high] = "upper"
  
  posts = cbind.data.frame(Model= as.factor(Model), posts)
  posts %>%
    ggplot(aes(y = forcats::fct_rev(term), x = estimate, xmin = lower, xmax = upper, color = Model)) +
    labs(y="term", x="\u03B8") +
    geom_point(size=1.6, position=ggstance::position_dodgev(height = adj+.001, preserve="total")) +
    geom_errorbarh(height = .5, size=.65, position=ggstance::position_dodgev(height = adj, preserve="total")) 
}
