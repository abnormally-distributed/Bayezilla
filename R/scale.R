#' Scale a data frame or matrix
#'
#' @details an improvement of the base R scale function. Unlike R's standard scale function, this
#' allows for factor columns or character columns to be present in the data, and will simply
#' leave those untouched without throwing an error.
#' @param data a data frame or vector
#' @param scale.type one of "sd" (default) for mean centered standard deviation scaled, "medmad" for
#' median centered median absolute deviation scaled, "medsd" for median centered standard
#' deviation scaled, "mean" for mean centered but unscaled, or "med" for median centered
#' but unscaled. Also accepts 1, 2, 3, 4, and 5 respectively as input types.
#' @export
#' @return A data frame
#' @examples
#' scale(data)
#'

scale = function (data, scale.type = 1)
{

  if (isTRUE(is.vector(data))) {
    Vector = "YES"
    data = cbind.data.frame(x = data)
  } else {
    Vector = "NO"
    data = as.data.frame(data)
  }

  if (scale.type == "sd" || scale.type == 1) {
    Scale1 = function(x) {
      Mean = mean(x)
      Sd = sd(x)
      x = x - Mean
      x = x/Sd
      return(x)
    }
    ind <- sapply(data, is.numeric)
    scaled.data = data
    scaled.data[ind] <- lapply(scaled.data[ind], Scale1)

    if (Vector=="YES") {
      return(as.vector(as.matrix(scaled.data)))
    }
    else {
      return(as.data.frame(scaled.data))
    }
  }
  if (scale.type == "medmad" || scale.type == 2) {
    Scale2 = function(x) {
      Median = median(x)
      Mad = mad(x)
      x = x - Median
      x = x/Mad
      return(x)
    }
    ind <- sapply(data, is.numeric)
    scaled.data = data
    scaled.data[ind] <- lapply(scaled.data[ind], Scale2)
    if (Vector=="YES") {
      return(as.vector(as.matrix(scaled.data)))
    }
    else {
      return(as.data.frame(scaled.data))
    }
  }
  if (scale.type == "medsd" || scale.type == 3) {
    Scale3 = function(x) {
      Median = median(x)
      Mad = sd(x)
      x = x - Median
      x = x/Mad
      return(x)
    }
    ind <- sapply(data, is.numeric)
    scaled.data = data
    scaled.data[ind] <- lapply(scaled.data[ind], Scale3)
    if (Vector=="YES") {
      return(as.vector(as.matrix(scaled.data)))
    }
    else {
      return(as.data.frame(scaled.data))
    }
  }
  if (scale.type == "mean" || scale.type == 4) {
    Scale4 = function(x) {
      Mean = mean(x)
      x = x - Mean
      return(x)
    }
    ind <- sapply(data, is.numeric)
    scaled.data = data
    scaled.data[ind] <- lapply(scaled.data[ind], Scale4)
    if (Vector=="YES") {
      return(as.vector(as.matrix(scaled.data)))
    }
    else {
      return(as.data.frame(scaled.data))
    }
  }
  if (scale.type == "med" || scale.type == 5) {
    Scale5 = function(x) {
      Median = median(x)
      x = x - Median
      return(x)
    }
    ind <- sapply(data, is.numeric)
    scaled.data = data
    scaled.data[ind] <- lapply(scaled.data[ind], Scale5)
    if (Vector=="YES") {
      return(as.vector(as.matrix(scaled.data)))
    }
    else {
      return(as.data.frame(scaled.data))
    }
  }
}
