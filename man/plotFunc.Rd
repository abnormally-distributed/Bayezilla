% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotFunc.R
\name{plotFunc}
\alias{plotFunc}
\title{plot a function}
\usage{
plotFunc(func, from = NULL, to = NULL, n = 1500,
  xname = "<U+03C7>", xlab = xname, ylab = NULL, xlim = NULL,
  args = list(), color = "#005d8bCC", size = 0.75, ...)
}
\arguments{
\item{func}{The name of a function}

\item{from}{the lower limit of the x-axis}

\item{to}{the upper limit of the x-axis}

\item{n}{the number of x values to be evaluated. Defaults to 1500.}

\item{xname}{character; argument of the function in func, the name of the x axis}

\item{xlab}{label of the x-axis}

\item{ylab}{label of the y-axis}

\item{xlim}{restrict the range of the function}

\item{args}{list of the additional argument of the function}

\item{color}{the color of the line. defaults to blue}

\item{size}{the size of the line. defaults to .75}

\item{...}{Additional arguments to stat_function()}
}
\description{
plot a function
}
\examples{
plotFunc(dnorm, from = -8 to = 8)
}
