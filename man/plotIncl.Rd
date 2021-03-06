% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plotIncl}
\alias{plotIncl}
\title{Use the image function to plot the variables included at each iteration of an SSVS model}
\usage{
plotIncl(model, col_labs = NULL, font = "serif",
  colors = c("#99004C", "#1c1c1c", "#0065CC"), ...)
}
\arguments{
\item{model}{the runjags object}

\item{col_labs}{column labels}

\item{font}{font for the axis labels. Defaults to "serif"}

\item{colors}{colors defaults to magenta, dark grey, and blue, c("#99004C", "#1c1c1c", "#0065CC")}

\item{...}{other arguments to image()}
}
\value{
an image
}
\description{
This lets you view the presence or absence of a coefficient
in the model sampeld over N MCMC chains. Positive values and negative 
values are colored differently, which lets you assess the stability
of the sign over different models.
}
\details{
An example of output: \cr
\cr
\if{html}{\figure{imageIncl.png}{}}
\if{latex}{\figure{imageIncl.png}{}}
}
\examples{
plotIncl()

}
