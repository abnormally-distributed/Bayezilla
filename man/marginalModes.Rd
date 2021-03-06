% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary_utilities.R
\name{marginalModes}
\alias{marginalModes}
\title{Get the marginal modes of a posterior distribution}
\usage{
marginalModes(x, keeppars = NULL, droppars = NULL)
}
\arguments{
\item{x}{set of posterior samples to be summarized in a runjags or stanfit object}

\item{keeppars}{list of specific variables to keep}

\item{droppars}{list of parameters to exclude}
}
\value{
a vector of modes
}
\description{
Get the marginal modes of a posterior distribution
}
\examples{
marginalModes(fit)
}
