% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cB.R
\name{pcbinom}
\alias{pcbinom}
\title{Distribution function for the continuous binomial distribution}
\usage{
pcbinom(x, size, prob, lower.tail = TRUE)
}
\arguments{
\item{x}{vector of quantiles}

\item{size}{number of trials}

\item{prob}{probability of success on each trial}
}
\value{
the CDF/SF values
}
\description{
Distribution function for the continuous binomial distribution
}
\examples{
pcbinom(x = 4, size = 8, prob = 0.3)
}
