% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genmet.R
\docType{methods}
\name{postACF}
\alias{postACF}
\alias{postACF,Data4M-method}
\title{Estimate the autocorrelation function before interpolating locations.}
\usage{
postACF(object, ...)

\S4method{postACF}{Data4M}(object, lag.max = 30)
}
\arguments{
\item{object}{An object of or inheriting from class \code{Data4M}.}

\item{lag.max}{The maximum lag to compute.}
}
\value{
An object of class \code{acf}, which is a list containing the following elements:
 \itemize{
     \item lag: A vector containing the lags at which the acf is estimated.
     \item acf: An vector with the same length as 'lag' containing the estimated acf.
     \item type: Always set to  "correlation".
     \item n.used A vector giving the sum of kernel values for each lag.
     \item series The name of the series 'object'.
     \item snames Always set to NULL.
 }
}
\description{
Estimate the acf for step length after interpolating the data.
}
