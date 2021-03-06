% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getset.R
\docType{methods}
\name{simulate4M}
\alias{simulate4M}
\alias{initialize,Simulate4M-method}
\alias{simulate4M,Model4M-method}
\title{Create or extract a \code{Simulate4M} object.}
\usage{
simulate4M(x, ...)

\S4method{initialize}{Simulate4M}(.Object, Model4M = new("Model4M"))

\S4method{simulate4M}{Model4M}(x)
}
\arguments{
\item{x}{An object to coerce to class \code{Simulate4M}, possible of class \code{Model4M}.}
}
\description{
Create or extract a \code{Simulate4M} object.
}
\details{
Consider using \code{\link[=simulate.4M]{simulate}} to create new objects of class \code{Simulate4M}.

The \code{initialize} method is not meant to be used by the user. Use \code{simulate4M} instead.
}
\seealso{
Other Construct4M: \code{\link{data4M}},
  \code{\link{model4M}}, \code{\link{setModel4M}},
  \code{\link{setSim4M}}
}
\concept{Construct4M}
