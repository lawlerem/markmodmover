% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit.R
\docType{methods}
\name{fit}
\alias{fit}
\alias{fit4M}
\alias{fit.4M}
\alias{fit,Data4M,SetModel4M-method}
\alias{fit,Data4M,missing-method}
\title{Fit a 4M model}
\usage{
\S4method{fit}{Data4M,SetModel4M}(object, model, max.tries = 10,
  DLL = "markmodmover", stationary.tolerance = 4, ...)

\S4method{fit}{Data4M,missing}(object, max.tries = 10,
  DLL = "markmodmover", stationary.tolerance = 4, ...)
}
\arguments{
\item{object}{An object of class \code{Data4M}.}

\item{model}{An object of class \code{SetModel4M}.}

\item{max.tries}{How many times shoud the fitting process be tried before stopping? Defaults to 10.}

\item{DLL}{Which DLL should be loaded? Should only be changed if the user modifies the C++ TMB template.}

\item{stationary.tolerance}{Number of decimal places that any entry of the stationary distribution is allowed to be zero without restarting the fitting procedure.}

\item{...}{Optional way to define a model without first creating an object. Also used to pass arguments to internal functions.}
}
\value{
Returns an object of class \code{Model4M} containing estimated model parameters, residuals, and Viterbi-decoded state estimates.
}
\description{
Fit a 4M model
}
\examples{
sealData<- data4M(greyseal)
sealData<- interpolate(sealData,Time.Step = 1)
seal4M2<- fit(sealData)
seal4M3<- fit(sealData,N.States = 3)

HMM3<- setModel4M(Use.HMM = T,N.States = 3)
sealHMM3<- fit(sealData,HMM3)

}
