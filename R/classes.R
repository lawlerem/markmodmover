#' An S4 class to hold and pre-process movement data.
#'
#' @family 4M-classes
#'
#' @slot Identification A string, list, number, etc. used for identification.
#' @slot Observed.Locations A data frame with columns \code{Date}, \code{Lon}, and \code{Lat}.
#' @slot Interpolation.Parameters A named vector with positive entries \code{Time.Step} and \code{Group.Cutoff}.
#' @slot Interpolated.Locations A data frame with columns \code{Date}, \code{Lon}, \code{Lat}, and \code{Group}.
#' @slot Movement.Data A named list holding the pre-processed data with the following entries:
#'   \describe{
#'     \item{Movement.Data}{A data frame with columns \code{Deflection.Angle}, \code{Step.Length}, and \code{Group}.}
#'     \item{Step.Length.Starting.Values}{A numeric vector giving the first step length for each group.}
#'     \item{Step.Length.Mean}{The mean of the unstandardized step lengths.}
#'   }
#'
#' @export
setClass(
        Class = "Data4M",
        slots = c(Identification = "ANY",
                  Observed.Locations = "data.frame",
                  Interpolation.Parameters = "numeric",
                  Interpolated.Locations = "data.frame",
                  Movement.Data = "list"
        )
)


#' An S4 class to hold model definition and fitting parameters.
#'
#' @family 4M-classes
#'
#' @slot N.States The number of states \code{n} used in the model.
#' @slot Use.HMM If \code{TRUE} fixes autocorrelation to zero.
#' @slot Distribution Distribution for step length. Can be one of \code{gamma} or \code{log-normal}.
#' @slot Zero.Inflation A named logical vector with entries \code{Step.Length} and \code{Deflection.Angle}. A \code{TRUE} entry enables zero inflation.
#' @slot Starting.Values A named list giving the starting values of the working parameters with the following entries:
#' \describe{
#'   \item{Tpm.Working.Pars}{An \code{n} by \code{n} matrix.}
#'   \item{Theta.Working.Pars}{An \code{n} by 2 matrix. The first column is untransformed center, the second is the logit of concentration.}
#'   \item{Step.Working.Pars}{An \code{n} by 3 matrix. The first column is the log of the intercept, the second is the log of the autocorrelation, the third is the log of the standard deviation.}
#'   \item{Logit.Step.Zero.Probs}{A length \code{n} vector.}
#'   \item{Logit.Angle.Zero.Probs}{A length \code{n} vector.}
#'   \item{Angle.Zero.Working.Pars}{An \code{n} by 2 matrix.}
#' }
#' @slot Parameter.Mapping A named list with the same structure as \code{Starting.Values}. Controls TMB parameter mapping. See details of \code{\link[TMB]{MakeADFun}}.
#'
#' @export
setClass(
        Class = "SetModel4M",
        slots = c(
                N.States = "integer",
                Use.HMM = "logical",
                Distribution = "character",
                Zero.Inflation = "logical",
                Starting.Values = "list",
                Parameter.Mapping = "list"
        )
)
.implemented.distributions<- c("gamma","log-normal")



#' An S4 class for fitted models from \linkS4class{Data4M} and \linkS4class{SetModel4M}.
#'
#' @family 4M-classes
#'
#' @section Extends:
#' \code{\link{Data4M-class}}, \code{\link{SetModel4M-class}}.
#'
#' @slot Parameter.Estimates A named list holding parameter estimates with the following entries:
#' \describe{
#'   \item{Deflection.Angle.Parameters}{A data.frame with columns \code{Center} and \code{Concentration}.}
#'   \item{Step.Length.Parameters}{A data.frame with columns \code{Intercept}, \code{Autocorrelation}, and \code{Standard.Deviation}.}
#'   \item{Transition.Probabilities}{A square matrix of transition probabilities.}
#'   \item{Stationary.Distribution}{A numeric vector giving the stationary distribution of \code{Transition.Probabilties}.}
#'   \item{Zero.Inflation}{A named list containing zero-inflation parameters.}
#' }
#' @slot AIC A numeric value giving the AIC of the fitted model.
#' @slot Residuals A data.frame with columns \code{Deflection.Angle} and \code{Step.Length} giving probability-scale forecast residuals.
#' @slot Viterbi.Path A integer vector giving the Viterbi decoded state estimates.
#' @slot Convergence A character vector giving the convergence message.
#' @slot TmbEnvironment The environment code for the TMB \code{MakeADFun} object.
#'
#' @export
setClass(
        Class = "Model4M",
        slots = c(Parameter.Estimates = "list",
                  AIC = "numeric",
                  Residuals = "data.frame",
                  Viterbi.Path = "integer",
                  Convergence = "character",
                  TmbEnvironment = "environment"
                  ),
        contains = c("Data4M",
                     "SetModel4M")
)

#' An S4 class for holding simulations from a \code{\link{Model4M}} class.
#'
#' @family 4M-classes
#'
#' @section Extends:
#' \code{\link{Model4M-class}}.
#'
#' @slot Simulated.Locations A 3d array holding the simulated locations with columns \code{Date}, \code{Lon}, \code{Lat}, and \code{Group}, and depth corresponding to individual simulations.
#' @slot Simulated.Data A 3d array holding the simulated movement data with columns \code{Deflection.Angle}, \code{Step.Length}, and \code{Group}, and depth corresponding to individual simulations.
#' @slot Simulated.Viterbi.Path A matrix holding the simulated state path with rows corresponding to individual simulations.
#' @slot Refit.Parameters A named list holding refit parameter estimates with the following entries:
#' \describe{
#'   \item{Deflection.Angle.Parameters}{A 3d array with rows corresponding to states, columns corresponding to parameters, and depth corresponding to individual simulations.}
#'   \item{Step.Length.Parameters}{A 3d array with rows corresponding to states, columns corresponding to parameters, and depth corresponding to individual simulations.}
#'   \item{Transition.Probabilities}{A 3d array with depth corresponding to individual simulations.}
#'   \item{Stationary.Distribution}{A matrix with rows corresponding to individual simulations and column corresponding to state.}
#'   \item{Zero.Inflation}{A named list containing zero-inflation parameters.}
#' }
#' @slot Refit.Residuals A 3d array with columns \code{Deflection.Angle} and \code{Step.Length}, and depth corresponding to individual simulations.
#' @slot Refit.Viterbi.Path A matrix holding the refit Viterbi-estimated state path with rows corresponding to individual simulations.
#' @slot Refit.Convergence A character vector holding the convergence messages from each simulation.
#' @slot Refit.Environment A list holding the TMB environments for each simulation.
#'
#' @export
setClass(Class = "Simulate4M",
         slots = c(Simulated.Locations = "array",
                   Simulated.Data = "array",
                   Simulated.Viterbi.Path = "array",
                   Refit.Parameters = "list",
                   Refit.Residuals = "array",
                   Refit.Viterbi.Path = "array",
                   Refit.Convergence = "character",
                   Refit.Environment = "list"
                   ),
         contains = "Model4M"
)



