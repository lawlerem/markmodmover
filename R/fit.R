#' Fit a 4M model
#'
#' @param object An object of class \code{Data4M}.
#' @param model An object of class \code{SetModel4M}.
#'
#' @name fit
#' @aliases fit4M fit.4M
NULL

#' @export
#' @noRd
setGeneric(name = "fit",
           def = function(object,
                          model,
                          ...) standardGeneric("fit")
)

#' @param max.tries How many times shoud the fitting process be tried before stopping? Defaults to 10.
#' @param DLL Which DLL should be loaded? Should only be changed if the user modifies the C++ TMB template.
#' @param stationary.tolerance Number of decimal places that any entry of the stationary distribution is allowed to be zero without restarting the fitting procedure.
#' @param ... Optional way to define a model without first creating an object. Also used to pass arguments to internal functions.
#'
#' @return Returns an object of class \code{Model4M} containing estimated model parameters, residuals, and Viterbi-decoded state estimates.
#'
#' @examples
#' sealData<- data4M(greyseal)
#' sealData<- interpolate(sealData,Time.Step = 1)
#' seal4M2<- fit(sealData)
#' seal4M3<- fit(sealData,N.States = 3)
#'
#' HMM3<- setModel4M(Use.HMM = T,N.States = 3)
#' sealHMM3<- fit(sealData,HMM3)
#'
#' @export
#' @rdname fit
setMethod(f = "fit",
          signature = c(object = "Data4M",
                        model = "SetModel4M"),
          definition = function(object,
                                model,
                                max.tries = 10,
                                DLL = "markmodmover",
                                stationary.tolerance = 4,
#                                convergence.error = T,
                                ...) {

                Turning.Angle.Data<- pi/180 * movementData(object)$Movement.Data$Deflection.Angle
                Step.Length.Data<- movementData(object)$Movement.Data$Step.Length
                Step.Length.Starting.Values<- movementData(object)$Step.Length.Starting.Values
                Step.Distribution<- grep(distribution(model),
                                         c("gamma",
                                           "log-normal"))
                Data.Grouping<- c(0,
                                  cumsum(rle(as.numeric(groups(object)))$lengths)
                                  )
                data<- list(theta = Turning.Angle.Data,
                            st_dist = Step.Length.Data,
                            st_dist_starts = Step.Length.Starting.Values,
                            grouping = Data.Grouping,
                            step_distribution = Step.Distribution,
                            step_zero_inflation = as.integer(zeroInflation(model)[["Step.Length"]]),
                            angle_zero_inflation = as.integer(zeroInflation(model)[["Deflection.Angle"]])
                            )


                Internal.Environment<- environment()



                for( i in seq(max.tries) ) {
                        Tpm.Working.Pars<- startingValues(model)$Tpm.Working.Pars
                        Tpm.Working.Pars[is.na(Tpm.Working.Pars)]<- runif(sum(is.na(Tpm.Working.Pars)),
                                                                          min = -2,
                                                                          max = 2)

                        Theta.Working.Pars<- startingValues(model)$Theta.Working.Pars
                        Theta.Working.Pars[is.na(Theta.Working.Pars)]<- runif(sum(is.na(Theta.Working.Pars)),
                                                                              min = -2,
                                                                              max = 2)

                        Step.Working.Pars<- startingValues(model)$Step.Working.Pars
                        Step.Working.Pars[is.na(Step.Working.Pars)]<- runif(sum(is.na(Step.Working.Pars)),
                                                                            min = -2,
                                                                            max = 2)

                        Logit.Step.Zero.Probs<- startingValues(model)$Logit.Step.Zero.Probs
                        Logit.Step.Zero.Probs[is.na(Logit.Step.Zero.Probs)]<- runif(sum(is.na(Logit.Step.Zero.Probs)),
                                                                                    min = -2,
                                                                                    max = 2)

                        Logit.Angle.Zero.Probs<- startingValues(model)$Logit.Angle.Zero.Probs
                        Logit.Angle.Zero.Probs[is.na(Logit.Angle.Zero.Probs)]<- runif(sum(is.na(Logit.Angle.Zero.Probs)),
                                                                                    min = -2,
                                                                                    max = 2)

                        Angle.Zero.Working.Pars<- startingValues(model)$Angle.Zero.Working.Pars
                        Angle.Zero.Working.Pars[is.na(Angle.Zero.Working.Pars)]<- runif(sum(is.na(Angle.Zero.Working.Pars)),
                                                                                    min = -2,
                                                                                    max = 2)


                        map<- parameterMapping(model)
                        if( useHMM(model) == T ){
                                if( distribution(model) == "gamma" ) {
                                        Step.Working.Pars[,2]<- -Inf   ### autocorrelation to 0
                                } else if( distribution(model) == "log-normal" ) {
                                        Step.Working.Pars[,2]<- 0
                                }
                                map[["Step.Working.Pars"]][,2]<- NA ### Fix autocorrelation
                        } else {}

                        if( !zeroInflation(model)[["Step.Length"]] ) {
                                Logit.Step.Zero.Probs[]<- -Inf
                                map[["Logit.Step.Zero.Probs"]][]<- NA
                        } else {}

                        if( !zeroInflation(model)[["Deflection.Angle"]] ) {
                                Logit.Angle.Zero.Probs[]<- -Inf
                                map[["Logit.Angle.Zero.Probs"]][]<- NA

                                Angle.Zero.Working.Pars[,1]<- 0
                                Angle.Zero.Working.Pars[,2]<- -Inf
                                map[["Angle.Zero.Working.Pars"]][]<- NA
                        }

                        parameters<- list(tpm_working_pars_matrix = Tpm.Working.Pars,
                                          theta_working_pars = Theta.Working.Pars,
                                          dist_working_pars = Step.Working.Pars[,-2],
                                          acf_working_pars = Step.Working.Pars[,2],
                                          logit_step_zero_probs = Logit.Step.Zero.Probs,
                                          logit_angle_zero_probs = Logit.Angle.Zero.Probs,
                                          angle_zero_working_pars = Angle.Zero.Working.Pars,
                                          dummy = 0)



                        map<- list(tpm_working_pars_matrix = map$Tpm.Working.Pars,
                                   theta_working_pars = map$Theta.Working.Pars,
                                   dist_working_pars = map$Step.Working.Pars[,-2],
                                   acf_working_pars = map$Step.Working.Pars[,2],
                                   logit_step_zero_probs = map$Logit.Step.Zero.Probs,
                                   logit_angle_zero_probs = map$Logit.Angle.Zero.Probs,
                                   angle_zero_working_pars = map$Angle.Zero.Working.Pars)

                        map<- lapply(map,as.factor)


                        Unordered.ADFun<- TMB::MakeADFun(data = data,
                                                  parameters = parameters,
                                                  map = c(map,list(dummy = as.factor(NA))),
                                                  DLL = DLL,
                                                  silent = T)

                        Unordered.Optim<- logical()
                        tryCatch({
                                Unordered.Optim<- nlminb(Unordered.ADFun$par,
                                                         Unordered.ADFun$fn,
                                                         Unordered.ADFun$gr)
                                assign(x = "Unordered.Optim",
                                       value = Unordered.Optim,
                                       envir = Internal.Environment)
                                if( Unordered.Optim$convergence != 0 ) {
                                        warning("Optimization was not convergent.")
                                } else {}
                                if( !(0 %in% round(Unordered.ADFun$report()$start_probs,
                                                   digits = stationary.tolerance)
                                     )
                                  ) {
                                  break()
                                }
                        },
                        warning = function(e) {},
                        error = function(e) {}
                        )
                }

                if( class(Unordered.Optim) == "logical" ) {
#                        if( convergence.error == T ) {
                                stop("Optimization did not work. Try increasing max.tries.")
#                        } else {
#                                Unordered.Optim<- nlminb(Unordered.ADFun$par,
#                                                         Unordered.ADFun$fn,
#                                                         Unordered.ADFun$gr)
#                                assign(x = "Unordered.Optim",
#                                       value = Unordered.Optim,
#                                       envir = Internal.Environment)
#                        }
                }


                Unordered.Report<- Unordered.ADFun$report()


                State.Order<- order(Unordered.Report$theta_pars[,2]) ### order by angle concentration parameter

                Tpm.Working.Pars<- Unordered.Report$tpm_working_pars_matrix[State.Order,State.Order]
                Theta.Working.Pars<- Unordered.Report$theta_working_pars[State.Order,]
                Step.Working.Pars[,c(1,3)]<- Unordered.Report$dist_working_pars[State.Order,]
                Step.Working.Pars[,2]<- Unordered.Report$acf_working_pars[State.Order]
                Logit.Step.Zero.Probs<- Unordered.Report$logit_step_zero_probs[State.Order]
                Logit.Angle.Zero.Probs<- Unordered.Report$logit_angle_zero_probs[State.Order]
                Angle.Zero.Working.Pars<- Unordered.Report$angle_zero_working_pars[State.Order,]

                parameters<- list(tpm_working_pars_matrix = Tpm.Working.Pars,
                                          theta_working_pars = Theta.Working.Pars,
                                          dist_working_pars = Step.Working.Pars[,-2],
                                          acf_working_pars = Step.Working.Pars[,2],
                                          logit_step_zero_probs = Logit.Step.Zero.Probs,
                                          logit_angle_zero_probs = Logit.Angle.Zero.Probs,
                                          angle_zero_working_pars = Angle.Zero.Working.Pars,
                                          dummy = 0)

                map<- lapply(map,function(x) {
                                        x[]<- NA
                                        x<- as.factor(x)
                                        x<- droplevels(x)
                                        return(x)
                                 }
                            )



                Ordered.ADFun<- TMB::MakeADFun(data = data,
                                          parameters = parameters,
                                          map = map,
                                          DLL = DLL,
                                          silent = T)
                Ordered.Optim<- nlminb(Ordered.ADFun$par,
                                       Ordered.ADFun$fn,
                                       Ordered.ADFun$gr)
                Report<- Ordered.ADFun$report()

                Parameter.Estimates<- list(
                        Deflection.Angle.Parameters = data.frame(Center = Report$theta_pars[,1],
                                                                 Concentration = Report$theta_pars[,2]),
                        Step.Length.Parameters = data.frame(Intercept = Report$dist_pars[,1],
                                                            Autocorrelation = Report$acf_pars,
                                                            Standard.Deviation = Report$dist_pars[,2]),
                        Transition.Probability.Matrix = Report$tpm,
                        Stationary.Distribution = Report$start_probs,
                        Zero.Inflation = list(
                                Step.Zero.Probs = Report$step_zero_probs,
                                Angle.Zero.Probs = Report$angle_zero_probs,
                                Angle.Zero.Pars = data.frame(
                                        Center = Report$angle_zero_pars[,1],
                                        Concentration = Report$angle_zero_pars[,2]
                                )
                        )
                )

                Model.4M<- model4M(object,model)

                parameterEstimates(Model.4M)<- Parameter.Estimates
                aic(Model.4M)<- 2*Ordered.Optim$objective + 2*(length(Unordered.Optim$par)-1)  ### Minus 1 to account for dummy = 0
                viterbiPath(Model.4M)<- as.integer(Report$viterbi_path)
                convergence(Model.4M)<- Unordered.Optim$message
                tmbEnvironment(Model.4M)<- Ordered.ADFun$env

                if( Unordered.Optim$convergence != 0 ) {
                        warning("Optimization was not convergent.")
                }


                Data.Index<- seq(nrow(movementData(Model.4M)$Movement.Data))
                Theta.Residuals<- sapply(Data.Index,function(k) {
                        quantile<- sum(Report$forecast[,k] %*% (Report$theta_cdf_array[,,k] %% 1))
                        return(2*quantile - 1)
                })
                Step.Residuals<- sapply(Data.Index,function(k) {
                        quantile<- sum(Report$forecast[,k] %*% Report$dist_cdf_array[,,k])
                        return(2*quantile-1)
                })

                residuals(Model.4M)<- data.frame(Deflection.Angle = Theta.Residuals,
                                                Step.Length = Step.Residuals)

                return(Model.4M)
          }
)
#' @export
#' @rdname fit
setMethod(f = "fit",
          signature = c(object = "Data4M",
                        model = "missing"),
          definition = function(object,
                                max.tries = 10,
                                DLL = "markmodmover",
                                stationary.tolerance = 4,
                                ...) {
                Set.Model.4M<- setModel4M(...)
                return(fit(object,Set.Model.4M))
          }
)
#' @export
#' @noRd
setMethod(f = "fit",
          signature = c(object = "Model4M",
                        model = "SetModel4M"),
          definition = function(object,
                                model,
                                ...) {
                return(fit(data4M(object),model,...))
          }
)
#' @export
#' @noRd
setMethod(f = "fit",
          signature = c(object = "Model4M",
                        model = "missing"),
          definition = function(object,
                                ...) {
                return(fit(data4M(object),SetModel4M(object),...))
          }
)
