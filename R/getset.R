#' @include classes.R geodesy.R
NULL



######################
###                ###
###  Construct4M   ###
###                ###
######################

#' Create or extract a \code{Data4M} object.
#'
#' Takes a data.frame of time-stamped locations and store it in an object of class \code{Data4M} for further processing. Otherwise takes an object and coerces it to class Data4M.
#' 
#' @param x An object to coerce to class \code{Data4M}. If a data.frame, it must contain at least the columns \code{Date}, \code{Longitude}, and \code{Latitude}. These column names can be shortened e.g., \code{Lat}.
#'
#' @family Construct4M
#' @export
setGeneric(name = "data4M",
           def = function(x,...) standardGeneric("data4M"))
#' @details The \code{initialize} function is not meant to be used by the user, use \code{data4M} instead.
#'
#' @export
#' @rdname data4M
setMethod(
        f = "initialize",
        signature = "Data4M",
        definition = function(.Object,
                              Identification = list(),
                              Data = data.frame(
                                        Date = .POSIXct(numeric(0)),
                                        Lon = numeric(0),
                                        Lat = numeric(0)
                                     )) {
                identification(.Object)<- Identification
                
                
                ### Find the column for Date
                ###
                Data.Date<- grep("Date",colnames(Data),ignore.case = T)
                if( length(Data.Date) == 0 ) {
                        stop("Could not find column ``Date``.")
                } else if ( length(Data.Date) > 1 ) {
                        warning(paste0("Multiple columns corresponding to ``Date`` found, using column ",
                                       Data.Date[[1]],"."))
                        Data.Date<- Data.Date[[1]]
                } else {}

                if( !("POSIXct" %in% class(Data[,Data.Date])) ) {
                        Data[,Data.Date]<- as.POSIXct(Data[,Data.Date])
                } else {}
                
                ###
                ### Find the column for Longitude
                ###
                Data.Lon<- grep("Lon",colnames(Data),ignore.case = T)
                if( length(Data.Lon) == 0 ) {
                        stop("Could not find column ``Longitude``.")
                } else if ( length(Data.Lon) > 1 ) {
                        warning(paste0("Multiple columns corresponding to ``Longitude`` found, using column ",
                                       Data.Lon[[1]],"."))
                        Data.Lon<- Data.Lon[[1]]
                } else {}
                
                ###
                ### Find the column for Latitude
                ###
                Data.Lat<- grep("Lat",colnames(Data),ignore.case = T)
                if( length(Data.Lat) == 0 ) {
                        stop("Could not find column ``Longitude``.")
                } else if ( length(Data.Lat) > 1 ) {
                        warning(paste0("Multiple columns corresponding to ``Latitude`` found, using column ",
                                       Data.Lat[[1]],"."))
                        Data.Lat<- Data.Lat[[1]]
                } else {}
                
                
                observedLocations(.Object)<- data.frame(
                        Date = Data[,Data.Date],
                        Lon = Data[,Data.Lon],
                        Lat = Data[,Data.Lat]
                )
                
                interpolationParameters(.Object)<- c(
                        Time.Step = numeric(0),
                        Group.Cutoff = numeric(0)
                )
                
                interpolatedLocations(.Object)<- data.frame(
                        Date = .POSIXct(numeric(0)),
                        Lon = numeric(0),
                        Lat = numeric(0),
                        Group = factor()
                )
                
                movementData(.Object)<- list(
                        Movement.Data = data.frame(
                                Deflection.Angle = numeric(0),
                                Step.Length = numeric(0),
                                Group = factor()
                        ),
                        Step.Length.Starting.Values = numeric(0),
                        Step.Length.Mean = numeric(0)
                )
              
              return(.Object)
        }                                          
)
#' @param Identification An optional argument of any type to help the user identify the data set.
#'
#' @export
#' @rdname data4M
setMethod(f = "data4M",
          signature = "data.frame",
          definition = function(x,Identification = list()) {
                return(new("Data4M",
                           Data = x,
                           Identification = Identification))
          }
)
#' @export
#' @rdname data4M
setMethod(f = "data4M",
          signature = "Model4M",
          definition = function(x) return(as(x,"Data4M"))
)






#' Create or extract a \code{SetModel4M} object.
#'
#' Creates or extracts an object used to hold a model definition.
#'
#' @param x An object to coerce to class \code{SetModel4M}. Missing in calls to create a new object.
#'
#' @family Construct4M
#' @export
setGeneric(name = "setModel4M",
           def = function(x,...) standardGeneric("setModel4M")
)
#' @details The \code{initialize} function is not meant to be used by the user, use \code{setModel4M} instead.
#'
#' @export
#' @rdname setModel4M
setMethod(
        f = "initialize",
        signature = "SetModel4M",
        definition = function(.Object,
                              N.States = 2,
                              Use.HMM = FALSE,
                              Distribution = "gamma",
                              Zero.Inflation = c(Step.Length = FALSE,
                                                 Deflection.Angle = FALSE),
                              Starting.Values = list(),
                              Parameter.Mapping = list()
                              ) {
                nStates(.Object)<- as.integer(N.States)
                useHMM(.Object)<- Use.HMM
                distribution(.Object)<- grep(Distribution,
                                             .implemented.distributions,
                                             ignore.case = T,
                                             value = T)
                zeroInflation(.Object)<- c(Step.Length = FALSE,
                                           Deflection.Angle = FALSE)
                ###
                ### Initialize everything to NA (random starting values)
                ###
                .Object@Starting.Values<- list(
                        Tpm.Working.Pars = matrix(NA,
                                                  nrow = N.States,
                                                  ncol = N.States),
                        Theta.Working.Pars = matrix(NA,
                                                    nrow = N.States,
                                                    ncol = 2),
                        Step.Working.Pars = matrix(NA,
                                                   nrow = N.States,
                                                   ncol = 3),
                        Logit.Step.Zero.Probs = rep(NA,N.States),
                        Logit.Angle.Zero.Probs = rep(NA,N.States),
                        Angle.Zero.Working.Pars = matrix(NA,
                                                         nrow = N.States,
                                                         ncol = 2)
                )
                if( length(Starting.Values) != 0 ) {
                        startingValues(.Object)<- Starting.Values
                }
                
                ###
                ### These need to be factors, but to keep the dimensions we'll convert later
                ###
                .Object@Parameter.Mapping<- list(
                        Tpm.Working.Pars = matrix(seq(N.States^2),
                                                  nrow = N.States,
                                                  ncol = N.States),
                        Theta.Working.Pars = matrix(seq(N.States*2),
                                                    nrow = N.States,
                                                    ncol = 2),
                        Step.Working.Pars = matrix(seq(N.States*3),
                                                   nrow = N.States,
                                                   ncol = 3),
                        Logit.Step.Zero.Probs = seq(N.States),
                        Logit.Angle.Zero.Probs = seq(N.States),
                        Angle.Zero.Working.Pars = matrix(seq(N.States*2),
                                                         nrow = N.States,
                                                         ncol = 2)
                )
                if( length(Parameter.Mapping) != 0 ) {
                        parameterMapping(.Object)<- Parameter.Mapping
                }
                
                return(.Object)
        }
)
#' @param N.States An integer giving number of states to use in the model.
#' @param Use.HMM A logical value which flags use of an HMM. If true, sets all step length autocorrelation parameters to 0.
#' @param Distribution Which distribution should be used for step length? Currently available options are \code{gamma} and \code{log-normal}.
#' @param Zero.Inflation A logical vector with names \code{Step.Length} and \code{Deflection.Angle}. \code{TRUE} entries enable zero inflation.
#' @param Starting.Values A list holding starting values. All starting values default to randomly picked values, see details and the \code{\link{SetModel4M}} class documentation.
#' @param Parameter.Mapping A list holding parameter mapping with the same structure as \code{Starting.Values}. See the TMB documentation for \code{\link{MakeADFun}}.
#'
#' @export
#' @rdname setModel4M
setMethod(f = "setModel4M",
          signature = "missing",
          definition = function(N.States = 2,
                                Use.HMM = FALSE,
                                Distribution = "gamma",
                                Zero.Inflation = c(Step.Length = FALSE,
                                                   Deflection.Angle = FALSE),
                                Starting.Values = list(),
                                Parameter.Mapping = list()) {
               Set.Model.4M<- new("SetModel4M",
                                  N.States = N.States,
                                  Use.HMM = Use.HMM,
                                  Distribution = Distribution,
                                  Zero.Inflation = Zero.Inflation,
                                  Starting.Values = Starting.Values,
                                  Parameter.Mapping = Parameter.Mapping)
                return(Set.Model.4M)                   
                }
)
#' @export
#' @rdname setModel4M
setMethod(f = "setModel4M",
          signature = "Model4M",
          definition = function(x) return(as(x,"SetModel4M"))
)

# Maybe I don't need to document this? roxygen2 complains with ``missing name".
setValidity(
        Class = "SetModel4M",
        method = function(object) {
                Errors<- NULL
                
                if( length(nStates(object)) != 1 ) {
                        Errors<- c(Errors,"N.States must have length 1.")
                } else {}
                if( length(useHMM(object)) != 1 ) {
                        Errors<- c(Errors,"Use.HMM must have length 1.")
                } else {}
                if( length(distribution(object)) != 1 ) {
                        Errors<- c(Errors,"Distribution must have length 1.")
                } else {}
                if( length(zeroInflation(object)) != 2 | 
                        !all(c("Step.Length","Deflection.Angle") %in% names(zeroInflation(object))) ) {
                        Errors<- c(Errors,"Zero.Inflation must be of length two with names Step.Length and Deflection.Angle.")
                } else {}
                ###
                ### Starting Values Validity
                ###
                if( length(startingValues(object)) != 6 |
                        !all(c("Tpm.Working.Pars",
                               "Theta.Working.Pars",
                               "Step.Working.Pars",
                               "Logit.Step.Zero.Probs",
                               "Logit.Angle.Zero.Probs",
                               "Angle.Zero.Working.Pars") %in% names(startingValues(object))) ) {
                         Errors<- c(Errors,"Starting.Values must be of length six with names Tpm.Working.Pars, Theta.Working.Pars, Step.Working.Pars, Logit.Step.Zero.Probs, Logit.Angle.Zero.Probs, and Angle.Zero.Working.Pars.")      
                } else if ( dim(startingValues(object)$Tpm.Working.Pars) != c(nStates(object),
                                                                             nStates(object)) |
                            dim(startingValues(object)$Theta.Working.Pars) != c(nStates(object),
                                                                               2) |
                            dim(startingValues(object)$Step.Working.Pars) != c(nStates(object),
                                                                               3) | 
                            length(startingValues(object)$Logit.Step.Zero.Probs) != nStates(object) |
                            length(startingValues(object)$Logit.Angle.Zero.Probs) != nStates(object) |
                            dim(startingValues(object)$Angle.Zero.Working.Pars) != c(nStates(object),
                                                                                     2) ) {
                        Errors<- c(Errors,"Dimensions of starting values cannot be changed. See startingValues(x).")
                }
                
                if( length(Errors) > 0 ) {
                        stop(paste(Errors,collapse = "\n"))
                }
        }
)










#' Create or extract a \code{Model4M} object.
#'
#' @details Consider using \code{fit} to create a new \code{Model4M} object.
#'
#' @param x An object to coerce to class \code{Model4M}, possibly of class \code{Data4M}.
#'
#' @family Construct4M
#' @export
setGeneric(name = "model4M",
           def = function(x,...) standardGeneric("model4M")
)
#' @export
#' @rdname model4M
setMethod(
        f = "initialize",
        signature = "Model4M",
        definition = function(.Object,
                              Data4M = new("Data4M"),
                              SetModel4M = new("SetModel4M")) {
                as(.Object,"Data4M")<- Data4M
                as(.Object,"SetModel4M")<- SetModel4M
                
                parameterEstimates(.Object)<- list(Deflection.Angle.Parameters = 0,
                                                   Step.Length.Parameters = 0,
                                                   Transition.Probability.Matrix = 0,
                                                   Stationary.Distribution = 0,
                                                   Zero.Inflation = list(
                                                        Step.Zero.Probs = 0,
                                                        Angle.Zero.Probs = 0,
                                                        Angle.Zero.Pars = data.frame(
                                                                Center = 0,
                                                                Concentration = 0
                                                        )
                                                   ))
                aic(.Object)<- 0
                residuals(.Object)<- data.frame(
                        Deflection.Angle<- 0,
                        Step.Length<- 0
                )
                viterbiPath(.Object)<- 0L
                convergence(.Object)<- ""
                tmbEnvironment(.Object)<- emptyenv()
                
                return(.Object)
        }
)
#' @param SetModel4M An object of class \code{SetModel4M}.
#'
#' @export
#' @rdname model4M
setMethod(f = "model4M",
          signature = "Data4M",
          definition = function(x,Set.Model) {
        return(new("Model4M",Data4M = x,SetModel4M = Set.Model))
        }
)
#' @export
#' @rdname model4M
setMethod(f = "model4M",
          signature = "Simulate4M",
          definition = function(x) return(as(x,"Model4M"))
)





#' Create or extract a \code{Simulate4M} object.
#'
#' @details Consider using \code{\link[=simulate.4M]{simulate}} to create new objects of class \code{Simulate4M}.
#'
#' @param x An object to coerce to class \code{Simulate4M}, possible of class \code{Model4M}.
#'
#' @family Construct4M
#' @export
setGeneric(name = "simulate4M",
           def = function(x,...) standardGeneric("simulate4M"))
#' @details The \code{initialize} method is not meant to be used by the user. Use \code{simulate4M} instead.
#'
#' @export
#' @rdname simulate4M
setMethod(f = "initialize",
          signature = "Simulate4M",
          definition = function(.Object,
                                Model4M = new("Model4M")) {
                
                as(.Object,"Model4M")<- Model4M
                
                simulatedLocations(.Object)<- array(0,dim = c(1,1,1))
                simulatedData(.Object)<- array(0,dim = c(1,1,1))
                simulatedViterbiPath(.Object)<- array(0,dim = c(1,1))
                refitParameters(.Object)<- list(Deflection.Angle.Parameters = 0,
                                                Step.Length.Parameters = 0,
                                                Transition.Probability.Matrix = 0,
                                                Stationary.Distribution = 0,
                                                Zero.Inflation = list(
                                                        Step.Zero.Probs = 0,
                                                        Angle.Zero.Probs = 0,
                                                        Angle.Zero.Pars = data.frame(
                                                                Center = 0,
                                                                Concentration = 0
                                                        )
                                                )
                )
                
                refitViterbiPath(.Object)<- array(0,dim = c(1,1))
                                                
                refitResiduals(.Object)<- array(0,dim = c(1,1,1))
                
                return(.Object)
          }
)
#' @export
#' @rdname simulate4M
setMethod(f = "simulate4M",
          signature = "Model4M",
          definition = function(x) {
                return(new("Simulate4M",Model4M = x))
          }
)






######################
###                ###
###  AccessData4M  ###
###                ###
######################

#' Get or set slots from an object of class \code{Data4M}.
#'
#' @param x An object of class \code{Data4M}.
#' @param value A replacement value.
#'
#' @family Access4M
#' @name AccessData4M
NULL
setMethod(f = "$",
                        signature = "Data4M",
                        definition = function(x,name) {
                                slotNames<- grep(name,
                                                 slotNames(x),
                                                 ignore.case =T,
                                                 value = T)
                                slots<- lapply(slotNames,function(slot) return(slot(x,slot)))
                                if( length(slots) == 1 ) {
                                        slots<- slots[[1]]
                                } else {
                                        warning("More than one slot accessed, returning list of all values.")
                                        names(slots)<- slotNames
                                }
                                return(slots)
                        }
)

#' @export
#' @rdname AccessData4M
setGeneric(name = "identification",
                        def = function(x) standardGeneric("identification")
)
#' @export
setMethod(f = "identification",
                        signature = "Data4M",
                        definition = function(x) return(x@Identification)
)
#' @export
#' @rdname AccessData4M
setGeneric(name = "identification<-",
                        def = function(x,value) standardGeneric("identification<-")
)
#' @export
setReplaceMethod(f = "identification",
                        signature = "Data4M",
                        definition = function(x,value) {
                           x@Identification<- value
                           return(x)
                        }
)

#' @export
#' @rdname AccessData4M
setGeneric(name = "observedLocations",
                        def = function(x) standardGeneric("observedLocations")
)
#' @export
setMethod(f = "observedLocations",
                        signature = "Data4M",
                        definition = function(x) return(x@Observed.Locations)
)
#' @export
#' @rdname AccessData4M
setGeneric(name = "observedLocations<-",
                        def = function(x,value) standardGeneric("observedLocations<-")
)
#' @export
setReplaceMethod(f = "observedLocations",
                        signature = "Data4M",
                        definition = function(x,value) {
                           x@Observed.Locations<- value
                           return(x)
                        }
)

#' @export
#' @rdname AccessData4M
setGeneric(name = "interpolationParameters",
                        def = function(x) standardGeneric("interpolationParameters")
)
#' @export
setMethod(f = "interpolationParameters",
                        signature = "Data4M",
                        definition = function(x) return(x@Interpolation.Parameters)
)
#' @export
#' @rdname AccessData4M
setGeneric(name = "interpolationParameters<-",
                        def = function(x,value) standardGeneric("interpolationParameters<-")
)
#' @export
setReplaceMethod(f = "interpolationParameters",
                        signature = "Data4M",
                        definition = function(x,value) {
                           x@Interpolation.Parameters<- value
                           return(x)
                        }
)

#' @export
#' @rdname AccessData4M
setGeneric(name = "interpolatedLocations",
                        def = function(x) standardGeneric("interpolatedLocations")
)
#' @export
setMethod(f = "interpolatedLocations",
                        signature = "Data4M",
                        definition = function(x) return(x@Interpolated.Locations)
)
#' @export
#' @rdname AccessData4M
setGeneric(name = "interpolatedLocations<-",
                        def = function(x,value) standardGeneric("interpolatedLocations<-")
)
#' @export
setReplaceMethod(f = "interpolatedLocations",
                        signature = "Data4M",
                        definition = function(x,value) {
                           x@Interpolated.Locations<- value
                           return(x)
                        }
)

#' @export
#' @rdname AccessData4M
setGeneric(name = "movementData",
                        def = function(x) standardGeneric("movementData")
)
#' @export
setMethod(f = "movementData",
                        signature = "Data4M",
                        definition = function(x) return(x@Movement.Data)
)
#' @export
#' @rdname AccessData4M
setGeneric(name = "movementData<-",
                        def = function(x,value) standardGeneric("movementData<-")
)
#' @export
setReplaceMethod(f = "movementData",
                        signature = "Data4M",
                        definition = function(x,value) {
                           x@Movement.Data<- value
                           return(x)
                        }
)

#' @export
#' @rdname AccessData4M
setGeneric(name = "groups",
           def = function(x) standardGeneric("groups")
)
#' @export
setMethod(f = "groups",
          signature = "Data4M",
          definition = function(x) {
                data<- movementData(x)
                groups<- data$Movement.Data$Group
                return(groups)
          }
)




#' Get or set slots of an object of class \code{SetModel4M}.
#'
#' @param x An object of class \code{SetModel4M}.
#' @param value A replacement value.
#'
#' @family Access4M
#' @name AccessSetModel4M
NULL

setMethod(f = "$",
          signature = "SetModel4M",
          definition = function(x,name) {
                slotNames<- grep(name,
                                 slotNames(x),
                                 ignore.case =T,
                                 value = T)
		slots<- lapply(slotNames,
		               function(slot) return(slot(x,slot)))
                if( length(slots) == 1 ) {
                        slots<- slots[[1]]
                } else {
                        warning("More than one slot accessed, returning list of all values.")
                        names(slots)<- slotNames
	        }
                return(slots)
	}
)

#' @export
#' @rdname AccessSetModel4M
setGeneric(name = "nStates",
           def = function(x) standardGeneric("nStates")
)
#' @export
setMethod(f = "nStates",
          signature = "SetModel4M",
          definition = function(x) return(x@N.States)
)
#' @export
#' @rdname AccessSetModel4M
setGeneric(name = "nStates<-",
           def = function(x,value) standardGeneric("nStates<-")
)
#' @export
setReplaceMethod(f = "nStates",
                 signature = "SetModel4M",
                 definition = function(x,value) {
                        x@N.States<- as.integer(value)
                        return(x)
                 }
)


#' @export
#' @rdname AccessSetModel4M
setGeneric(name = "useHMM",
           def = function(x) standardGeneric("useHMM")
)
#' @export
setMethod(f = "useHMM",
          signature = "SetModel4M",
          definition = function(x) return(x@Use.HMM)
)
#' @export
#' @rdname AccessSetModel4M
setGeneric(name = "useHMM<-",
           def = function(x,value) standardGeneric("useHMM<-")
)
#' @export
setReplaceMethod(f = "useHMM",
                 signature = "SetModel4M",
                 definition = function(x,value) {
                        x@Use.HMM<- value
                        return(x)
                 }
)


#' @export
#' @rdname AccessSetModel4M
setGeneric(name = "distribution",
           def = function(x) standardGeneric("distribution")
)
#' @export
setMethod(f = "distribution",
          signature = "SetModel4M",
          definition = function(x) return(x@Distribution)
)
#' @export
#' @rdname AccessSetModel4M
setGeneric(name = "distribution<-",
           def = function(x,value) standardGeneric("distribution<-")
)
#' @export
setReplaceMethod(f = "distribution",
                 signature = "SetModel4M",
                 definition = function(x,value) {
                        x@Distribution<- grep(value,
                                              .implemented.distributions,
                                              ignore.case = T,
                                              value = T)
                        return(x)
                 }
)


#' @export
#' @rdname AccessSetModel4M
setGeneric(name = "zeroInflation",
           def = function(x) standardGeneric("zeroInflation")
)
#' @export
setMethod(f = "zeroInflation",
          signature = "SetModel4M",
          definition = function(x) return(x@Zero.Inflation)
)
#' @export
#' @rdname AccessSetModel4M
setGeneric(name = "zeroInflation<-",
           def = function(x,value) standardGeneric("zeroInflation<-")
)
#' @export
setReplaceMethod(f = "zeroInflation",
                 signature = "SetModel4M",
                 definition = function(x,value) {
                        x@Zero.Inflation<- value
                        return(x)
                 }
)


#' @export
#' @rdname AccessSetModel4M
setGeneric(name = "startingValues",
           def = function(x) standardGeneric("startingValues")
)
#' @export
setMethod(f = "startingValues",
          signature = "SetModel4M",
          definition = function(x) return(x@Starting.Values)
)
#' @export
#' @rdname AccessSetModel4M
setGeneric(name = "startingValues<-",
           def = function(x,value) standardGeneric("startingValues<-")
)
#' @export
setReplaceMethod(f = "startingValues",
                 signature = "SetModel4M",
                 definition = function(x,value) {
                        ###
                        ### Replaces only the supplied elements in value, leaves the rest the same.
                        ###
                        Init.Starting.Values<- startingValues(x)
                        
                        Given.Starting.Value.Names<- lapply(names(Init.Starting.Values), function(name) {
                                Given.Names<- grep(names(value),
                                                   name,
                                                   ignore.case = T,
                                                   value = T)
                                                   
                                if(length(Given.Names) == 0) {
                                        return(NULL)
                                } else {
                                        return(Given.Names)
                                }
                        })
                        
                        Given.Starting.Value.Names<- do.call(c,Given.Starting.Value.Names)
                        if( length(Given.Starting.Value.Names) != length(value) ) {
                                stop("Names of replacement must be in Tpm.Working.Pars, Theta.Working.Pars, Step.Working.Pars, Logit.Step.Zero.Probs, Logit.Angle.Zero.Probs, or Angle.Zero.Working.Pars.") 
                        } else {
                                names(value)<- Given.Starting.Value.Names
                                Init.Starting.Values[names(value)]<- value
                        }
                        
                        x@Starting.Values<- Init.Starting.Values
                        return(x)
                 }
)


#' @export
#' @rdname AccessSetModel4M
setGeneric(name = "parameterMapping",
           def = function(x) standardGeneric("parameterMapping")
)
#' @export
setMethod(f = "parameterMapping",
          signature = "SetModel4M",
          definition = function(x) return(x@Parameter.Mapping)
)
#' @export
#' @rdname AccessSetModel4M
setGeneric(name = "parameterMapping<-",
           def = function(x,value) standardGeneric("parameterMapping<-")
)
#' @export
setReplaceMethod(f = "parameterMapping",
                 signature = "SetModel4M",
                 definition = function(x,value) {
                        ###
                        ### Replaces only the supplied elements in value, leaves the rest the same.
                        ###
                        Init.Parameter.Mapping<- parameterMapping(x)
                        
                        Given.Parameter.Names<- lapply(names(Init.Parameter.Mapping),function(name) {
                                Given.Names<- grep(names(value),
                                                   name,
                                                   ignore.case = T,
                                                   value = T)
                                
                                if(length(Given.Names) == 0) {
                                        return(NULL)
                                } else {
                                        return(Given.Names)
                                }
                        })
                        
                        if( length(Given.Parameter.Names) != length(value) ) {
                                stop("Names of replacement must be in Tpm.Working.Pars, Theta.Working.Pars, Step.Working.Pars, Logit.Step.Zero.Probs, Logit.Angle.Zero.Probs, or Angle.Zero.Working.Pars.") 
                        } else {
                                names(value)<- Given.Parameter.Names
                                Init.Parameter.Mapping[names(value)]<- value
                        }
                        
                        x@Parameter.Mapping<- Init.Parameter.Mapping
                        return(x)
                 }
)





#' Get or set slots of an object of class \code{Model4M}.
#'
#' @param x/object An object of class \code{Model4M}.
#' @param value A replacement value.
#'
#' @family Access4M
#' @name AccessModel4M
NULL

#' @export
#' @rdname AccessModel4M
setGeneric(name = "SetModel4M",
           def = function(x) standardGeneric("SetModel4M")
)
#' @export
setMethod(f = "SetModel4M",
          signature = "Model4M",
          definition = function(x) return(as(x,"SetModel4M"))
)

#' @export
#' @rdname AccessModel4M
setGeneric(name = "parameterEstimates",
           def = function(x) standardGeneric("parameterEstimates")
)
#' @export
setMethod(f = "parameterEstimates",
          signature = "Model4M",
          definition = function(x) return(x@Parameter.Estimates)
)
#' @export
#' @rdname AccessModel4M
setGeneric(name = "parameterEstimates<-",
           def = function(x,value) standardGeneric("parameterEstimates<-")
)
#' @export
setReplaceMethod(f = "parameterEstimates",
          signature = "Model4M",
          definition = function(x,value) {
                x@Parameter.Estimates<- value
                return(x)
          }
)


#' @export
#' @rdname AccessModel4M
setGeneric(name = "aic",
           def = function(x) standardGeneric("aic")
)
#' @export
setMethod(f = "aic",
          signature = "Model4M",
          definition = function(x) return(x@AIC)
)
#' @export
#' @rdname AccessModel4M
setGeneric(name = "aic<-",
           def = function(x,value) standardGeneric("aic<-")
)
#' @export
setReplaceMethod(f = "aic",
          signature = "Model4M",
          definition = function(x,value) {
                x@AIC<- value
                return(x)
          }
)

#setGeneric(name = "residuals",
#           def = function(x) standardGeneric("residuals")
#)
#' @export
#' @rdname AccessModel4M
setMethod(f = "residuals",
          signature = "Model4M",
          definition = function(object) return(object@Residuals)
)
#' @export
#' @rdname AccessModel4M
setGeneric(name = "residuals<-",
           def = function(x,value) standardGeneric("residuals<-")
)
#' @export
setReplaceMethod(f = "residuals",
          signature = "Model4M",
          definition = function(x,value) {
                x@Residuals<- value
                return(x)
          }
)


#' @export
setGeneric(name = "viterbiPath",
           def = function(x) standardGeneric("viterbiPath")
)
#' @export
setMethod(f = "viterbiPath",
          signature = "Model4M",
          definition = function(x) return(x@Viterbi.Path)
)
#' @export
#' @rdname AccessModel4M
setGeneric(name = "viterbiPath<-",
           def = function(x,value) standardGeneric("viterbiPath<-")
)
#' @export
setReplaceMethod(f = "viterbiPath",
          signature = "Model4M",
          definition = function(x,value) {
                x@Viterbi.Path<- value
                return(x)
          }
)


#' @export
#' @rdname AccessModel4M
setGeneric(name = "convergence",
           def = function(x) standardGeneric("convergence")
)
#' @export
setMethod(f = "convergence",
          signature = "Model4M",
          definition = function(x) return(x@Convergence)
)
#' @export
#' @rdname AccessModel4M
setGeneric(name = "convergence<-",
           def = function(x,value) standardGeneric("convergence<-")
)
#' @export
setReplaceMethod(f = "convergence",
          signature = "Model4M",
          definition = function(x,value) {
                x@Convergence<- value
                return(x)
          }
)


#' @export
#' @rdname AccessModel4M
setGeneric(name = "tmbEnvironment",
           def = function(x) standardGeneric("tmbEnvironment")
)
#' @export
setMethod(f = "tmbEnvironment",
          signature = "Model4M",
          definition = function(x) return(x@TmbEnvironment)
)
#' @export
#' @rdname AccessModel4M
setGeneric(name = "tmbEnvironment<-",
           def = function(x,value) standardGeneric("tmbEnvironment<-")
)
#' @export
setReplaceMethod(f = "tmbEnvironment",
          signature = "Model4M",
          definition = function(x,value) {
                x@TmbEnvironment<- value
                return(x)
          }
)









#' Get or set slots of an object of class \code{Simulate4M}.
#'
#' @param x An object of class \code{Simulate4M}.
#' @param value A replacement value.
#'
#' @family Access4M
#' @name AccessSimulate4M
NULL

#' @export
#' @rdname AccessSimulate4M
setGeneric(name = "simulatedLocations",
           def = function(x) standardGeneric("simulatedLocations")
)
#' @export
setMethod(f = "simulatedLocations",
          signature = "Simulate4M",
          def = function(x) return(x@Simulated.Locations)
)
#' @export
#' @rdname AccessSimulate4M
setGeneric(name = "simulatedLocations<-",
           def = function(x,value) standardGeneric("simulatedLocations<-")
)
#' @export
setReplaceMethod(f = "simulatedLocations",
          signature = "Simulate4M",
          def = function(x,value) {
                x@Simulated.Locations<- value
                return(x)
          }
)

#' @export
#' @rdname AccessSimulate4M
setGeneric(name = "simulatedData",
           def = function(x) standardGeneric("simulatedData")
)
#' @export
setMethod(f = "simulatedData",
          signature = "Simulate4M",
          definition = function(x) return(x@Simulated.Data)
)
#' @export
#' @rdname AccessSimulate4M
setGeneric(name = "simulatedData<-",
           def = function(x,value) standardGeneric("simulatedData<-")
)
#' @export
setReplaceMethod(f = "simulatedData",
          signature = "Simulate4M",
          definition = function(x,value) {
                x@Simulated.Data<- value
                return(x)
          }
)


#' @export
#' @rdname AccessSimulate4M
setGeneric(name = "simulatedViterbiPath",
           def = function(x) standardGeneric("simulatedViterbiPath")
)
#' @export
setMethod(f = "simulatedViterbiPath",
          signature = "Simulate4M",
          definition = function(x) return(x@Simulated.Viterbi.Path)
)
#' @export
#' @rdname AccessSimulate4M
setGeneric(name = "simulatedViterbiPath<-",
           def = function(x,value) standardGeneric("simulatedViterbiPath<-")
)
#' @export
setReplaceMethod(f = "simulatedViterbiPath",
          signature = "Simulate4M",
          definition = function(x,value) {
                x@Simulated.Viterbi.Path<- value
                return(x)
          }
)


#' @export
#' @rdname AccessSimulate4M
setGeneric(name = "refitParameters",
           def = function(x) standardGeneric("refitParameters")
)
#' @export
setMethod(f = "refitParameters",
          signature = "Simulate4M",
          definition = function(x) return(x@Refit.Parameters)
)
#' @export
#' @rdname AccessSimulate4M
setGeneric(name = "refitParameters<-",
           def = function(x,value) standardGeneric("refitParameters<-")
)
#' @export
setReplaceMethod(f = "refitParameters",
          signature = "Simulate4M",
          definition = function(x,value) {
                x@Refit.Parameters<- value
                return(x)
          }
)


#' @export
#' @rdname AccessSimulate4M
setGeneric(name = "refitResiduals",
           def = function(x) standardGeneric("refitResiduals")
)
#' @export
setMethod(f = "refitResiduals",
          signature = "Simulate4M",
          definition = function(x) return(x@Refit.Residuals)
)
#' @export
#' @rdname AccessSimulate4M
setGeneric(name = "refitResiduals<-",
           def = function(x,value) standardGeneric("refitResiduals<-")
)
#' @export
setReplaceMethod(f = "refitResiduals",
          signature = "Simulate4M",
          definition = function(x,value) {
                x@Refit.Residuals<- value
                return(x)
          }
)


#' @export
#' @rdname AccessSimulate4M
setGeneric(name = "refitViterbiPath",
           def = function(x) standardGeneric("refitViterbiPath")
)
#' @export
setMethod(f = "refitViterbiPath",
          signature = "Simulate4M",
          definition = function(x) return(x@Refit.Viterbi.Path)
)
#' @export
#' @rdname AccessSimulate4M
setGeneric(name = "refitViterbiPath<-",
           def = function(x,value) standardGeneric("refitViterbiPath<-")
)
#' @export
setReplaceMethod(f = "refitViterbiPath",
                 signature = "Simulate4M",
                 definition = function(x,value) {
                        x@Refit.Viterbi.Path<- value
                        return(x)
                 }
)


#' @export
#' @rdname AccessSimulate4M
setGeneric(name = "refitConvergence",
           def = function(x) standardGeneric("refitConvergence")
)
#' @export
setMethod(f = "refitConvergence",
          signature = "Simulate4M",
          definition = function(x) return(x@Refit.Convergence)
)
#' @export
#' @rdname AccessSimulate4M
setGeneric(name = "refitConvergence<-",
           def = function(x,value) standardGeneric("refitConvergence<-")
)
#' @export
setReplaceMethod(f = "refitConvergence",
                 signature = "Simulate4M",
                 definition = function(x,value) {
                        x@Refit.Convergence<- value
                        return(x)
                 }
)


#' @export
#' @rdname AccessSimulate4M
setGeneric(name = "refitEnvironment",
           def = function(x) standardGeneric("refitEnvironment")
)
#' @export
setMethod(f = "refitEnvironment",
          signature = "Simulate4M",
          definition = function(x) return(x@Refit.Environment)
)
#' @export
#' @rdname AccessSimulate4M
setGeneric(name = "refitEnvironment<-",
           def = function(x,value) standardGeneric("refitEnvironment<-")
)
#' @export
setReplaceMethod(f = "refitEnvironment",
                 signature = "Simulate4M",
                 definition = function(x,value) {
                        x@Refit.Environment<- value
                        return(x)
                 }
)





#' Subset or combine objects of class \code{Simulate4M}.
#'
#' @param x An object of class \code{Simulate4M}.
#' @param i Simulations to extract.
#'
#' @name SubsetCombineSimulate4M
NULL

#' @export
#' @rdname SubsetCombineSimulate4M
setMethod(f = "[",
          signature = c(x = "Simulate4M",
                        i = "numeric",
#                        j = "missing",
                        drop = "missing"),
          definition = function(x,i,j,drop) {
                
                if( length(i) == 0 || (length(i) == 1 && i == 0) ) {
                        return(x)
                } else {}
          
                simulatedLocations(x)<- simulatedLocations(x)[,,i]
                simulatedData(x)<- simulatedData(x)[,,i]
                simulatedViterbiPath(x)<- simulatedViterbiPath(x)[i,]
                refitParameters(x)$Deflection.Angle.Parameters<- refitParameters(x)$Deflection.Angle.Parameters[,,i]
                refitParameters(x)$Step.Length.Parameters<- refitParameters(x)$Step.Length.Parameters[,,i]
                refitParameters(x)$Transition.Probability.Matrix<- refitParameters(x)$Transition.Probability.Matrix[,,i]
                refitParameters(x)$Stationary.Distribution<- refitParameters(x)$Stationary.Distribution[i,]
                refitParameters(x)$Zero.Inflation$Step.Zero.Probs<- refitParameters(x)$Zero.Inflation$Step.Zero.Probs[i,]
                refitParameters(x)$Zero.Inflation$Angle.Zero.Probs<- refitParameters(x)$Zero.Inflation$Angle.Zero.Probs[i,]
                refitParameters(x)$Zero.Inflation$Angle.Zero.Pars<- refitParameters(x)$Zero.Inflation$Angle.Zero.Pars[,,i]
                refitResiduals(x)<- refitResiduals(x)[,,i]
                refitViterbiPath(x)<- refitViterbiPath(x)[i,]
                refitConvergence(x)<- refitConvergence(x)[i]
                refitEnvironment(x)<- refitEnvironment(x)[i]
                
                return(x)
          }
)

#' @export
#' @rdname SubsetCombineSimulate4M
setMethod(f = "[[",
          signature = c(x = "Simulate4M",
                        i = "numeric",
                        j = "missing"),
          definition = function(x,i,j) {
                if( length(i) != 1 ) {
                        stop("Index i must be a single value.")
                } else {}
                
                new.Model4M<- as(x,"Model4M")
                
                observedLocations(new.Model4M)<- cbind(Date = interpolationParameters(x)[["Time.Step"]]*seq(nrow(simulatedLocations(x)[,,i])),
                                                       data.frame(simulatedLocations(x)[,,i]))
                observedLocations(new.Model4M)$Date<- as.POSIXct(observedLocations(new.Model4M)$Date,
                                                                 origin = observedLocations(x)$Date[[1]])
                                                                 
                interpolatedLocations(new.Model4M)<- data.frame(simulatedLocations(x)[,,i])
                movementData(new.Model4M)$Movement.Data<- data.frame(simulatedData(x)[,,i])
                
                interpolatedLocations(new.Model4M)$Group<- as.factor(interpolatedLocations(new.Model4M)$Group)
                movementData(new.Model4M)$Movement.Data$Group<- as.factor(movementData(new.Model4M)$Movement.Data$Group)
                
                Parameter.Estimates<- parameterEstimates(new.Model4M)
                Parameter.Estimates$Deflection.Angle.Parameters<- refitParameters(x)$Deflection.Angle.Parameters[,,i]
                Parameter.Estimates$Step.Length.Parameters<- refitParameters(x)$Step.Length.Parameters[,,i]
                Parameter.Estimates$Transition.Probability.Matrix<- refitParameters(x)$Transition.Probability.Matrix[,,i]
                Parameter.Estimates$Zero.Inflation$Step.Zero.Probs<- refitParameters(x)$Zero.Inflation$Step.Zero.Probs[i,]
                Parameter.Estimates$Zero.Inflation$Angle.Zero.Probs<- refitParameters(x)$Zero.Inflation$Angle.Zero.Probs[i,]
                Parameter.Estimates$Zero.Inflation$Angle.Zero.Pars<- as.data.frame(refitParameters(x)$Zero.Inflation$Angle.Zero.Pars[,,i])
                
                parameterEstimates(new.Model4M)<- Parameter.Estimates
                aic(new.Model4M)<- 0
                residuals(new.Model4M)<- data.frame(refitResiduals(x)[,,i])
                viterbiPath(new.Model4M)<- as.integer(refitViterbiPath(x)[i,])
                convergence(new.Model4M)<- refitConvergence(x)[[i]]
                tmbEnvironment(new.Model4M)<- refitEnvironment(x)[[i]]
                
                return(new.Model4M)
          }
)

#' @param ... Additional objects of class Simulate4M with the same structure as \code{x}.
#'
#' @export
#' @rdname SubsetCombineSimulate4M
setMethod(f = "c",
          signature = "Simulate4M",
          definition = function(x,...) {
                
                if( !requireNamespace("abind",quietly = T) ) {
                        stop("Package abind is required to combine, please install it.")
                }
                
                others<- list(...,x)
                classes<- do.call(c,lapply(others,class))
                if( !all(classes == "Simulate4M") ) {
                        stop("All arguments must have class Simulate4M.")
                } else {}
                
                Simulated.Locations<- lapply(others,simulatedLocations)
                simulatedLocations(x)<- do.call(abind::abind,
                                                c(Simulated.Locations,
                                                  list(along = 3)))
                                                  
                Simulated.Data<- lapply(others,simulatedData)
                simulatedData(x)<- do.call(abind::abind,
                                           c(Simulated.Data,
                                             list(along = 3)))
                                             
                Simulated.Viterbi.Path<- lapply(others,simulatedViterbiPath)
                simulatedViterbiPath(x)<- do.call(rbind,
                                                  Simulated.Viterbi.Path)
                                                  
                Angle.Pars<- lapply(others,function(y) return(refitParameters(y)$Deflection.Angle.Parameters))
                refitParameters(x)$Deflection.Angle.Parameters<- do.call(abind::abind,
                                                                         c(Angle.Pars,
                                                                           list(along = 3)))
                                                                           
                Step.Pars<- lapply(others,function(y) return(refitParameters(y)$Step.Length.Parameters))
                refitParameters(x)$Step.Length.Parameters<- do.call(abind::abind,
                                                                    c(Step.Pars,
                                                                      list(along = 3)))
                                                                      
                TPM<- lapply(others,function(y) return(refitParameters(y)$Transition.Probability.Matrix))
                refitParameters(x)$Transition.Probability.Matrix<- do.call(abind::abind,
                                                                           c(TPM,
                                                                             list(along = 3)))
                                                                             
                Stationaries<- lapply(others,function(y) return(refitParameters(y)$Stationary.Distribution))
                refitParameters(x)$Stationary.Distribution<- do.call(rbind,
                                                                     Stationaries)
                                                                     
                Step.Zero.Probs<- lapply(others,function(y) return(refitParameters(y)$Zero.Inflation$Step.Zero.Probs))
                refitParameters(x)$Zero.Inflation$Step.Zero.Probs<- do.call(rbind,
                                                                            Step.Zero.Probs)
                                                                            
                Angle.Zero.Probs<- lapply(others,function(y) return(refitParameters(y)$Zero.Inflation$Angle.Zero.Probs))
                refitParameters(x)$Zero.Inflation$Angle.Zero.Probs<- do.call(rbind,
                                                                             Angle.Zero.Probs)
                                                                             
                Angle.Zero.Pars<- lapply(others,function(y) return(refitParameters(y)$Zero.Inflation$Angle.Zero.Pars))
                refitParameters(x)$Zero.Inflation$Angle.Zero.Pars<- do.call(abind::abind,
                                                                            c(Angle.Zero.Pars,
                                                                              list(along = 3)))
                                                                              
                Refit.Residuals<- lapply(others,refitResiduals)
                refitResiduals(x)<- do.call(abind::abind,
                                            c(Refit.Residuals,
                                              list(along = 3)))
                                              
                Refit.Viterbi.Path<- lapply(others,refitViterbiPath)
                refitViterbiPath(x)<- do.call(rbind,
                                              Refit.Viterbi.Path)
                                              
                Refit.Convergence<- lapply(others,refitConvergence)
                refitConvergence(x)<- do.call(c,
                                              Refit.Convergence)
                                              
                Refit.Environment<- lapply(others,refitEnvironment)
                refitEnvironment(x)<- do.call(c,
                                              Refit.Environment)
                                              
                return(x)
          }
)




