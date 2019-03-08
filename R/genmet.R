#' @include classes.R getset.R
NULL

#' @export
#' @noRd
setGeneric(name = "difftime")

#' Compute time differences in the observed locations of an object of class \code{Data4M}.
#'
#' @param time1 An object of class \code{Data4M}.
#' @param units Units in which the result is desired. See \code{\link{difftime}}.
#'
#' @return Returns a numeric vector of time differences between observed locations.
#'
#' @examples
#' sealdata<- data4M(greyseal)
#' summary(difftime(sealdata))
#'
#' @export
#' @name difftime4M
setMethod(f = "difftime",
          signature = "Data4M",
          definition = function(time1,
                                units = "hours") {
                return(as.numeric(difftime(time1 = tail(observedLocations(time1)$Date,n = -1),
                                           time2 = head(observedLocations(time1)$Date,n = -1),
                                           units = units)))
          }
)

#' @export
#' @noRd
setGeneric(name = "interpolateSummary",
           def = function(object,...) standardGeneric("interpolateSummary")
)
#' @export
#' @noRd
setMethod(f = "interpolateSummary",
          signature = "Data4M",
          definition = function(object,
                                Time.Step,
                                Group.Cutoff = 2*Time.Step,
                                units = "hrs") {

                units<- grep(units,
                             c("days","hrs","mins","secs"),
                             ignore.case = T,
                             value = T)
                if( length(units) != 1 ) {
                        stop("Units must be one of days, hrs, mins, or secs.")
                }

                invisible(Group.Cutoff) ### Forces lazy evaluation

                if( units == "days" ) {
                        Time.Step<- Time.Step * 24*60*60
                        Group.Cutoff<- Group.Cutoff * 24*60*60
                } else if( units == "hrs" ) {
                        Time.Step<- Time.Step * 60*60
                        Group.Cutoff<- Group.Cutoff * 60*60
                } else if( units == "mins") {
                        Time.Step<- Time.Step * 60
                        Group.Cutoff<- Group.Cutoff * 60
                } else {}

                interpolationParameters(object)<- c("Time.Step" = Time.Step,
                                                    "Group.Cutoff" = Group.Cutoff)

                if( Time.Step <= 0 ) {
                        stop("Time.Step must be positive.")
                } else if ( Group.Cutoff < Time.Step ) {
                        stop("Group.Cutoff must be greater than or equal to Time.Step.")
                } else {}

                Observed.Locations<- observedLocations(object)
                Observed.Locations<- Observed.Locations[!duplicated(Observed.Locations$Date),]

                Group.Streaks<- rle(diff(as.numeric(Observed.Locations$Date))>Group.Cutoff)
                Group.Lengths<- c(1, Group.Streaks$lengths)
                Group.Starts<- c(TRUE,Group.Streaks$values)
                if( !tail(Group.Starts,1) ) {
                        Group.Starts<- c(Group.Starts,TRUE)
                }
                Group.Starts<- which(Group.Starts)  ### Picks out the starting times for each contiguous block

                Group.Indices<- sapply(seq_along(head(Group.Starts,-1)),
                                       function(i) {   ### Creates subsetting indices
                                                return(c(cumsum(Group.Lengths)[Group.Starts[[i]]],
                                                         cumsum(Group.Lengths)[Group.Starts[[i+1]]-1]))
                                       }
                )
                ###
                ### The first row of Group.Indices is the starting time for each group
                ### second The row is the ending time for each group
                ###

                Observed.Locations$Group<- NA
                for(i in seq(ncol(Group.Indices)) ) {
                        Observed.Locations$Group[ Group.Indices[1,i]:Group.Indices[2,i] ]<- i
                }

               Observed.Locations<- Observed.Locations[!(is.na(Observed.Locations$Group)),]
               Observed.Locations$Group<- as.factor(Observed.Locations$Group)


               Location.Groups<- split(Observed.Locations,Observed.Locations$Group)

               Interpolated.Groups<- lapply(Location.Groups,
                                            function(group) {
                                                if( nrow(group) < 2 ) {
                                                        return(NULL)
                                                } else {}

                                                t0<- as.numeric(head(group$Date,1))
                                                tT<- as.numeric(tail(group$Date,1))
                                                t<- seq(t0, tT, by = Time.Step)

                                                Interpolated.Locations<- data.frame(Date = as.POSIXct(t,origin="1970-01-01"))
                                                if( nrow(Interpolated.Locations) < 3 ) {
                                                        return(NULL)
                                                } else {}

                                                return(Interpolated.Locations)
                                            }
                                     )

                Null.Indices<- unlist(lapply(Interpolated.Groups,is.null))
                Interpolated.Groups<- Interpolated.Groups[!Null.Indices]
                names(Interpolated.Groups)<- seq_along(Interpolated.Groups)
                Interpolated.Groups<- lapply(names(Interpolated.Groups),function(name) {
                                                group<- Interpolated.Groups[[name]]
                                                group$Group<- as.integer(name)
                                                return(group)
                                      }
                )

                Interpolated.Locations<- do.call(rbind.data.frame,Interpolated.Groups)
                Interpolated.Locations$Group<- as.factor(Interpolated.Locations$Group)

                Interpolation.Summary<- c(Nobs = nrow(Interpolated.Locations),
                                         NGroup = length(unique(Interpolated.Locations$Group)))

                return(Interpolation.Summary)
    }
)

#' @title Compute diagnostic statistics for a grid of time steps and group cutoffs.
#'
#' @export
#' @rdname findStep
setGeneric(name = "findStep",
           def = function(object,...) standardGeneric("findStep")
)

#' @param object An object of or inheriting from class \code{Data4M}.
#' @param step.grid A numeric vector giving the values of Time.Step to consider.
#' @param cutoff.grid A numeric vector giving the multiples of Time.Step to use as the Group.Cutoff (all entries should be between 1 and 2).
#'
#' @return A list giving the values of n_prop and n_adj, along with the values which maximize each statistic.
#'
#' @examples
#' sealdata<- data4M(greyseal)
#' findStep(sealdata,cutoff.grid = seq(1,2,length.out = 5)
#'
#' @export
#' @rdname findStep
setMethod(f = "findStep",
          signature = "Data4M",
          definition = function(object,
                                step.grid = seq(summary(difftime(object))[[3]],
                                                summary(difftime(object))[[5]],
                                                length.out = 10),
                                cutoff.grid = seq(1,2,length.out = 11),
                                units = "hrs") {
                nAdj<- nProp<- matrix(0, nrow = length(step.grid), ncol = length(cutoff.grid))
                rownames(nAdj)<- rownames(nProp)<- step.grid
                colnames(nAdj)<- colnames(nProp)<- cutoff.grid

                for( X in step.grid ) {
                    for( Y in cutoff.grid ) {
                        interp.Summ<- interpolateSummary(object = object,
                                                        Time.Step = X,
                                                        Group.Cutoff = X*Y,
                                                        units = units)

                        nProp[paste(X),paste(Y)]<- interp.Summ[["Nobs"]] / nrow(observedLocations(object))
                        nAdj[paste(X),paste(Y)]<- (interp.Summ[["Nobs"]] - 2*interp.Summ[["NGroup"]]) / (nrow(observedLocations(object)) - 2)
                    }
                }

                whichNProp<- which.max(nProp)
                whichNAdj<- which.max(nAdj)

                nProp.max<- nProp[whichNProp %% nrow(nProp), whichNProp %/% nrow(nProp) + 1,drop = FALSE]
                nAdj.max<- nAdj[whichNAdj %% nrow(nAdj), whichNAdj %/% nrow(nAdj) + 1,drop = FALSE]

                return(list(nProp,nProp.max,nAdj,nAdj.max))
          }
)

#' @title Estimate the autocorrelation function before interpolating locations.
#'
#' @param object An object of or inheriting from class \code{Data4M}.
#'
#' @export
#' @rdname preACF
setGeneric(name = "preACF",
           def = function(object,...) standardGeneric("preACF")
)

#' @description Estimate the acf for step length before interpolating the data.
#'
#' @param Time.Step The time interval used for the regular time grid.
#' @param units Can be one of \code{days}, \code{hrs}, \code{mins}, or \code{secs}.
#' @param lag.max The maximum lag to compute.
#' @param kernel.type Which kernel to use? Can be one of "rectangle", "sinc", or "gaussian". Can also supply a user-defined function that takes the arguments d, kernel.width, and n.
#' @param kernel.width The bandwidth of the kernel used in weighting pairs of observations. Note that the time differences are standardized to have mean 1. The default setting depends on the kernel used.
#'
#' @return An object of class \code{acf}, which is a list containing the following elements:
#'  \itemize{
#'      \item lag: A vector containing the lags at which the acf is estimated.
#'      \item acf: An vector with the same length as 'lag' containing the estimated acf.
#'      \item type: Always set to  "correlation".
#'      \item n.used A vector giving the sum of kernel values for each lag.
#'      \item series The name of the series 'object'.
#'      \item snames Always set to NULL.
#'  }
#'
#' @export
#' @rdname preACF
setMethod(f = "preACF",
          signature = "Data4M",
          definition = function(object,
                                Time.Step,
                                units = "hrs",
                                lag.max = 30,
                                kernel = "gaussian",
                                kernel.width = "default") {

            units<- grep(units,
                         c("days","hrs","mins","secs"),
                         ignore.case = T,
                         value = T)
            if( length(units) != 1 ) {
                    stop("Units must be one of days, hrs, mins, or secs.")
            }

            if( units == "days" ) {
                    Time.Step<- Time.Step * 24*60*60
            } else if( units == "hrs" ) {
                    Time.Step<- Time.Step * 60*60
            } else if( units == "mins") {
                    Time.Step<- Time.Step * 60
            } else {}

            Locations<- observedLocations(object)
            Step.Length<- sapply(seq(nrow(Locations)-1),
                                  function(t) .distance.formula(x.last = Locations[t,],
                                                                x.curr = Locations[t+1,],
                                                                lat = 3,
                                                                lon = 2) / diff(as.numeric(Locations[c(t,t+1),1]))
                            )
            Step.Standard<- (Step.Length - mean(Step.Length))/sd(Step.Length)

            Date<- sapply(seq(nrow(Locations)-1),
                            function(t) mean(Locations[c(t,t+1),1])
            )
            Date<- as.numeric(Date)
            Time.Mean<- mean(diff(Date))

            Date<- Date / Time.Mean
            Time.Step<- Time.Step/Time.Mean

            if( identical(kernel,"gaussian") ) {
                kernel<- function(d,kernel.width,n) dnorm(x = d,mean = 0,sd = kernel.width)
                if( kernel.width == "default" ) {
                    kernel.width<- 0.25
                }
            } else if( identical(kernel,"sinc") ) {
                kernel<- function(d,kernel.width,n) ifelse(d != 0,sin(pi*kernel.width*d)/(pi*kernel.width*d*n),1)
                if( kernel.width == "default" ) {
                    kernel.width<- 1
                }
            } else if( identical(kernel,"rectangle") ) {
                kernel<- function(d,kernel.width,n) 1*(abs(d)<= 0.5*kernel.width)
                if( kernel.width == "default" ) {
                    kernel.width<- 0.5
                }
            } else if( class(kernel) != "function" || any(!(c("d","kernel.width,n") %in% names(formals(kernel)))) ) {
                stop("kernel must be either a function taking arguments d, kernel.width and n, or one of rectangle, sinc, or gaussian.")
            } else {}

            Acf.Kern<- Acf.Vals<- lag<- array(seq(lag.max+1)-1,dim = c(lag.max+1,1,1))
            Indices<- expand.grid(seq_along(Step.Standard),seq_along(Step.Standard))
            Obs.Lags<- data.frame(Step1 = Step.Standard[Indices[,1]],
                                  Step2 = Step.Standard[Indices[,2]],
                                  Time.Diff = abs(Date[Indices[,1]] - Date[Indices[,2]]),
                                  Kernel.Vals = 0)

            for( k in seq_along(lag) ) {
                Obs.Lags$Kernel.Vals<- kernel(abs( lag[[k]]*Time.Step - Obs.Lags$Time.Diff ),kernel.width = kernel.width,n = nrow(Obs.Lags))
                Acf.Kern[[k]]<- sum(Obs.Lags$Kernel.Vals)
                Acf.Vals[[k]]<- with(Obs.Lags,sum(Step1*Step2*Kernel.Vals)/Acf.Kern[[k]])
            }

            Acf.Object<- structure(list(lag = lag,
                                        acf = Acf.Vals,
                                        type = "correlation",
                                        n.used = Acf.Kern,
                                        series = "stepLength",
                                        snames = NULL),
                                   class = "acf")

            return(Acf.Object)
        }
)


#' @title Estimate the autocorrelation function before interpolating locations.
#'
#' @param object An object of or inheriting from class \code{Data4M}.
#'
#' @export
#' @rdname postACF
setGeneric(name = "postACF",
           def = function(object,...) standardGeneric("postACF")
)

#' @description Estimate the acf for step length after interpolating the data.
#'
#' @param lag.max The maximum lag to compute.
#'
#' @return An object of class \code{acf}, which is a list containing the following elements:
#'  \itemize{
#'      \item lag: A vector containing the lags at which the acf is estimated.
#'      \item acf: An vector with the same length as 'lag' containing the estimated acf.
#'      \item type: Always set to  "correlation".
#'      \item n.used A vector giving the sum of kernel values for each lag.
#'      \item series The name of the series 'object'.
#'      \item snames Always set to NULL.
#'  }
#'
#' @export
#' @rdname postACF
setMethod(f = "postACF",
          signature = "Data4M",
          definition = function(object,
                                lag.max = 30) {

            Movement.Data<- movementData(object)
            Step.Lengths<- Movement.Data$Movement.Data[,c("Group","Step.Length")]
            Starts<- cbind.data.frame(Group = factor(as.numeric(names(Movement.Data$Step.Length.Starting.Values))),
                                      Step.Length = Movement.Data$Step.Length.Starting.Values)

            Step.Lengths<- rbind.data.frame(Starts,Step.Lengths)
            Step.Lengths$Step.Length<- with(Step.Lengths,(Step.Length - mean(Step.Length))/sd(Step.Length))

            Acf.Num<- Acf.Vals<- lag<- array(seq(lag.max+1)-1,dim = c(lag.max+1,1,1))

            for( k in seq_along(lag) ) {
                groupAcf<- tapply(Step.Lengths$Step.Length,Step.Lengths$Group,function(group) {
                    if( length(group) <= lag[[k]] ) {
                        prods<- 0
                        num<- 0
                    } else if( lag[[k]] == 0 ) {
                        prods<- group*group
                        num<- length(group)
                    } else {
                        prods<- head(group,-lag[[k]]) * tail(group,-lag[[k]])
                        num<- length(group)-lag[[k]]
                    }
                    return(data.frame(prod = sum(prods),num = num))
                })
                fullAcf<- colSums(do.call(rbind.data.frame,groupAcf))
                Acf.Vals[[k]]<- fullAcf[[1]] / fullAcf[[2]]
                Acf.Num[[k]]<- fullAcf[[2]]
            }

            Acf.Object<- structure(list(lag = lag,
                                        acf = Acf.Vals,
                                        type = "correlation",
                                        n.used = Acf.Num,
                                        series = "stepLength",
                                        snames = NULL),
                                   class = "acf")

            return(Acf.Object)
        }
)


#' @title Linearly interpolate and group locations to a regular time step.
#'
#' @param object An object of or inheriting from class \code{Data4M}.
#'
#' @export
#' @rdname interpolate
setGeneric(name = "interpolate",
           def = function(object,...) standardGeneric("interpolate")
)

#' @description Provide a time interval and cutoff level to linearly interpolate a movement track.
#'
#' @param method Which interpolation method to use? Can be linear, fmm, or natural (the last two of which are spline methods).
#' @param Time.Step The time interval used in the regular time grid.
#' @param Group.Cutoff Observations which are spaced further apart than this will be separated into different group.
#' @param units Can be one of \code{days}, \code{hrs}, \code{mins}, or \code{secs}.
#'
#' @return Returns a copy of \code{object} with the Interpolation.Parameters, Interpolated.Locations, and Movement.Data slots filled in.
#'
#' @examples
#' sealdata<- data4M(greyseal)
#' sealdata<- interpolate(sealdata,Time.Step = 66, units = "mins")
#'
#' @export
#' @rdname interpolate
setMethod(f = "interpolate",
          signature = "Data4M",
          definition = function(object,
                                method = "linear",
                                Time.Step,
                                Group.Cutoff = 2*Time.Step,
                                units = "hrs") {

                units<- grep(units,
                             c("days","hrs","mins","secs"),
                             ignore.case = T,
                             value = T)
                if( length(units) != 1 ) {
                        stop("Units must be one of days, hrs, mins, or secs.")
                }

                invisible(Group.Cutoff) ### Forces lazy evaluation

                if( units == "days" ) {
                        Time.Step<- Time.Step * 24*60*60
                        Group.Cutoff<- Group.Cutoff * 24*60*60
                } else if( units == "hrs" ) {
                        Time.Step<- Time.Step * 60*60
                        Group.Cutoff<- Group.Cutoff * 60*60
                } else if( units == "mins") {
                        Time.Step<- Time.Step * 60
                        Group.Cutoff<- Group.Cutoff * 60
                } else {}

                interpolationParameters(object)<- c("Time.Step" = Time.Step,
                                                    "Group.Cutoff" = Group.Cutoff)

                if( Time.Step <= 0 ) {
                        stop("Time.Step must be positive.")
                } else if ( Group.Cutoff < Time.Step ) {
                        stop("Group.Cutoff must be greater than or equal to Time.Step.")
                } else {}

                Observed.Locations<- observedLocations(object)
                Observed.Locations<- Observed.Locations[!duplicated(Observed.Locations$Date),]

                Group.Streaks<- rle(diff(as.numeric(Observed.Locations$Date))>Group.Cutoff)
                Group.Lengths<- c(1, Group.Streaks$lengths)
                Group.Starts<- c(TRUE,Group.Streaks$values)
                if( !tail(Group.Starts,1) ) {
                        Group.Starts<- c(Group.Starts,TRUE)
                }
                Group.Starts<- which(Group.Starts)  ### Picks out the starting times for each contiguous block

                Group.Indices<- sapply(seq_along(head(Group.Starts,-1)),
                                       function(i) {   ### Creates subsetting indices
                                                return(c(cumsum(Group.Lengths)[Group.Starts[[i]]],
                                                         cumsum(Group.Lengths)[Group.Starts[[i+1]]-1]))
                                       }
                )
                ###
                ### The first row of Group.Indices is the starting time for each group
                ### second The row is the ending time for each group
                ###

                Observed.Locations$Group<- NA
                for(i in seq(ncol(Group.Indices)) ) {
                        Observed.Locations$Group[ Group.Indices[1,i]:Group.Indices[2,i] ]<- i
                }

               Observed.Locations<- Observed.Locations[!(is.na(Observed.Locations$Group)),]
               Observed.Locations$Group<- as.factor(Observed.Locations$Group)


               Location.Groups<- split(Observed.Locations,Observed.Locations$Group)

               Interpolated.Groups<- lapply(Location.Groups,
                                            function(group) {
                                                if( nrow(group) < 2 ) {
                                                        return(NULL)
                                                } else {}

                                                t0<- as.numeric(head(group$Date,1))
                                                tT<- as.numeric(tail(group$Date,1))
                                                t<- seq(t0, tT, by = Time.Step)

                                                if( method == "linear") {
                                                    Interpolated.Lon<- approx(group$Date,
                                                                              group$Lon,
                                                                              xout = t)$y
                                                    Interpolated.Lat<- approx(group$Date,
                                                                              group$Lat,
                                                                              xout = t)$y
                                                } else if( method %in% c("fmm", "natural") ) {
                                                    Interpolated.Lon<- stats::spline(group$Date,
                                                                                     group$Lon,
                                                                                     xout = t,
                                                                                     method = method)$y
                                                    Interpolated.Lat<- stats::spline(group$Date,
                                                                                     group$Lat,
                                                                                     xout = t,
                                                                                     method = method)$y
                                                } else {
                                                    stop("Method must be either linear, fmm, or natural.")
                                                }

                                                Interpolated.Locations<- data.frame(Date = as.POSIXct(t,origin="1970-01-01"),
                                                                                    Lon = Interpolated.Lon,
                                                                                    Lat = Interpolated.Lat
                                                                                )
                                                if( nrow(Interpolated.Locations) < 3 ) {
                                                        return(NULL)
                                                } else {}

                                                return(Interpolated.Locations)
                                            }
                                     )

                Null.Indices<- unlist(lapply(Interpolated.Groups,is.null))
                Interpolated.Groups<- Interpolated.Groups[!Null.Indices]
                names(Interpolated.Groups)<- seq_along(Interpolated.Groups)
                Interpolated.Groups<- lapply(names(Interpolated.Groups),function(name) {
                                                group<- Interpolated.Groups[[name]]
                                                group$Group<- as.integer(name)
                                                return(group)
                                      }
                )

                Interpolated.Locations<- do.call(rbind.data.frame,Interpolated.Groups)
                Interpolated.Locations$Group<- as.factor(Interpolated.Locations$Group)

                interpolatedLocations(object)<- Interpolated.Locations


                Movement.Data<- by(Interpolated.Locations,
                                   INDICES = Interpolated.Locations$Group,
                                   function(group) {
                                        Deflection.Angle<- sapply(2:(nrow(group)-1),
                                                                  function(t) .deflection.angle.formula(x.last = group[t-1,],
                                                                                                        x.curr = group[t,],
                                                                                                        x.next = group[t+1,],
                                                                                                        lat = 3,
                                                                                                        lon = 2)
                                                                 )

                                       Step.Length<- sapply(1:(nrow(group)-1),
                                                            function(t) .distance.formula(x.last = group[t,],
                                                                                          x.curr = group[t+1,],
                                                                                          lat = 3,
                                                                                          lon = 2))

                                       return(data.frame(Deflection.Angle = c(NA,Deflection.Angle),
                                                         Step.Length = Step.Length,
                                                         Group = group$Group[[1]]))
                                  }
                )

                Step.Length.Mean<- mean(unlist(lapply(Movement.Data,
                                                      function(group) {
                                                        return(group$Step.Length)
                                                      }
                )))
                Step.Length.Starting.Values<- unlist(lapply(Movement.Data,
                                                            function(group) {
                                                                return(group$Step.Length[[1]] / Step.Length.Mean)
                                                            }
                ))

                Movement.Data<- do.call(rbind.data.frame,
                                        lapply(Movement.Data,
                                               function(group) {
                                                return(group[-1,])
                                               })
                )
                Movement.Data[,"Step.Length"]<- Movement.Data[,"Step.Length"] / Step.Length.Mean

                movementData(object) <- list(Movement.Data = Movement.Data,
                                             Step.Length.Starting.Values = Step.Length.Starting.Values,
                                             Step.Length.Mean = Step.Length.Mean)

                return(object)
        }
)


#' Simulate data from an object of class \code{Model4M}.
#'
#' @param object An object of class \code{Model4M}.
#'
#' @name simulate.4M
NULL

#setGeneric(name = "simulate")
#' @param nsim How many simulations should there be?
#' @param seed Sets the RNG seed. If NULL, does not alter the current RNG settings.
#' @param refit Should a model be refit to each set of simulated data?
#' @param keep.locations Should the simulated locations be kept?
#' @param keep.data Should the simulated movement data be kept?
#' @param keep.residuals If \code{refit = T}, should the residuals be kept?
#' @param Use.HMM If \code{refit = T}, should an HMM be fit? Overrides the value stored in \code{model} and \code{Set.Model.4M}.
#' @param Set.Model.4M If \code{refit = T}, which SetModel4M object should be used to refit the data? Use.HMM is controlled separately.
#' @param ... Used to pass arguments to internal functions.
#'
#' @return Returns an object of class \code{Simulate4M} containing simulations and/or results of refitting the simulations.
#'
#' @examples
#' sealData<- data4M(greyseal)
#' sealData<- interpolate(sealData,Time.Step = 1)
#' seal4M2<- fit(sealData)
#' sim4M2<- simulate(seal4M2,nsim = 5)
#'
#' @export
#' @rdname simulate.4M
setMethod(f = "simulate",
          signature = "Model4M",
          definition = function(object,
                                nsim = 1,
                                seed = NULL,
                                refit = T,
                                keep.locations = T,
                                keep.data = T,
                                keep.residuals = T,
                                Use.HMM = useHMM(object),
                                Set.Model.4M = SetModel4M(object),
                                ...) {
                if( !is.null(seed) ) {
                        set.seed(seed)
                }

                scale<- movementData(object)$Step.Length.Mean
                useHMM(Set.Model.4M)<- Use.HMM


                simulations<- as.list(seq(nsim))
                simulations<- lapply(simulations,
                                    function(foo) {
                        Simulated.Data<- tmbEnvironment(object)$simulate()
                        Sim.Data.4M<- new("Data4M")

                        Sim.Angles<- Simulated.Data$sim_theta_series
                        Sim.Angles<- Sim.Angles - 360*floor((Sim.Angles + 180)/360)

                        interpolationParameters(Sim.Data.4M)<- interpolationParameters(object)
                        movementData(Sim.Data.4M)<- list(Movement.Data = data.frame(Deflection.Angle = Sim.Angles,
                                                                                    Step.Length = Simulated.Data$sim_dist_series[-1],
                                                                                    Group = factor(1)),
                                                         Step.Length.Starting.Values = Simulated.Data$sim_dist_series[[1]],
                                                         Step.Length.Mean = movementData(object)$Step.Length.Mean)

                        if( keep.locations == T ) {
                                interpolatedLocations(Sim.Data.4M)<- data.frame(Lon = c(1,2,seq_along(Simulated.Data$sim_state_path)),
                                                                                Lat = 0)


                                interpolatedLocations(Sim.Data.4M)[1:2,]<- interpolatedLocations(object)[1:2,2:3]

                                f.heading<- 180 + .heading.formula(x.last = interpolatedLocations(Sim.Data.4M)[2,],
                                                                   x.curr = interpolatedLocations(Sim.Data.4M)[1,],
                                                                   lat = 2,
                                                                   lon = 1)
                                curr.distance<- movementData(Sim.Data.4M)$Step.Length.Starting.Values * movementData(Sim.Data.4M)$Step.Length.Mean

                                for( k in 3:nrow(interpolatedLocations(Sim.Data.4M)) ) {
                                        interpolatedLocations(Sim.Data.4M)[k,1:2]<- .dest.Point(x.curr = interpolatedLocations(Sim.Data.4M)[k-1,],
                                                                                                f.heading.curr = f.heading,
                                                                                                distance.curr = curr.distance,
                                                                                                deflection.angle = movementData(Sim.Data.4M)$Movement.Data$Deflection.Angle[[k-2]],
                                                                                                distance = movementData(Sim.Data.4M)$Movement.Data$Step.Length[[k-2]] * movementData(Sim.Data.4M)$Step.Length.Mean,
                                                                                                lat = 2,
                                                                                                lon = 1)

                                        f.heading<- 180 + .heading.formula(x.last = interpolatedLocations(Sim.Data.4M)[k,],
                                                                           x.curr = interpolatedLocations(Sim.Data.4M)[k-1,],
                                                                           lat = 2,
                                                                           lon = 1)

                                        curr.distance<- movementData(Sim.Data.4M)$Movement.Data$Step.Length[[k-2]] * movementData(Sim.Data.4M)$Step.Length.Mean
                                }

                                interpolatedLocations(Sim.Data.4M)[,"Group"]<- as.factor(1)
                        } else {}

                        temp.env<- environment()
                        if( refit == T ) {
                                tryCatch({
                                        Sim.Data.4M<- fit(Sim.Data.4M,
                                                          Set.Model.4M,
                                                          #convergence.error = F)
                                                          stationary.tolerance = 10)
                                        assign("Sim.Data.4M",
                                               Sim.Data.4M,
                                               envir = temp.env)
                                        },
                                        error = function(e) {
                                                Sim.Data.4M<- model4M(Sim.Data.4M,
                                                                      Set.Model.4M)
                                                convergence(Sim.Data.4M)<- "ERROR"

                                                assign("Sim.Data.4M",
                                                       Sim.Data.4M,
                                                       envir = temp.env)
                                        })
                        }

                        if( keep.data == F ) {
                                if( refit == F & keep.locations == F ) {
                                        stop("Must re-fit the model, keep the locations, or keep the data.")
                                } else {}

                                movementData(Sim.Data.4M)<- list()
                        }

                        if( keep.residuals == F ) {
                                residuals(Sim.Data.4M)<- data.frame(numeric(0))
                        }

                        return(list(Sim.Data.4M = Sim.Data.4M,
                                    Simulated.Viterbi.Path = as.integer(Simulated.Data$sim_state_path)
                                    )
                              )
                })

                Sim.4M<- simulate4M(object)
                useHMM(Sim.4M)<- Use.HMM


                if( keep.locations == T ) {
                       simulatedLocations(Sim.4M)<- array(0,dim = c(dim(interpolatedLocations(simulations[[1]]$Sim.Data.4M)),nsim),
                                                          dimnames = list(NULL,c("Lat","Lon","Group"),NULL))

                        for( i in seq(nsim) ) {
                                simulatedLocations(Sim.4M)[,,i]<- data.matrix(interpolatedLocations(simulations[[i]]$Sim.Data.4M))
                        }
                } else {
                    simulatedLocations(Sim.4M)<- array(0,dim = c(1,1,nsim))
                }



                if( keep.data == T ) {
                        simulatedData(Sim.4M)<- array(0,dim = c(dim(movementData(simulations[[1]]$Sim.Data.4M)$Movement.Data),nsim),
                                                      dimnames = list(NULL,c("Deflection.Angle","Step.Length","Group"),NULL))

                        for( i in seq(nsim) ) {
                                simulatedData(Sim.4M)[,,i]<- data.matrix(movementData(simulations[[i]]$Sim.Data.4M)$Movement.Data)
                        }
                } else {
                    simulatedData(Sim.4M)<- array(0,dim = c(1,1,nsim))
                }



                simulatedViterbiPath(Sim.4M)<- do.call(rbind,lapply(simulations,function(x) return(x$Simulated.Viterbi.Path))
                )

                Convergence.ERROR<- do.call(c,lapply(simulations,function(i) return(convergence(i$Sim.Data.4M) == "ERROR")))
                correct<- min(which( !Convergence.ERROR ))

                if( refit == T ) {
                        Parameter.Dims<- parameterEstimates(simulations[[correct]]$Sim.Data.4M)
                        Angle.Parameters<- array(0,dim = c(dim(Parameter.Dims$Deflection.Angle.Parameters),nsim))
                        Step.Parameters<- array(0,dim = c(dim(Parameter.Dims$Step.Length.Parameters),nsim))

                        Transition.Probs<- array(0,dim = c(dim(Parameter.Dims$Transition.Probability.Matrix),nsim))
                        Stationary.Distribution<- do.call(rbind,lapply(simulations,function(x) return(parameterEstimates(x$Sim.Data.4M)$Stationary.Distribution)))


                        Step.Zero.Probs<- do.call(rbind,lapply(simulations,function(x) return(parameterEstimates(x$Sim.Data.4M)$Zero.Inflation$Step.Zero.Probs)))
                        Angle.Zero.Probs<- do.call(rbind,lapply(simulations,function(x) return(parameterEstimates(x$Sim.Data.4M)$Zero.Inflation$Angle.Zero.Probs)))
                        Angle.Zero.Pars<- array(0,dim = c(dim(Parameter.Dims$Zero.Inflation$Angle.Zero.Pars),nsim))

                        for( i in which( !Convergence.ERROR) ) {
                                Angle.Parameters[,,i]<- data.matrix(parameterEstimates(simulations[[i]]$Sim.Data.4M)$Deflection.Angle.Parameters)
                                Step.Parameters[,,i]<- data.matrix(parameterEstimates(simulations[[i]]$Sim.Data.4M)$Step.Length.Parameters)
                                Transition.Probs[,,i]<- data.matrix(parameterEstimates(simulations[[i]]$Sim.Data.4M)$Transition.Probability.Matrix)
                                Angle.Zero.Pars[,,i]<- data.matrix(parameterEstimates(simulations[[i]]$Sim.Data.4M)$Zero.Inflation$Angle.Zero.Pars)
                        }

                        dimnames(Angle.Parameters)[[2]]<- c("Center","Concentration")
                        dimnames(Step.Parameters)[[2]]<- c("Reversion.Level","Autocorrelation","Standard.Deviation")
                        dimnames(Angle.Zero.Pars)[[2]]<- c("Center","Concentration")


                        refitParameters(Sim.4M)<- list(Deflection.Angle.Parameters = Angle.Parameters,
                                                       Step.Length.Parameters = Step.Parameters,
                                                       Transition.Probability.Matrix = Transition.Probs,
                                                       Stationary.Distribution = Stationary.Distribution,
                                                       Zero.Inflation = list(Step.Zero.Probs = Step.Zero.Probs,
                                                                             Angle.Zero.Probs = Angle.Zero.Probs,
                                                                             Angle.Zero.Pars = Angle.Zero.Pars)
                                                      )
                        if( keep.residuals == T ) {
                                Refit.Residuals<- array(0,dim = c(dim(residuals(simulations[[correct]]$Sim.Data.4M)),nsim),
                                                        dimnames = list(NULL,c("Deflection.Angle","Step.Length",NULL)))
                                for( i in which( !Convergence.ERROR) ) {
                                        Refit.Residuals[,,i]<- data.matrix(residuals(simulations[[i]]$Sim.Data.4M))
                                }

                                refitResiduals(Sim.4M)<- Refit.Residuals
                        }

                        refitViterbiPath(Sim.4M)<- do.call(rbind,lapply(simulations,function(x) return(viterbiPath(x$Sim.Data.4M))))

                        Refit.Convergence<- rep(NA,nsim)
                        for( i in seq(nsim) ) {
                                Refit.Convergence[[i]]<- convergence(simulations[[i]]$Sim.Data.4M)
                        }
                        refitConvergence(Sim.4M)<- Refit.Convergence

                        Refit.Environment<- list(seq(nsim))
                        for( i in which( !Convergence.ERROR) ) {
                                Refit.Environment[[i]]<- tmbEnvironment(simulations[[i]]$Sim.Data.4M)
                        }
                        refitEnvironment(Sim.4M)<- Refit.Environment
                }

                return(Sim.4M)
          }
)

#' @examples
#'
#' setSim<- setSim4M()
#' parameters(setSim)<- list(Deflect = cbind(c(0,0),c(0.2,0.8)))
#' parameters(setSim)<- list(Step = rbind(c(0.2,0.05,0.3),c(0.02,0.8,0.5)))
#' parameters(setSim)<- list(Trans = rbind(c(0.8,0.2),c(0.3,0.7)))
#' nsim = 1;seed = NULL;refit = T;keep.data=T;keep.residuals=T;DLL="markmodmover"
setMethod(f = "simulate",
          signature = "SetSim4M",
          definition = function(object,
                                nsim = 1,
                                seed = NULL,
                                refit = T,
                                keep.data = T,
                                keep.residuals = T,
                                DLL = "markmodmover",
                                ...) {

                Data.Grouping<- c(0,
                                  cumsum(rle(as.numeric(rep(1,simLength(object))))$lengths)
                                  )

                data<- list(theta = rnorm(simLength(object)),
                            st_dist = abs(rnorm(simLength(object)))+0.05,
                            st_dist_starts = 1,
                            #grouping = rep(1,simLength(object)),
                            grouping = Data.Grouping,
                            step_distribution = grep(distribution(object),
                                                     c("gamma",
                                                       "log-normal")),
                            step_zero_inflation = 0,
                            angle_zero_inflation = 0
                            )

                Tpm.Working.Pars<- parameters(object)$Transition.Probabilities
                for( i in seq(nrow(Tpm.Working.Pars)) ) {
                    Tpm.Working.Pars[i,]<- Tpm.Working.Pars[i,] / Tpm.Working.Pars[i,i]
                }
                Tpm.Working.Pars<- -log(1/Tpm.Working.Pars - 1)


                Theta.Working.Pars<- parameters(object)$Deflection.Angle.Parameters
                Theta.Working.Pars[,2]<- -log(1/Theta.Working.Pars[,2] - 1)

                Step.Working.Pars<- parameters(object)$Step.Length.Parameters
                if( distribution(object) == "gamma" ) {
                    Step.Working.Pars[,1]<- log(Step.Working.Pars[,1])
                    Step.Working.Pars[,2]<- -log(1/Step.Working.Pars[,2] - 1)
                    Step.Working.Pars[,3]<- log(Step.Working.Pars[,3]-0.01)
                } else if( distribution(object) == "log-normal" ) {
                    Step.Working.Pars[,2]<- -log(1/Step.Working.Pars[,2] - 1)
                    Step.Working.Pars[,3]<- log(Step.Working.Pars[,3]-0.01)
                } else {}

                model<- setModel4M(N.States = nStates(object))
                Logit.Step.Zero.Probs<- startingValues(model)$Logit.Step.Zero.Probs
                Logit.Step.Zero.Probs[is.na(Logit.Step.Zero.Probs)]<- runif(sum(is.na(Logit.Step.Zero.Probs)),
                                                                            min = -3,
                                                                            max = 3)

                Logit.Angle.Zero.Probs<- startingValues(model)$Logit.Angle.Zero.Probs
                Logit.Angle.Zero.Probs[is.na(Logit.Angle.Zero.Probs)]<- runif(sum(is.na(Logit.Angle.Zero.Probs)),
                                                                            min = -3,
                                                                            max = 3)

                Angle.Zero.Working.Pars<- startingValues(model)$Angle.Zero.Working.Pars
                Angle.Zero.Working.Pars[is.na(Angle.Zero.Working.Pars)]<- runif(sum(is.na(Angle.Zero.Working.Pars)),
                                                                            min = -3,
                                                                            max = 3)

                if( useHMM(object) == T ){
                        if( distribution(object) == "gamma" ) {
                                Step.Working.Pars[,2]<- -Inf   ### autocorrelation to 0
                        } else if( distribution(object) == "log-normal" ) {
                                Step.Working.Pars[,2]<- 0
                        }
                } else {}

                Logit.Step.Zero.Probs[]<- -Inf

                Logit.Angle.Zero.Probs[]<- -Inf

                Angle.Zero.Working.Pars[,1]<- 0
                Angle.Zero.Working.Pars[,2]<- -Inf

                parameters<- list(tpm_working_pars_matrix = Tpm.Working.Pars,
                                  theta_working_pars = Theta.Working.Pars,
                                  dist_working_pars = Step.Working.Pars[,-2],
                                  acf_working_pars = Step.Working.Pars[,2],
                                  logit_step_zero_probs = Logit.Step.Zero.Probs,
                                  logit_angle_zero_probs = Logit.Angle.Zero.Probs,
                                  angle_zero_working_pars = Angle.Zero.Working.Pars,
                                  dummy = 0)

                map<- parameterMapping(model)
                if( useHMM(object) == T ){
                        if( distribution(object) == "gamma" ) {
                                Step.Working.Pars[,2]<- -Inf   ### autocorrelation to 0
                        } else if( distribution(object) == "log-normal" ) {
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

                map<- list(tpm_working_pars_matrix = map$Tpm.Working.Pars,
                           theta_working_pars = map$Theta.Working.Pars,
                           dist_working_pars = map$Step.Working.Pars[,-2],
                           acf_working_pars = map$Step.Working.Pars[,2],
                           logit_step_zero_probs = map$Logit.Step.Zero.Probs,
                           logit_angle_zero_probs = map$Logit.Angle.Zero.Probs,
                           angle_zero_working_pars = map$Angle.Zero.Working.Pars)

                map<- lapply(map,as.factor)

                ADFun<- TMB::MakeADFun(data = data,
                                       parameters = parameters,
                                       DLL = DLL,
                                       map = c(map,list(dummy = as.factor(NA))),
                                       silent = T)
                Report<- ADFun$report()
                TmbParameters<- list(
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

                Fab.Model4M<- model4M(new("Model4M"),setModel4M(N.States = nStates(object)))
                movementData(Fab.Model4M)$Movement.Data<- data.frame(Deflection.Angle = data$theta,
                                                                     Step.Length = data$st_dist,
                                                                     Group = 1)
                movementData(Fab.Model4M)$Step.Length.Starting.Values<- 1
                movementData(Fab.Model4M)$Step.Length.Mean<- 1
                useHMM(Fab.Model4M)<- useHMM(object)
                nStates(Fab.Model4M)<- nStates(object)
                interpolationParameters(Fab.Model4M)<- c(Time.Step = 1,Group.Cutoff = 2)
                parameterEstimates(Fab.Model4M)<- TmbParameters
                tmbEnvironment(Fab.Model4M)<- ADFun$env

                assign("test",Fab.Model4M,globalenv())

                Sim.4M<- simulate(Fab.Model4M,
                                  nsim = nsim,
                                  seed = seed,
                                  refit = refit,
                                  keep.locations = F,
                                  keep.data = keep.data,
                                  keep.residuals = keep.residuals,
                                  ...)

                return(Sim.4M)
            }
)



#' Compute simulation parameter errors.
#'
#' @param x An object of class \code{Simulate4M}.
#' @param center.fun What function should be used to comptue the center? Defaults to mean.
#'
#' @return Functions \code{rmse} and \code{mae} return a list of parameter errors. Function \code{stateError} returns a numeric vector containing the percentage of state estimates that are not correct.
#'
#' @examples
#' sealData<- data4M(greyseal)
#' sealData<- interpolate(sealData,Time.Step = 1)
#' seal4M2<- fit(sealData)
#' sim4M2<- simulate(seal4M2,nsim = 5)
#' rmse(sim4M2)
#' stateError(sim4M2)
#'
#' @name simulationError
NULL

#' @noRd
setGeneric(name = "simEstErrors",
           def = function(x,type,center.fun) standardGeneric("simEstErrors")
)
#' @noRd
setMethod(f = "simEstErrors",
          signature = "Simulate4M",
          definition = function(x,type,center.fun) {
                if( !(type %in% c("rmse","mae","bias")) ) {
                        stop("Error type must be either rmse, mae, or bias.")
                } else {}

                True.Parameters<- parameterEstimates(x)
                Refit.Parameters<- refitParameters(x)

                sweeper<- function(refit,true,diffs.margin,mean.margin,type = type) {
                        diffs<- sweep(refit,MARGIN = diffs.margin, STATS = true)
                        if( type == "rmse" ) {
                                diffs<- diffs^2
                        } else if( type == "mae" ) {
                                diffs<- abs(diffs)
                        } else if( type == "bias" ) {
                                diffs<- diffs
                        } else {}

                        error<- apply(diffs, MARGIN = mean.margin,center.fun)
                        if( type == "rmse" ) {
                                error<- sqrt(error)
                        } else {}

                        return(error)
                }

                Angle.Error<- sweeper(refit = Refit.Parameters$Deflection.Angle.Parameters,
                                      true = data.matrix(True.Parameters$Deflection.Angle.Parameters),
                                      diffs.margin = c(1,2),
                                      mean.margin = c(1,2),
                                      type = type)
                Step.Error<- sweeper(refit = Refit.Parameters$Step.Length.Parameters,
                                     true = data.matrix(True.Parameters$Step.Length.Parameters),
                                     diffs.margin = c(1,2),
                                     mean.margin = c(1,2),
                                     type = type)
                TPM.Error<- sweeper(refit = Refit.Parameters$Transition.Probability.Matrix,
                                    true = data.matrix(True.Parameters$Transition.Probability.Matrix),
                                    diffs.margin = c(1,2),
                                    mean.margin = c(1,2),
                                    type = type)
                Stat.Dist.Error<- sweeper(refit = Refit.Parameters$Stationary.Distribution,
                                          true = True.Parameters$Stationary.Distributio,
                                          diffs.margin = 2,
                                          mean.margin = 2,
                                          type = type)

                return(list(Deflection.Angle.Parameters = Angle.Error,
                            Step.Length.Parameters = Step.Error,
                            Transition.Probability.Matrix = TPM.Error,
                            Stationary.Distribution = Stat.Dist.Error))
          }
)

#' @export
#' @rdname simulationError
setGeneric(name = "rmse",
           def = function(x,...) standardGeneric("rmse")
)
#' @export
setMethod(f = "rmse",
          signature = "Simulate4M",
          def = function(x,center.fun=mean) return(simEstErrors(x,type = "rmse",center.fun))
)

#' @export
#' @rdname simulationError
setGeneric(name = "mae",
           def = function(x,...) standardGeneric("mae")
)
#' @export
setMethod(f = "mae",
          signature = "Simulate4M",
          def = function(x,center.fun=mean) return(simEstErrors(x,type = "mae",center.fun))
)

#' @export
#' @rdname simulationError
setGeneric(name = "bias",
           def = function(x,...) standardGeneric("bias")
)
#' @export
setMethod(f = "bias",
          signature = "Simulate4M",
          def = function(x,center.fun=mean) return(simEstErrors(x,type = "bias",center.fun))
)

#' @export
#' @rdname simulationError
setGeneric(name = "stateError",
           def = function(x) standardGeneric("stateError")
)
#' @export
setMethod(f = "stateError",
          signature = "Simulate4M",
          def = function(x) {
                return(rowMeans(simulatedViterbiPath(x) != refitViterbiPath(x)))
          }
)


#' Check simulations for undetected non-convergence.
#'
#' @param x An object of class Simulate4M
#'
#' @return Returns an index vector which flags simulations which violate conditions.
#' @rdname simCheck
setGeneric(name = "simCheck",
           def = function(x,...) standardGeneric("simCheck")
)
#' @param Stationary.Distribution.Check Any simulation which has a stationary distribution having any element less than this value will be flagged.
#' @param Angle.Concentration.Check Any simulation which has a deflection angle concentration less than this value will be flagged.
#' @param Tpm.Tolerance Any simulation which has a transition probability matrix which contains a row with at least two entries which
#'     are equal up to this many decimal places will be flagged.
#'
#' @export
#' @rdname simCheck
setMethod(f = "simCheck",
          signature = "Simulate4M",
          def = function(x,
                         Stationary.Distribution.Check = 0.05,
                         Angle.Concentration.Check = 0.05,
                         Tpm.Tolerance = 4) {
                Convergence.Index<- which(refitConvergence(x) != "relative convergence (4)")

                Stationary.Distribution.Index<- sapply(seq_len(nStates(x)), function(i) {
                        return(which(refitParameters(x)$Stationary.Distribution[,i] < Stationary.Distribution.Check))
                })
                if( is.list(Stationary.Distribution.Index) ) {
                Stationary.Distribution.Index<- do.call(c,Stationary.Distribution.Index)
                }

                Concentration.Index<- which(refitParameters(x)$Deflection.Angle.Parameters[1,2,] < Angle.Concentration.Check)

                Tpm<- refitParameters(x)$Transition.Probability.Matrix
                Tpm<- round(Tpm,Tpm.Tolerance)

                Tpm.Index<- apply(Tpm,MARGIN = 3,function(i) {
                  Row.Test<- apply(i,MARGIN = 1,function(j) {
                    if( anyDuplicated(j) ) {
                      return(TRUE)
                    } else {
                      return(FALSE)
                    }
                  })
                  return(any(Row.Test))
                })
                Tpm.Index<- which(Tpm.Index)

                Index<- c(Convergence.Index,Stationary.Distribution.Index,Concentration.Index)

                if( length(Index) > 0 ) {
                        Index<- sort(unique(Index))
                }

                return(Index)
          }
)
