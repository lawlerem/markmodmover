#' @include classes.R getset.R
NULL

#' Print method for objects of class \code{Data4M}.
#'
#' @export
#' @noRd
setMethod(
        f = "show",
        signature = "Data4M",
        definition = function(object) {
                Observed.Locations<- observedLocations(object)
                N.Obs.Locs<- nrow(Observed.Locations)
                Time.Frame<- range(Observed.Locations$Date)
                Lon.Extent<- range(Observed.Locations$Lon)
                Lat.Extent<- range(Observed.Locations$Lat)

                N.Int.Locs<- nrow(interpolatedLocations(object))
                Time.Step<- interpolationParameters(object)["Time.Step"]
                Group.Cutoff<- interpolationParameters(object)["Group.Cutoff"]
                N.Groups<- length(levels(groups(object)))

                cat("",
                    "Observed Locations",
                    "------------------",
                    paste("Length:",N.Obs.Locs,"locations."),
                    paste("Date Range:",Time.Frame[[1]],"to",Time.Frame[[2]]),
                    paste("Spatial Range:",
                          round(Lon.Extent[[1]],4),
                          " to ",
                          round(Lon.Extent[[2]],4),
                          " longitude \n",
                          "              ",
                          round(Lat.Extent[[1]],4),
                          " to ",
                          round(Lat.Extent[[2]],4),
                          " latitude."),
                    "",
                    "Time Differences",
                    sep = "\n")

                print(summary(difftime(object)))

                cat("",
                    "",
                    "Interpolated Locations",
                    "----------------------",
                    paste("Length:",N.Int.Locs,"locations."),
                    paste("Contiguous Groups:",N.Groups,"groups."),
                    paste("Time step:",Time.Step/(60*60),"hrs             Group cutoff:",
                          Group.Cutoff/(60*60),"hrs"),
                    "",
                    sep = "\n"
                    )

                return(invisible())
      }
)


#' Print method for objects of class \code{SetModel4M}.
#'
#' @export
#' @noRd
setMethod(
        f = "show",
        signature = "SetModel4M",
        definition = function(object) {
                cat("",
                    "Model Definition",
                    "----------------",
                    paste("Number of States:",nStates(object)),
                    paste("HMM:",useHMM(object)),
                    paste("Step Length Distribution:",distribution(object)),
                    paste("Step Zero Inflation:",zeroInflation(object)[["Step.Length"]]),
                    paste("Angle Zero Inflation:",zeroInflation(object)[["Deflection.Angle"]]),
                    "",
                    sep = "\n")

               return(invisible())
        }
)


#' Print method for objects of class \code{Model4M}.
#'
#' @export
#' @noRd
setMethod(f = "show",
          signature = "Model4M",
          definition = function(object) {
                show(data4M(object))
                show(SetModel4M(object))

                Parameters<- parameterEstimates(object)

                cat("",
                    "",
                    paste("Convergence:",convergence(object)),
                    "",
                    "",
                    "Deflection Angle Parameter Estimates",
                    "------------------------------------",
                    sep = "\n")

                print(round(Parameters$Deflection.Angle.Parameters,4))

                if( zeroInflation(object)[["Deflection.Angle"]] == T ) {
                        cat("",
                            paste("Zero-Inflation Probability:",
                                  paste(round(Parameters$Zero.Inflation$Angle.Zero.Probs,4),
                                        collapse = " ")
                                  ),
                            "",
                            sep = "\n")

                        print(round(Parameters$Zero.Inflation$Angle.Zero.Pars,4))
                }

                cat("",
                    "",
                    "Step Length Parameter Estimates",
                    "-------------------------------",
                    "",
                    paste("Mean step length:",
                          round(movementData(object)$Step.Length.Mean,4),
                          "km"),
                    "",
                    sep = "\n")

                print(round(Parameters$Step.Length.Parameters,4))

                if( zeroInflation(object)[["Step.Length"]] == T ) {
                        cat("",
                            paste("Zero-Inflation Probability:",
                                  paste(round(Parameters$Zero.Inflation$Step.Zero.Probs,4),
                                        collapse = " ")
                                 ),
                           sep = "\n")
                }

                cat("",
                    "",
                    "Transition Probabilities",
                    "------------------------",
                    "",
                    "Transition Probability Matrix",
                    "",
                    sep = "\n")

                print(round(Parameters$Transition.Probability.Matrix,4))

                cat("",
                    paste("Stationary Distribution / Activity Budget:",
                          paste(round(Parameters$Stationary.Distribution,4),
                                collapse = " ")
                         ),
                    sep = "\n")

                return(invisible())

          }
)




#' Plot objects of class \code{Data4M} or \code{Model4M}.
#'
#' @param x An object of class \code{Data4M} or \code{Model4M}.
#' @param y A character vector choosing which plots to draw. Can be any combination of \code{all}, \code{satellite}, \code{locations}, \code{data}, and \code{residuals}. Defaults to \code{all}.
#'
#' @name plot4M
NULL

#' @param ask Wait for input to switch to next plot?
#' @param Viterbi.Path State path used to color locations. Not typically entered by the user.
#' @param N.States The number of states used in the model. Not typically entered by the user.
#'
#' @examples
#' sealdata<- data4M(greyseal)
#' sealdata<- interpolate(sealdata,Time.Step = 1)
#' plot(sealdata,ask = F)
#' plot(sealdata,"data",log = T,lag = 2,N.Grid.Lines = 150,N.Greys = 30)
#' plot(sealdata,"locations",pch = 11,lty = 4,col = "red")
#'
#' seal4M3<- fit(sealdata,N.States = 3)
#' plot(seal4M3,ask = F)
#' plot(seal4M3,"residuals")
#'
#' library(RColorBrewer)
#' plot(seal4M3,"locations",palette = "Dark2")
#'
#' @export
#' @rdname plot4M
setMethod(
        f = "plot",
        signature = c(x = "Data4M",
                      y = "character"),
        definition = function(x,
                              y,
                              ask = interactive(),
                              Viterbi.Path = NA,
                              N.States = 0,
                              h = 0.25,
                              ...) {

                y<- sapply(y,function(elem) {
                        grep(elem,
                             c("all","satellite","locations","data"),
                             ignore.case = T,
                             value = T)
                })

                if( "all" %in% y ) {
                        y<- c("satellite","locations","data")
                } else {}

                Optional.Args<- list(...)

                if( !("lag" %in% names(Optional.Args)) ) {
                        Optional.Args$lag<- 1
                } else {}
                Lag<- Optional.Args$lag; Optional.Args$lag<- NULL

                if( !("log" %in% names(Optional.Args)) ) {
                        Optional.Args$log<- F
                } else {}
                Log<- Optional.Args$log; Optional.Args$log<- NULL

                if( !("N.Grid.Lines" %in% names(Optional.Args)) ) {
                        Optional.Args$N.Grid.Lines<- 100
                } else {}
                N.Grid.Lines<- Optional.Args$N.Grid.Lines; Optional.Args$N.Grid.Lines<- NULL

                if( !("N.Greys" %in% names(Optional.Args)) ) {
                        Optional.Args$N.Greys<- 25
                } else {}
                N.Greys<- Optional.Args$N.Greys; Optional.Args$N.Greys<- NULL





                if( "satellite" %in% y ) {
                        par(mfrow=c(1,1))

                        Plot.Args<- list(main = "Satellite Locations",
                                         xlab = "Longitude",
                                         ylab = "Latitue",
                                         type = "l",
                                         pch = 20)
                        Plot.Args[names(Optional.Args)]<- Optional.Args

                        do.call(plot,
                                c(list(x = observedLocations(x)$Lon,
                                       y = observedLocations(x)$Lat),
                                  Plot.Args))

                        y[which(y == "satellite")]<- NA
                        if( any(!is.na(y)) & ask == T ) {
                                readline("Press <Enter> for next plot. ")
                        } else {}
                }


                if( "locations" %in% y ) {
                        par(mfrow=c(1,1))

                        Plot.Args<- list(main = "Interpolated Locations",
                                         xlab = "Longitude",
                                         ylab = "Latitude",
                                         col = "black")
                        Plot.Args[names(Optional.Args)]<- Optional.Args
                        Plot.Args$type<- "n"

                        Locations<- interpolatedLocations(x)

                        do.call(plot,
                                c(list(x = range(Locations$Lon),
                                       y = range(Locations$Lat) + c(0,0.2*diff(range(Locations$Lat)))),
                                       Plot.Args))

                        if( !("pch" %in% names(Plot.Args)) ) {
                                Plot.Args$pch<- 20
                        }
                        Plot.Args$pch<- rep_len(Plot.Args$pch,
                                                length.out = N.States+1)

                        Plot.Args$col<- rep_len(Plot.Args$col,
                                                length.out = N.States+1)

                        if( !("cex" %in% names(Plot.Args)) ) {
                                Plot.Args$cex<- 1.0
                        }

                        invisible(
                                by(Locations,
                                   Locations$Group,
                                   function(Loc) {
                                        lines(Loc$Lon,
                                              Loc$Lat)

                                        if( !((length(Viterbi.Path) == 1) || anyNA(Viterbi.Path)) ) {
                                                Color.Path<- c(0,
                                                               Viterbi.Path[Viterbi.Path[,2] == Loc$Group[[1]],1],
                                                               0)+1
                                        } else {
                                                Color.Path<- 1
                                        }

                                        points(Loc$Lon,
                                               Loc$Lat,
                                               pch = Plot.Args$pch[Color.Path],
                                               cex = Plot.Args$cex,
                                               col = Plot.Args$col[Color.Path])
                                               #col = "black")
                                   })
                        )

                        if( N.States != 0 ) {
                                legend(x = "top",
                                       legend = paste("State",seq(N.States)),
                                       col = Plot.Args$col[-1],
                                       pch = Plot.Args$pch[-1],
                                       horiz = T)
                        }

                        y[which(y == "locations")]<- NA

                        if( any(!is.na(y)) & ask == T ) {
                                readline("Press <Enter> for next plot. ")
                        }
                }

                if( "data" %in% y ) {

                        Plot.Args<- list(main = "",
                                         xlab = "",
                                         ylab = "")
                        Plot.Args[names(Optional.Args)]<- Optional.Args

                        Movement.Data<- movementData(x)$Movement.Data

                        Step.Length<- Movement.Data$Step.Length

                        Deflection.Angle<- Movement.Data$Deflection.Angle
                        Cartesian.Data.X<- Step.Length*cos(pi/180 * Deflection.Angle)
                        Cartesian.Data.Y<- Step.Length*sin(pi/180 * Deflection.Angle)


                        par(mfrow = c(2,1))

                        Plot.Args$col<- grey(seq(1,0,length.out = N.Greys))
                        Plot.Args$main<- "Deflection Angle and Step Length"

                        Cartesian.Density.Values<- MASS::kde2d(Cartesian.Data.X,
                                                         Cartesian.Data.Y,
                                                         n = N.Grid.Lines)
                        do.call(image,
                                c(list(x = Cartesian.Density.Values),
                                  Plot.Args)
                               )



                        if( Log == T ) {
                                Step.Length<- log(Step.Length)
                        }

                        Step.Min<- min(Step.Length)
                        Step.Max<- max(Step.Length)

                        Plot.Args$main<- "Step Length Lag Plot"
                        if( Log == T ) {
                                Plot.Args$main<- paste("(Log)",Plot.Args$main)
                        }
                        Plot.Args$xlab<- "Step length at time t-1"
                        Plot.Args$ylab<- "Step length at time t"

                        Step.Density.Values<- MASS::kde2d(head(Step.Length,-Lag),
                                                    tail(Step.Length,-Lag),
                                                    n = N.Grid.Lines,
                                                    lims = c(Step.Min,
                                                             Step.Max,
                                                             Step.Min,
                                                             Step.Max)
                                                    )

                        do.call(image,
                                c(list(x = Step.Density.Values),
                                  Plot.Args)
                                )

                        if( Log == T ) {
                                Step.Length<- exp(Step.Length)
                        }
                }
        }
)
#' @export
#' @noRd
setMethod(
        f = "plot",
        signature = c(x = "Data4M",
                      y = "missing"),
        definition = function(x,
                              y,
                              ask = interactive(),
                              ...) {
                y<- "all"
                plot(x,y,ask = ask,...)
        }
)


#' @param h Vector of bandwidth for residual density plots. See \code{\link[MASS]{kde2d}}.
#' @param palette If using package \code{RColorBrewer} and not supplying \code{col}, the palette to use for colors.
#' @param ... Optional arguments used to control internal plotting functions, such as \code{pch} or \code{main}. Other useful arguments include...
#'   \describe{
#'     \item{lag}{For data lag-plot, which time lag should be used? Defaults to 1.}
#'     \item{log}{For data lag-plot, should the log of step lengths be used? Defaults to FALSE.}
#'     \item{N.Grid.Lines}{How many grid lines should be used for each axis in the mesh for the data plots? Defaults to 100.}
#'     \item{N.Greys}{What greyscale depth to use for density plots of the data. Defaults to 25.}
#'   }
#'
#' @export
#' @rdname plot4M
setMethod(f = "plot",
          signature = c(x = "Model4M",
                        y = "character"),
          definition = function(x,
                                y,
                                h = 0.25,
                                ask = interactive(),
                                palette = "Set1",
                                ...) {

                y<- sapply(y,function(elem) {
                        grep(elem,
                        c("all","satellite","locations","data","residuals"),
                        ignore.case = T,
                        value = T)
                })

                if( "all" %in% y ) {
                        y<- c("satellite","locations","data","residuals")
                } else{}

                Plot.Args<- list(...)

                if( "satellite" %in% y ) {

                        plot(x = data4M(x),
                             y = "satellite")

                        y[which(y == "satellite")]<- NA
                        if( any(!is.na(y)) & ask == T ) {
                                readline("Press <Enter> for next plot. ")
                        } else {}
                } else {}



                if( "locations" %in% y ) {
                        if( "col" %in% names(Plot.Args) ) {
                                colors<- rep_len(Plot.Args$col,
                                                 length.out = nStates(x)+1)
                        } else if( requireNamespace("RColorBrewer",quietly = T) ) {
                                colors<- c("black",RColorBrewer::brewer.pal(max(3,nStates(x)),palette))
#                        } else if( "RColorBrewer" %in% .packages(all.available = T) ) {
#                               answer<- ""
#                               while( !(answer %in% c("y","n")) ) {
#                                       answer<- readline("Load package RColorBrewer? [y/n] ")
#                              }
#
#                                if( answer == "y" ) {
#                                        require(RColorBrewer)
#                                        colors<- c("black",brewer.pal(max(3,nStates(x)),palette))
#                                } else {
#                                        colors<- c("black",palette()[which(palette() != "black")])
#                                }
                        } else {
                                warning("Consider installing RColorBrewer for color-coding the state path.")
                                colors<- "black"
                        }

                        do.call(plot,
                                c(list(x = data4M(x),
                                     y = "locations",
                                     ask = F,
                                     Viterbi.Path = data.frame(viterbiPath(x),
                                                               movementData(x)$Movement.Data$Group),
                                     N.States = nStates(x),
                                     col = colors),
                                  Plot.Args[-which(names(Plot.Args) != "col")]
                                  )
                               )

                        y[which(y == "locations")]<- NA

                        if( any(!is.na(y)) & ask == T ) {
                                readline("Press <Enter> for next plot. ")
                        } else {}
                } else {}


                if( "data" %in% y ) {
                        plot(data4M(x),
                             y = "data",
                             ...)

                        y[which(y == "data")]<- NA
                        if( any(!is.na(y)) & ask == T ) {
                                readline("Press <Enter> for next plot. ")
                        } else {}
                } else {}


                if( "residuals" %in% y ) {
                        par(mfrow = c(2,2))

                        ###
                        ### acf for dist      --   lagplot for dist
                        ### qqplot for theta  --    qqplot for dist
                        ###
                        Args<- list(main = "Step Length ACF",
                                    ylab = "Correlation",
                                    xlab = "Lag",
                                    pch = 20)
                        Args[names(Plot.Args)]<- Plot.Args
                        do.call(acf,
                                c(list(x = residuals(x)$Step.Length),
                                  Args)
                               )

                        if( !("lag" %in% names(Args) ) ) {
                                Args$lag<- 1
                        } else {}
                        Lag<- Args$lag; Args$lag<- NULL

                        if( !("log" %in% names(Args)) ) {
                                Args$log<- NULL
                        } else {}

                        if( !("N.Grid.Lines" %in% names(Args)) ) {
                                Args$N.Grid.Lines<- 100
                        } else {}
                        N.Grid.Lines<- Args$N.Grid.Lines; Args$N.Grid.Lines<- NULL

                        if( !("N.Greys" %in% names(Args)) ) {
                                Args$N.Greys<- 10
                        } else {}
                        N.Greys<- Args$N.Greys; Args$N.Greys<- NULL




                        Args<- list(main = "Step Length Residuals",
                                    xlab = "Residual at time t-1",
                                    ylab = "Residual at time t")
                        Args[names(Plot.Args)]<- Plot.Args

                        Residuals<- residuals(x)$Step.Length

                        Args$col<- grey(seq(1,0,length.out = N.Greys))

                        Residual.Density.Values<- MASS::kde2d(head(Residuals,-Lag),
                                                        tail(Residuals,-Lag),
                                                        n = N.Grid.Lines,
                                                        h = h,
                                                        lims = c(-1,1,-1,1))
                        do.call(image,
                                c(list(x = Residual.Density.Values),
                                  Args)
                                )


                        Args$main<- "Deflection Angle Q-Q Plot"
                        Args$xlab<- "Uniform(-1,1) Distribution"
                        Args$ylab<- "Residuals"
                        Args$pch<- 20
                        do.call(qqplot,
                                c(list(y = residuals(x)$Deflection.Angle,
                                       x = seq(-1,1,2/nrow(residuals(x)))),
                                  Args)
                               )
                        abline(a = 0, b = 1, col = "red")

                        Args$main<- "Step Length Q-Q Plot"
                        do.call(qqplot,
                                c(list(y = residuals(x)$Step.Length,
                                       x = seq(-1,1,2/nrow(residuals(x)))),
                                  Args)
                               )
                        abline(a = 0, b = 1, col = "red")

                        par(mfrow = c(1,1))

                }

                return(invisible())
        }
)
