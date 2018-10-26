.r.Earth<- 6371

.distance.formula<- function(x.last,
                             x.curr,
                             lat = 1,
                             lon = 2,
                             convert.to.radians = TRUE) {
   if(lat == lon) warning("Latitude index same as Longitude index")

   if(convert.to.radians == TRUE) {
      x.curr[[lat]]<- (pi/180)*x.curr[[lat]]
      x.curr[[lon]]<- (pi/180)*x.curr[[lon]]
      x.last[[lat]]<- (pi/180)*x.last[[lat]]
      x.last[[lon]]<- (pi/180)*x.last[[lon]]
   }

   inner.first<- sin(0.5*(x.curr[[lat]]-x.last[[lat]]))^2
   inner.second<- cos(x.curr[[lat]])*cos(x.last[[lat]])
   inner.third<- sin(0.5*(x.curr[[lon]]-x.last[[lon]]))^2
   dist<- 2*.r.Earth*asin(sqrt(inner.first+inner.second*inner.third))
   return(dist)
}

.heading.formula<- function(x.last,
                            x.curr,
                            lat = 1,
                            lon = 2,
                            convert.to.radians = TRUE) {
   if(lat == lon) warning("Latitude index same as Longitude index")

   if(convert.to.radians == TRUE) {
      x.curr[[lat]]<- (pi/180)*x.curr[[lat]]
      x.curr[[lon]]<- (pi/180)*x.curr[[lon]]
      x.last[[lat]]<- (pi/180)*x.last[[lat]]
      x.last[[lon]]<- (pi/180)*x.last[[lon]]
   }

   atan.arg1<- sin(x.curr[[lon]] - x.last[[lon]])*cos(x.curr[[lat]])
   atan.arg2<- cos(x.last[[lat]])*sin(x.curr[[lat]]) - sin(x.last[[lat]])*cos(x.curr[[lat]])*cos(x.curr[[lon]] - x.last[[lon]])

   heading<- 180/pi*atan2(atan.arg1, atan.arg2)
   return(heading)
}

.deflection.angle.formula<- function(x.last,
                                     x.curr,
                                     x.next,
                                     lat = 1,
                                     lon = 2,
                                     convert.to.radians = TRUE) {
   if(lat == lon) warning("Latitude index same as Longitude index")

   if(convert.to.radians == TRUE) {
      x.next[[lat]]<- (pi/180)*x.next[[lat]]
      x.next[[lon]]<- (pi/180)*x.next[[lon]]
      x.curr[[lat]]<- (pi/180)*x.curr[[lat]]
      x.curr[[lon]]<- (pi/180)*x.curr[[lon]]
      x.last[[lat]]<- (pi/180)*x.last[[lat]]
      x.last[[lon]]<- (pi/180)*x.last[[lon]]
   }

   heading.initial<- .heading.formula(x.last = x.curr,
                                      x.curr = x.next,
                                      lat = lat,
                                      lon = lon,
                                      convert.to.radians = FALSE)
   heading.final<- 180 + .heading.formula(x.last = x.curr,
                                          x.curr = x.last,
                                          lat = lat,
                                          lon = lon,
                                          convert.to.radians = FALSE)
   #  heading.final<- heading.final - 360*floor((heading.final+180)/360)

   deflection.angle<- heading.final - heading.initial
   deflection.angle<- deflection.angle - 360*floor((deflection.angle+180)/360)

   return(deflection.angle)
}

.dest.Point<- function(x.curr,
                       f.heading.curr,
                       distance.curr,
                       deflection.angle,
                       distance,
                       del = NULL,
                       lat = 1,
                       lon = 2,
                       convert.to.radians = TRUE){
   if(convert.to.radians) {
      f.heading.curr<- pi/180*f.heading.curr
      deflection.angle<- pi/180*deflection.angle
      x.curr[[lat]]<- pi/180*x.curr[[lat]]
      x.curr[[lon]]<- pi/180*x.curr[[lon]]
   }

   if( !is.null(del) ) {
      new.dist<- del*distance.curr
   } else {
      new.dist<- distance
   }


   new.heading<- f.heading.curr - deflection.angle

   x.next.lat<- asin( sin(x.curr[[lat]])*cos(new.dist/.r.Earth) + cos(x.curr[[lat]])*sin(new.dist/.r.Earth)*cos(new.heading) )
   x.next.lon<- x.curr[[lon]] + atan2( sin(new.heading)*sin(new.dist/.r.Earth)*cos(x.curr[[lat]]),
                                       cos(new.dist/.r.Earth) - sin(x.curr[[lat]])*sin(x.next.lat))
   x.next<- c(0,0)
   x.next[[lat]]<- x.next.lat
   x.next[[lon]]<- x.next.lon
   return(180/pi*x.next)
}
