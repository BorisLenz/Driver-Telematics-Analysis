TripMatchingFun <- function(driver){
  
  require("data.table")
  require("sp")
  require("maptools")
  require("rgeos")
  
  realignFun <- function(trip, driverN){
    trip <- fread(file.path(driversDirectory, driverN, paste0(trip, ".csv")))
    coordinatesDf <- as.data.frame(trip)
    lastPoint <- coordinatesDf[nrow(coordinatesDf), ]     
    coordinates(coordinatesDf) <- ~x+y
    originCoordinates <- elide(coordinatesDf, rotate = atan2(lastPoint$y, lastPoint$x) * (180/pi))
    return(originCoordinates)
  }
   realignedTrajectories <- mclapply(seq(1, 200), realignFun, mc.cores = numCores, driverN = driver)
  
  
}