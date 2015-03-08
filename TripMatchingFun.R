TripMatchingFun <- function(driver){
  
  #This function aligns and centers trajectories to later perform trajectory
  #matching based on qDistances or hausdorff distances
  
  require("data.table")
  require("sp")
  require("maptools")
  require("rgeos")
  
  realignFun <- function(tripNum, driverN){
    
    trip <- fread(file.path(driversDirectory, driverN, paste0(tripNum, ".csv")))
    coordinatesDf <- as.data.frame(trip)
    lastPoint <- coordinatesDf[nrow(coordinatesDf), ]     
    coordinates(coordinatesDf) <- ~x+y
    originCoordinates <- elide(coordinatesDf, rotate = atan2(lastPoint$y, lastPoint$x) * (180/pi), center = c(0, 0))
    
    return(originCoordinates)
  }
   realignedTrajectories <- mclapply(seq(1, 200), realignFun, mc.cores = numCores, driverN = driver)
  
  calculateFactors <- function(trip, tripToBeMatched){
    #calculate mins and maxs of centered and rotated trajetories
    #Original Trajectory
    trip2MatchMaxX <- max(tripToBeMatched$x)
    trip2MatchMinX <- min(tripToBeMatched$x)
    trip2MatchMaxY <- max(tripToBeMatched$y)
    trip2MatchMinY <- min(tripToBeMatched$y)
    #Other Original Trajectories
    nthTrip <- realignedTrajectories[[trip]]  
    #Mins and Maxs
    maxTrip <- max(nthTrip$x)
    
    
    
    
#     distanceRaw <- gDistance(tripToBeMatched, realignedTrajectories[[trip]], hausdorff = TRUE)
#     #reflected trajectories
#     reflectedDistance1 <- gDistance(tripToBeMatched, elide(realignedTrajectories[[trip]], reflect = c(FALSE, TRUE),
#                                                            center = c(0, 0)), hausdorff = TRUE)
#     reflectedDistance2 <- gDistance(tripToBeMatched, elide(realignedTrajectories[[trip]], reflect = c(TRUE, FALSE), 
#                                                            center = c(0, 0)), hausdorff = TRUE)
#     reflectedDistance3 <- gDistance(tripToBeMatched, elide(realignedTrajectories[[trip]], reflect = c(TRUE, TRUE),
#                                                            center = c(0, 0)), hausdorff = TRUE)
#     
#     return(c(distanceRaw, reflectedDistance1, reflectedDistance2, reflectedDistance3))
  }  
  
  distances <- mclapply(seq(1, 200), calculateFactors, mc.cores = numCores,
                                    tripToBeMatched = realignedTrajectories[[141]])
  
  
}