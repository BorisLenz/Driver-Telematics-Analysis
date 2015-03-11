TripMatchingFun <- function(driver, returnMatchedTrips = FALSE){
  
  #This function aligns and centers trajectories to later perform trajectory
  #matching based on qDistances or hausdorff distances
  
  require("data.table")
  require("parallel")
  require("sp")
  require("maptools")
  require("rgeos")
  
  #Rotate trajectories so that the first and last points are aligned with the x axis
  realignFun <- function(tripNum, driverN){
    
    trip <- fread(file.path(driversDirectory, driverN, paste0(tripNum, ".csv")))
    coordinatesDf <- as.data.frame(trip)
    lastPoint <- coordinatesDf[nrow(coordinatesDf), ]     
    coordinatesSp <- SpatialPoints(coordinatesDf)
    originCoordinates <- elide(coordinatesSp, rotate = atan2(lastPoint$y, lastPoint$x) * (180/pi), center = c(0, 0))
    
    return(originCoordinates)
  }
  
  TrajectoryParams <- function(xyTrip){
    #This function calculates: max, min, which.max, which.min of x and y parameters
    maxX <- max(coordinates(xyTrip)[, 1])
    minX <- min(coordinates(xyTrip)[, 1])
    maxXIdx <- which.max(coordinates(xyTrip)[, 1])
    minXIdx <- which.min(coordinates(xyTrip)[, 1])
    maxY <- max(coordinates(xyTrip)[, 2])
    minY <- min(coordinates(xyTrip)[, 2])    
    maxYIdx <- which.max(coordinates(xyTrip)[, 2])
    minYIdx <- which.min(coordinates(xyTrip)[, 2])    
    length <- length(coordinates(xyTrip))
    return(c(maxX, minX, maxY, minY, maxXIdx, minXIdx, maxYIdx, minYIdx, length))
  }
    
  realignedTrajectories <- mclapply(seq(1, 200), realignFun, mc.cores = numCores, driverN = driver)
  realignedTrajectoriesFactors <- lapply(realignedTrajectories, TrajectoryParams)
  
  #Reflect realigned trajectories
  ReflectTrajectories <- function(tripNum){
    realignedDf <- as.data.frame(realignedTrajectories[[tripNum]])
    realignedDf$y <- realignedDf$y * -1
    realignedSp <- SpatialPoints(realignedDf)
    return(realignedSp)
  }
  
  reflectedTrajectories <- mclapply(seq(1, 200), ReflectTrajectories, mc.cores = numCores)
  reflectedTrajectoriesFactors <- lapply(reflectedTrajectories, TrajectoryParams)
  
  #Trajectory Differences Calculation
  tripsMatched <- lapply(seq(1, 200), function(idx){
    selfDifferences <- abs(realignedTrajectoriesFactors[[idx]] - as.data.frame(realignedTrajectoriesFactors))
    selfDifferencesIdx <- selfDifferences < 80
    matchedSelfIdx <- which(colSums(selfDifferencesIdx) >= 8)
    matchedSelfIdx <- matchedSelfIdx[matchedSelfIdx != idx]
    reflectedDifferences <- abs(realignedTrajectoriesFactors[[idx]] - as.data.frame(reflectedTrajectoriesFactors))    
    reflectedDifferencesIdx <- reflectedDifferences < 80
    matchedReflectedIdx <- which(colSums(reflectedDifferencesIdx) >= 8)
    return(c(matchedSelfIdx, matchedReflectedIdx))
    })      
  
  numberTripsMatched <- sapply(tripsMatched, length)
  repeatedTrips <- as.data.frame(numberTripsMatched > 0)
  rownames(repeatedTrips) <- as.character(paste0(driver, "_", seq(1, 200)))
  colnames(repeatedTrips) <- "RepeatedTrip"
  
  #Write as .csv
  write.csv(repeatedTrips, file = file.path(tripMatchingDir, paste0(driver, ".csv")), row.names = TRUE)
  #Report progress
  print(paste0(which(list.files(driversDirectory) == driver), "/", length(list.files(driversDirectory))))  
    
  #Give an ID to the trips matched (They are the same as th submission file)
  #Return Data
  return(TRUE)  
}