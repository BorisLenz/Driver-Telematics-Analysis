TripMatchingFun <- function(driver, returnMatchedTrips = FALSE){
  
  #This function aligns and centers trajectories to later perform trajectory
  #matching based on qDistances or hausdorff distances
  
  require("data.table")
  require("parallel")
  require("sp")
  require("maptools")
  require("rgeos")
  
  #windows libraries
  library(foreach)  
  library(doParallel)  
  
  realignFun <- function(tripNum, driverN){
    
    trip <- fread(file.path(driversDirectory, driverN, paste0(tripNum, ".csv")))
    coordinatesDf <- as.data.frame(trip)
    lastPoint <- coordinatesDf[nrow(coordinatesDf), ]     
    coordinates(coordinatesDf) <- ~x+y
    originCoordinates <- elide(coordinatesDf, rotate = atan2(lastPoint$y, lastPoint$x) * (180/pi), center = c(0, 0))
    
    return(originCoordinates)
  }
  
  realignedTrajectories <- mclapply(seq(1, 200), realignFun, mc.cores = numCores, driverN = driver)
    
  TrajectoryParams <- function(xyTrip){
    #This function calculates: max, min, which.max, which.min of x and y parameters
    maxX <- max(xyTrip$x)
    minX <- min(xyTrip$x)
    maxXIdx <- which.max(xyTrip$x)
    minXIdx <- which.min(xyTrip$x)
    maxY <- max(xyTrip$y)
    minY <- min(xyTrip$y)
    maxYIdx <- which.max(xyTrip$y)
    minYIdx <- which.min(xyTrip$y)
  }
  
  #Similarity Measuring
  are.similar <- function(numeric1, numeric2, threshold){        
    #this function calculates if two parameters are close to one another         
    if(numeric1 > numeric2 - (numeric2 * threshold) &
         numeric1 < numeric2 + (numeric2 * threshold)){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }
  
  #Matching Function
  MatchingFun <- function(factors1, factors2){
    factorsDf <- as.data.frame(cbind(factors1, factors2, c(100, 100, 10, 10, 100, 100, 10, 10)))
    similarValues <- apply(factorsDf, 1, are.similar)
    if(sum(similarValues) == 8){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }    
  
  calculateFactors <- function(tripNu, tripToBeMatched){
    #calculate mins and maxs of centered and rotated trajetories
    #Mins and Maxs Original Trajectory
    trip2MatchFactors <- TrajectoryParams(tripToBeMatched)
    
    #Other Original Trajectories
    nthTrip <- realignedTrajectories[[tripNu]]  
    #Mins and Maxs 
    nthTripFactors <- TrajectoryParams(nthTrip)
    
    #Reflected trajetories
    reflectedTrip1 <- elide(nthTrip, reflect = c(FALSE, TRUE), center = c(0, 0))    
    #Mins and Maxs
    reflectedTrip1Factors <- TrajectoryParams(reflectedTrip1)
    
    reflectedTrip2 <- elide(nthTrip, reflect = c(TRUE, FALSE), center = c(0, 0))    
    #Mins and Maxs
    reflectedTrip2Factors <- TrajectoryParams(reflectedTrip2)
    
    reflectedTrip3 <- elide(nthTrip, reflect = c(TRUE, TRUE), center = c(0, 0))    
    #Mins and Maxs
    reflectedTrip3Factors <- TrajectoryParams(reflectedTrip3)         
    
    #Determine whether there is a matching trajectory or not (return a true or false statement)    
    if (MatchingFun(trip2MatchFactors, reflectedTrip1Factors) == TRUE |
          MatchingFun(trip2MatchFactors, reflectedTrip1Factors) == TRUE |
          MatchingFun(trip2MatchFactors, reflectedTrip2Factors) == TRUE |
          MatchingFun(trip2MatchFactors, reflectedTrip3Factors) == TRUE){
      return(TRUE)
    }else{
      return(FALSE) 
    }
  }
  
  #Return matched trips
  tripsMatched <- lapply(1:200, function(tripNumber){
    matches <- mclapply(seq(1, 200), calculateFactors, mc.cores = numCores,
                        tripToBeMatched = realignedTrajectories[[tripNumber]])
    if (returnMatchedTrips == FALSE){
      if(sum(matches > 0)){
        return(TRUE) 
      }else{
        return(FALSE)
      }           
    }else{
      if(sum(matches > 0)){
        return(list(TRUE, which(matches)))   #here for debugging only 
      }else{
        return(FALSE)
      }      
    }    
  })  
  
  #Give an ID to the trips matched (They are the same as th submission file)
  names(tripsMatched) <- as.character(paste0(driver, seq(1, 200)))
  #Return Data
  return(tripsMatched)  
}