#AXA Driver Telematics Analysis
#Ver 0.8.8 #  debugged EDA #7

#Init-----------------------------------------------
rm(list=ls(all=TRUE))

#Libraries, directories, options and extra functions----------------------
require("data.table")
require("parallel")
require("h2o")
require("DMwR")
require("prospectr")
require("adehabitatLT")
require("ggplot2")
#require("plotly")

#Set Working Directory
workingDirectory <- "/home/wacax/Wacax/Kaggle/AXA-Driver-Telematics-Analysis/"
setwd(workingDirectory)
driversDirectory <- "/home/wacax/Wacax/Kaggle/AXA-Driver-Telematics-Analysis/Data/drivers"
otherDataDirectory <- "/home/wacax/Wacax/Kaggle/AXA-Driver-Telematics-Analysis/Data/"
outputDirectory <- "/home/wacax/Wacax/Kaggle/AXA-Driver-Telematics-Analysis/Data/Output"
logsDirectory <- "/home/wacax/Wacax/Kaggle/AXA-Driver-Telematics-Analysis/Data/Logs"
vw77Dir = "/home/wacax/vowpal_wabbit-7.7/vowpalwabbit/"
#h2o location
h2o.jarLoc <- "/home/wacax/R/x86_64-pc-linux-gnu-library/3.1/h2o/java/h2o.jar"

#List all possible drivers identities
drivers <- list.files(driversDirectory)
#Detect available cores
numCores <- detectCores()

#Extra Functions
source(paste0(workingDirectory, "pretrainh2oNN.R"))
source(paste0(workingDirectory, "TrajectoryMatrixDistances.R"))

#Data Mining (Functions)------------------------
#Outlier Removal and transformation to quantiles
quantSigma <- function(vector, sigma = 5, returnVector = TRUE, movavWindow = 7){
  #n-sigma removal
  #above n sigmas
  vector <- vector[!vector < mean(vector) - ifelse(sum(vector > 0) > 1, sd(vector), 0) * sigma]
  #below n sigmas
  vectorWithoutOutliers <- vector[!vector > mean(vector) + ifelse(sum(vector > 0) > 1, sd(vector), 0) * sigma]
  if (length(vectorWithoutOutliers) <= movavWindow){
    vectorWithoutOutliers <- rep(0, movavWindow + 1)
  }
  smoothVectorWithoutOutliers <- movav(vectorWithoutOutliers, movavWindow)
    
  if (returnVector == TRUE){
    return(list(quantile(smoothVectorWithoutOutliers, seq(0.1, 1, by=0.1)), smoothVectorWithoutOutliers))
  }else{
    return(quantile(smoothVectorWithoutOutliers, seq(0.1, 1, by=0.1)))    
  }
}
#Transform data to speed distributions
velocitiesKH <- function(coordinates){
  speed <- 3.6 * sqrt(diff(coordinates$x)^2 + diff(coordinates$y)^2)
  return(speed)
}
#Angles Mining with Angle Rotation
#angleVector <- function(coordinatesAng){
#  coorDT <- as.data.table(cbind(diff(coordinatesAng$x), diff(coordinatesAng$y)))
  #Returns the angles between vectors
  #coordinatesAng is a 2D-array of shape (N,M) representing N vectors in M-dimensional space.
  #It returns a vector which is a 1D-array of values of shape (N-1,), with each value between 0 and pi.
#  dir2 <- coorDT[2:nrow(coorDT)]
#  dir1 <- coorDT[1:(nrow(coorDT) - 1)]
#  rawAngles <- acos(rowSums(dir1 * dir2) / ((sqrt(rowSums(dir1^2) * rowSums(dir2^2))) + 1e-06)) * (180/pi) #1e-06 is included to avoid a division by zero
  #NAs removal
#  angles <- na.omit(rawAngles)
  #set orthogonal values to zero
#  angles[angles >= 89] <- 0  
  #Outliers removal
#  anglesWoOutliers <- quantSigma(angles, returnVector = TRUE)
#  return(list(anglesWoOutliers[[1]], rawAngles, anglesWoOutliers[[2]]))
#}
#Define function to be passed as parallel
transform2Percentiles <- function(file, driverID){  
  trip <- fread(file.path(driversDirectory, driverID, paste0(file, ".csv")))
  #Velocities
  rawVelocities <- velocitiesKH(trip)
  #n-sigma removal
  velocityData <- quantSigma(rawVelocities)
  speedDist <- velocityData[[1]]
  speedSd <- sd(velocityData[[2]])
  #Velocitiy without stops
  speedDistWOStopsDist <- quantSigma(velocityData[[2]][velocityData[[2]] > 0.5], returnVector = FALSE)
  #Velocity without stops standard deviation
  speedDistWOStopsSd <- ifelse(sum(velocityData[[2]] > 0.5) > 1, sd(velocityData[[2]][velocityData[[2]] > 0.5]), 0)
  #Accelerations
  accelerationData <- quantSigma(diff(velocityData[[2]]), returnVector = TRUE)
  accelerationsDist <- accelerationData[[1]]
  #Accelerations Standard Deviations
  positiveAccelerationSd <- ifelse(sum(accelerationData[[2]] > 0) > 1, sd(accelerationData[[2]][accelerationData[[2]] > 0]), 0)
  negativeAccelerationSd <- ifelse(sum(accelerationData[[2]] < 0) > 1, sd(accelerationData[[2]][accelerationData[[2]] < 0]), 0)
  #Time spent stopping
  timeStopped <- sum(rawVelocities < 0.5)
  
  #Coordinates as ltraj object
  tripLtraj <- as.ltraj(trip, id = rep(file, nrow(trip)), typeII = FALSE)
  tripLtrajDf <- ld(tripLtraj)
  
  #Bad Data Detector
  badData <- ifelse(speedSd < 5, 1, 0)
  #rawDist <- quantile(rawVelocities, seq(0.1, 1, by=0.1)) #debug only
  
  #Distances
  distances <- na.omit(tripLtrajDf$dist)
  #Distances outlier removal
  distances[distances > mean(distances) + sd(distances) * 4] <- NA
  distances[is.na(distances)] <- mean(distances, na.rm = TRUE)
  distanceTrip <- sum(distances)
  #partition the trajectory using the method of Lavielle
  lav <- lavielle(tripLtraj, Lmin = 2, Kmax = 20)
  chooseLav <- chooseseg(lav, draw = FALSE)
  trajectoryPartitions <- tail(which(abs(chooseLav$D) > 0.70), 1)
  #kk <- findpath(lav, trajectoryPartitions) #here for debugging purposes only
  #plot(kk) #here for debugging purposes only

  #Angles
  #Using old arc cosine function
  #anglesData <- angleVector(trip)
  #turningAnglesDist <- anglesData[[1]]
  #turningAnglesSd <- sd(anglesData[[3]])    
  #Absolute angles using the adehabitatLT package
  tripLtrajDf$abs.angle[is.na(tripLtrajDf$abs.angle)] <- 0
  absAnglesData <- quantSigma(tripLtrajDf$abs.angle, returnVector = TRUE)  
  absAnglesDist <- absAnglesData[[1]]
  absAnglesSd <- sd(absAnglesData[[2]])
  absAnglesPositiveData <- quantSigma(abs(tripLtrajDf$abs.angle), returnVector = TRUE)   
  absAnglesPositiveDist <- absAnglesPositiveData[[1]]
  absAnglesPositiveSd <- sd(absAnglesPositiveData[[2]])
  
  #Relative angles using the adehabitatLT package
  tripLtrajDf$rel.angle[is.na(tripLtrajDf$rel.angle)] <- 0
  relAnglesData <- quantSigma(tripLtrajDf$rel.angle, returnVector = TRUE)  
  relAnglesDist <- relAnglesData[[1]]
  relAnglesSd <- sd(relAnglesData[[2]])
  relAnglesPositiveData <- quantSigma(abs(tripLtrajDf$rel.angle), returnVector = TRUE)   
  relAnglesPositiveDist <- relAnglesPositiveData[[1]]
  relAnglesPositiveSd <- sd(relAnglesPositiveData[[2]])
  
  #Relative Angles times speed
  anglesTimesSpeedData <- quantSigma(tripLtrajDf$rel.angle[c(-1, -(length(tripLtrajDf$rel.angle) - 1))] 
                                     * rawVelocities[-1], returnVector = TRUE)
  anglesTimesSpeedDist <- anglesTimesSpeedData[[1]]
  anglesTimesSpeedSd <- sd(anglesTimesSpeedData[[2]])
  
  #Positive Relative Angles times speed
  posAnglesTimesSpeedData <- quantSigma(abs(tripLtrajDf$rel.angle)[c(-1, -(length(tripLtrajDf$rel.angle) - 1))] 
                                     * rawVelocities[-1], returnVector = TRUE)
  posAnglesTimesSpeedDist <- posAnglesTimesSpeedData[[1]]
  posAnglesTimesSpeedSd <- sd(posAnglesTimesSpeedData[[2]])
  
  #Generalized Procrustes analysis of points / shapes
  #GPATrip <- gpagen(trip)  
  
  return(c(speedDist, speedSd, speedDistWOStopsDist, speedDistWOStopsSd, accelerationsDist, positiveAccelerationSd, 
           negativeAccelerationSd, nrow(trip), timeStopped, distanceTrip, trajectoryPartitions,
           absAnglesDist, absAnglesSd, absAnglesPositiveDist, absAnglesPositiveSd,
           relAnglesDist, relAnglesSd, relAnglesPositiveDist, relAnglesPositiveSd,
           anglesTimesSpeedDist, anglesTimesSpeedSd, posAnglesTimesSpeedDist, posAnglesTimesSpeedSd, 
           badData))
           #badData, sd(rawVelocities), rawDist)) #BAD DATA DEBUG ONLY
}

#EDA----------------------------------------
##EDA Pt. 1 Visualization of Trajectories
driverViz <- sample(drivers, 1)
results <- mclapply(seq(1, 200), function(file, driverID){
  tripCoordinates <- fread(file.path(driversDirectory, driverID, paste0(file, ".csv")))
  return(cbind(tripCoordinates, rep(file, nrow(tripCoordinates))))
}, mc.cores = numCores, driverID = driverViz)
print(paste0("Driver number: ", driverViz, " processed"))
dataTableReady2Plot <- do.call(rbind, results)
qplot(x, y, data = dataTableReady2Plot, colour = V2, geom = "point")
#py <- plotly()
#py$ggplotly()

##EDA Pt. 2 Visualization of Speeds with and without outlier replacement
driverViz <- sample(drivers, 1)
fileViz <- sample(1:200, 1)
tripCoordinates <- fread(file.path(driversDirectory, driverViz, paste0(fileViz, ".csv")))
ggplot(as.data.frame(tripCoordinates), aes(x = x, y = y)) + geom_point()
print(paste0("Driver number: ", driverViz, " trip number ", fileViz, " processed"))
speed2Plot <- 3.6 * sqrt(diff(tripCoordinates$x)^2 + diff(tripCoordinates$y)^2)
ggplot(as.data.frame(speed2Plot), aes(x = speed2Plot)) + geom_density()
#Remove data above 5 sigmas (5 standard deviations)
speed2Plot2 <- speed2Plot[!speed2Plot > mean(speed2Plot) + sd(speed2Plot) * 5]
ggplot(as.data.frame(speed2Plot2), aes(x = speed2Plot2)) + geom_density()

##EDA Pt. 3 Visualization of Turning Angles with and without outlier replacement
driverViz <- sample(drivers, 1)
fileViz <- sample(1:200, 1)
tripCoordinates <- fread(file.path(driversDirectory, driverViz, paste0(fileViz, ".csv")))
ggplot(as.data.frame(tripCoordinates), aes(x = x, y = y)) + geom_point()
print(paste0("Driver number: ", driverViz, " trip number ", fileViz, " processed"))
dir <- as.data.table(cbind(diff(tripCoordinates$x), diff(tripCoordinates$y)))
setnames(dir, old = c("V1", "V2"), new = c("x", "y"))
dir2 <- dir[2:nrow(dir)]
dir1 <- dir[1:(nrow(dir) - 1)]
angles <- acos(rowSums(dir1 * dir2) / ((sqrt(rowSums(dir1^ 2) * rowSums(dir2 ^2))) + 1e-06)) * (180/pi) #1e-06 is included to avoid a division by zero
ggplot(as.data.frame(angles), aes(x = angles)) + geom_density()
#Remove data above 4 sigmas (4 standard deviations)
#NAs removal
anglesWoOutliers <- na.omit(angles)
#Remove orthogonal values
anglesWoOutliers[anglesWoOutliers >= 90] <- 0
#Outliers Removal
anglesWoOutliers <- anglesWoOutliers[!anglesWoOutliers > mean(anglesWoOutliers) + sd(anglesWoOutliers) * 4]
ggplot(as.data.frame(anglesWoOutliers), aes(x = anglesWoOutliers)) + geom_density()

##EDA Pt. 4 LOF Algorithm visualization
driverViz <- sample(drivers, 1)
results <- unlist(mclapply(seq(1, 200), transform2Percentiles, mc.cores = numCores, driverID = driverViz))
results <- scale(matrix(results, nrow = 200, byrow = TRUE))
#Plot
par(mfrow=c(1, 2))
plot(density(lofactor(results, k=5)))
biplot(prcomp(results), cex=.8)

## EDA Pt. 5 Bad Data Detector
driverViz <- sample(drivers, 1)
fileViz <- sample(1:200, 1)
#Bad Data Detector
is.badData <- function(theTrip){
  #Velocities
  rawVelocitiesT <- velocitiesKH(theTrip)
  #n-sigma removal
  velocityDataR <- quantSigma(rawVelocitiesT)
  speedSdR <- ifelse(sd(velocityDataR[[2]]) < 5, 1, 0)
  return(speedSdR)
}

par(mfrow=c(2, 2))
randTrip <- fread(file.path(driversDirectory, driverViz, paste0(fileViz, ".csv")))
print(is.badData(randTrip))
plot(randTrip, main = "Random Trip")

#Known bad data
badTrip <- fread(file.path(driversDirectory, 1, paste0(53, ".csv")))
print(is.badData(badTrip))
plot(badTrip, main = "Bad Trip")

#known split good data
splitTrip <- fread(file.path(driversDirectory, 1, paste0(136, ".csv")))
print(is.badData(splitTrip))
plot(splitTrip, main = "Split Trip")

#good data with lots of stationary points
stillTrip <- fread(file.path(driversDirectory, 1, paste0(179, ".csv")))
print(is.badData(stillTrip))
plot(stillTrip, main = "Still Trip")

## EDA Pt. 6 Determine the minimal PCAs / number of neurons in the middle layer
#Begin with a randomly selected driver to start the PCA calculation
numberOfDrivers <- 100
initialDrivers <- sample(drivers, numberOfDrivers)
driversProcessed <- sapply(initialDrivers, function(driver){
  results <- unlist(mclapply(seq(1, 200), transform2Percentiles, mc.cores = numCores, driverID = driver))
  print(paste0("Driver number: ", driver, " processed"))
  return(results)  
})
driversProcessed <- signif(matrix(unlist(driversProcessed), nrow = numberOfDrivers * 200, byrow = TRUE), digits = 4)
#Bad Data Removal
driversProcessed <- driversProcessed[!driversProcessed[, ncol(driversProcessed)] == 1, -ncol(driversProcessed)]
#Scale Data
#driversProcessed <- scale(driversProcessed)

#Init h2o Server
#Start h2o directly from R
h2oServer <- h2o.init(ip = "localhost", max_mem_size = "5g", nthreads = -1)
h2oResult <- as.h2o(h2oServer, driversProcessed)
print(h2o.ls(h2oServer))
rm(driversProcessed)

#PCA
PCAModel <- h2o.prcomp(h2oResult)
plot(PCAModel@model$sdev)
#ggplot(data.frame(X = PCAModel@model$sdev), aes(x = X)) + geom_density()
h2o.shutdown(h2oServer, prompt = FALSE)

## EDA Pt. 7 Determine the best number of drivers and trips per driver using random forests (fastest algorithm available)
#Begin creating a grid
set.seed(10014)
driverGridVal <- sample(drivers, 15)
RFGrid <- expand.grid(.numberOfNegativeDrivers = c(3, 5, 10, 30, 60),
                      .numberOfNegativeRows = c(50, 200, 650, 1500, 3000))
RFGrid <-  RFGrid[c(-11, -16, -17, -21, -22, -23), ] #Remove row where sample is larger than the population when 'replace = FALSE' 

#Start h2o directly from R
h2oServer <- h2o.init(ip = "localhost", max_mem_size = "5g", nthreads = -1)

modelsGenerated <- apply(RFGrid, 1, function(driversSplit){  
  errorNthDriver <- sapply(driverGridVal, function(nthDriver){
    #Parallel processing of each driver data
    resultsFull <- unlist(mclapply(seq(1, 200), transform2Percentiles, mc.cores = numCores, driverID = nthDriver))
    resultsFull <- signif(matrix(resultsFull, nrow = 200, byrow = TRUE), digits = 3)
    #Bad Data Removal
    idxBadData <- resultsFull[, ncol(resultsFull)] == 1
    results <- resultsFull[!resultsFull[, ncol(resultsFull)] == 1, -ncol(resultsFull)]
    print(paste0("Driver number ", nthDriver, " processed"))  
    
    #Sample data from other drivers  
    numberOfDrivers <- as.numeric(driversSplit[1])
    initialDrivers <- sample(drivers[!drivers %in% nthDriver], numberOfDrivers)
    ExtraDrivers <- sapply(initialDrivers, function(driverExtra){
      results <- unlist(mclapply(seq(1, 200), transform2Percentiles, mc.cores = numCores, driverID = driverExtra))
      print(paste0("Driver number: ", driverExtra, " processed"))
      return(results)
    })
    ExtraDrivers <- signif(matrix(unlist(ExtraDrivers), nrow = numberOfDrivers * 200, byrow = TRUE), digits = 3)
    #Bad Data Removal
    ExtraDrivers <- ExtraDrivers[!ExtraDrivers[, ncol(ExtraDrivers)] == 1, -ncol(ExtraDrivers)]
    #Extra drivers sampling
    ExtraDrivers <- ExtraDrivers[sample(seq(1, nrow(ExtraDrivers)), as.numeric(driversSplit[2])), ]  
    
    #R matrix conversion to h2o object and stored in the server
    h2oResultPlusExtras <- as.h2o(h2oServer, cbind(c(rep(1, nrow(results)), rep(0, nrow(ExtraDrivers))), 
                                                   signif(rbind(results, ExtraDrivers), digits = 4)))
    print(h2o.ls(h2oServer))   
    
    #h2o.ai RF algorithm
    #Shuffle indexes
    #set.seed(1001001)
    randIdxs <- sample(seq(1, nrow(results) + nrow(ExtraDrivers)), nrow(results) + nrow(ExtraDrivers))
    #Cross Validation + Modelling
    driverRFModelCV <- h2o.randomForest(x = seq(2, ncol(h2oResultPlusExtras)), y = 1,
                                        data = h2oResultPlusExtras[randIdxs, ],
                                        nfolds = 5,
                                        classification = TRUE,
                                        ntree = c(50, 75, 100),
                                        depth = c(20, 50, 75), 
                                        verbose = TRUE)
    aucError <- driverRFModelCV@model[[1]]@model$auc
    h2o.rm(object = h2oServer, keys = h2o.ls(h2oServer)[, 1])       
    return(aucError)
  })
  return(c(as.numeric(driversSplit), errorNthDriver))
})
#h2o Shutdown
h2o.shutdown(h2oServer, prompt = FALSE)

modelsGenerated <- t(modelsGenerated)

averageAUC <- apply(modelsGenerated[, c(-1, -2)], 1, mean)
bestDriverNumbers <- which.max(averageAUC)
driversParameters <- modelsGenerated[bestDriverNumbers, c(1, 2)]
  
#Pre-Train Model using h2o.ai deeplearning


#Modelling---------------------------
#load(file.path(workingDirectory, "deepNetPath.RData"))

#Interrupting Remote Instance (in case this is running on an AWS spot instance)
SpotInstance <- TRUE
if (SpotInstance == TRUE){
  driversProcessed <- gsub(pattern = ".csv", replacement = "", x = list.files(outputDirectory))
  drivers <- c(driversProcessed[length(driversProcessed)],
               list.files(driversDirectory)[!(list.files(driversDirectory) %in% driversProcessed)])
}

#Init h2o Server
#Start from R
#h2oServer <- h2o.init(ip = "localhost", port = 54321, max_mem_size = '13g', startH2O = TRUE, nthreads = -1)

#Start h2o from command line
system(paste0("java -Xmx20G -jar ", h2o.jarLoc, " -port 54333 -name AXA &"))
#Connect R to h2o
h2oServer <- h2o.init(ip = "localhost", port = 54333, nthreads = -1)  
#checkpointModelKey <- ""
#checkpointModel <- ""

if ("deepNetPath" %in% ls()){
  #Load pretrained model
  checkpointModel <- h2o.loadModel(h2oServer, deepNetPath)
  checkpointModelKey <- h2o.ls(h2oServer)[, 1]  
  print(h2o.ls(h2oServer))  
}

driversPredictions <- lapply(drivers, function(driver){
  #Parallel processing of each driver data
  resultsFull <- unlist(mclapply(seq(1, 200), transform2Percentiles, mc.cores = numCores, driverID = driver))
  resultsFull <- signif(matrix(resultsFull, nrow = 200, byrow = TRUE), digits = 3)
  #Bad Data Removal
  idxBadData <- resultsFull[, ncol(resultsFull)] == 1
  results <- resultsFull[!resultsFull[, ncol(resultsFull)] == 1, -ncol(resultsFull)]
  print(paste0("Driver number ", driver, " processed"))  
    
  #Sample data from other drivers  
  numberOfDrivers <- 3
  initialDrivers <- sample(drivers[!drivers %in% driver], numberOfDrivers)
  ExtraDrivers <- sapply(initialDrivers, function(driver){
    results <- unlist(mclapply(seq(1, 200), transform2Percentiles, mc.cores = numCores, driverID = driver))
    print(paste0("Driver number: ", driver, " processed"))
    return(results)
  })
  ExtraDrivers <- signif(matrix(unlist(ExtraDrivers), nrow = numberOfDrivers * 200, byrow = TRUE), digits = 3)
  #Bad Data Removal
  ExtraDrivers <- ExtraDrivers[!ExtraDrivers[, ncol(ExtraDrivers)] == 1, -ncol(ExtraDrivers)]
  #Extra drivers sampling
  ExtraDrivers <- ExtraDrivers[sample(seq(1, nrow(ExtraDrivers)), 100), ]  
  
  #LOF Algorithm
  lofDriver <- signif(-(lofactor(resultsFull, k=10)) + 10, digits = 4)
  lofDriverRanking <- rank(lofDriver)
  print(paste0("Driver number ", driver, " processed with the LOF Algorithm"))  
  
  #Shuffle indexes
  #set.seed(1001001)
  randIdxs <- sample(seq(1, nrow(results) + nrow(ExtraDrivers)), nrow(results) + nrow(ExtraDrivers))
  
  #R matrix conversion to h2o object and stored in the server
  h2oResultPlusExtras <- as.h2o(h2oServer, cbind(c(rep(1, nrow(results)), rep(0, nrow(ExtraDrivers))), 
                                                 signif(rbind(results, ExtraDrivers), digits = 4)))
  h2oResultsNthDriver <- as.h2o(h2oServer, cbind(rep(1, nrow(resultsFull)),
                                                 signif(resultsFull[, -ncol(resultsFull)], digits = 4)))
  print(h2o.ls(h2oServer))
  
  #h2o.ai RF algorithm
  #Cross Validation + Modelling
  driverRFModelCV <- h2o.randomForest(x = seq(2, ncol(h2oResultPlusExtras)), y = 1,
                                      data = h2oResultPlusExtras[randIdxs, ],
                                      nfolds = 5,
                                      classification = TRUE,
                                      ntree = c(50, 75, 100, 150),
                                      depth = c(20, 50, 75), 
                                      verbose = TRUE)
  #Log Info
  aucRF <- driverRFModelCV@model[[1]]@model$auc
  ntreeRF <- driverRFModelCV@model[[1]]@model$params$ntree
  depthRF <- driverRFModelCV@model[[1]]@model$params$depth
  
  #Best Model
  bestCVRF <- driverRFModelCV@model[[1]]
  
  driverRFModel <- h2o.randomForest(x = seq(2, ncol(h2oResultPlusExtras)), y = 1,
                                    data = h2oResultPlusExtras[randIdxs, ],
                                    classification = TRUE,
                                    type = "BigData",
                                    ntree = driverRFModelCV@model[[1]]@model$params$ntree,
                                    depth = driverRFModelCV@model[[1]]@model$params$depth, 
                                    verbose = TRUE)
  
  print(driverRFModel)
  
  #probability Prediction of trips in Nth driver 
  predictionRF <- signif(as.data.frame(h2o.predict(driverRFModel, newdata = h2oResultsNthDriver)[, 3]), digits = 4)
  #Bad Data rounding up to one
  predictionRF[idxBadData, 1] <- 1
  predictionRFRank <- rank(predictionRF[, 1])
  
  predictionRFCV <- signif(as.data.frame(h2o.predict(bestCVRF, newdata = h2oResultsNthDriver)[, 3]), digits = 4)
  #Bad Data rounding up to one
  predictionRFCV[idxBadData, 1] <- 1
  predictionRFCVRank <- rank(predictionRFCV[, 1])
  
  print(h2o.ls(h2oServer))  
  h2o.rm(object = h2oServer, keys = h2o.ls(h2oServer)[, 1])   
  print(paste0("Driver number ", driver, " processed with RFs")) 
  
  #h2o.ai GBM algorithm
  #R matrix conversion to h2o object and stored in the server
  h2oResultPlusExtras <- as.h2o(h2oServer, cbind(c(rep(1, nrow(results)), rep(0, nrow(ExtraDrivers))), 
                                                 signif(rbind(results, ExtraDrivers), digits = 4)))
  h2oResultsNthDriver <- as.h2o(h2oServer, cbind(rep(1, nrow(resultsFull)),
                                                 signif(resultsFull[, -ncol(resultsFull)], digits = 4)))
  print(h2o.ls(h2oServer))
  
  #Cross Validation
  driverGBMModelCV <-  h2o.gbm(x = seq(2, ncol(h2oResultPlusExtras)), y = 1,
                               data = h2oResultPlusExtras[randIdxs, ],
                               nfolds = 5,
                               distribution = "bernoulli",
                               interaction.depth = c(2, 4, 7),
                               shrinkage = c(0.001, 0.003), 
                               n.trees = 150)
  #Log Info  
  aucGBM <- driverGBMModelCV@model[[1]]@model$auc
  interaction.depthGBM <- driverGBMModelCV@model[[1]]@model$params$interaction.depth
  shrinkageGBM <- driverGBMModelCV@model[[1]]@model$params$shrinkage
  
  #Best Model
  bestCVGBM <- driverGBMModelCV@model[[1]]
  
  #Modelling
  driverGBMModel <-  h2o.gbm(x = seq(2, ncol(h2oResultPlusExtras)), y = 1,
                             data = h2oResultPlusExtras[randIdxs, ],
                             distribution = "bernoulli",
                             interaction.depth = driverGBMModelCV@model[[1]]@model$params$interaction.depth,
                             shrinkage = driverGBMModelCV@model[[1]]@model$params$shrinkage, 
                             n.trees = 3000)  
  
  print(driverGBMModel)
  
  #probability Prediction of trips in Nth driver 
  predictionGBM <- signif(as.data.frame(h2o.predict(driverGBMModel, newdata = h2oResultsNthDriver)[, 3]), digits = 4)
  #Bad Data rounding up to one
  predictionGBM[idxBadData, 1] <- 1
  predictionGBMRank <- rank(predictionGBM[, 1])
  
  predictionGBMCV <- signif(as.data.frame(h2o.predict(bestCVGBM, newdata = h2oResultsNthDriver)[, 3]), digits = 4)
  #Bad Data rounding up to one
  predictionGBMCV[idxBadData, 1] <- 1
  predictionGBMCVRank <- rank(predictionGBMCV[, 1])
  
  print(h2o.ls(h2oServer))  
  h2o.rm(object = h2oServer, keys = h2o.ls(h2oServer)[, 1])     
  print(paste0("Driver number ", driver, " processed with GBMs"))  
  
  #h20.ai deep learning algorithm 
  #R matrix conversion to h2o object and stored in the server
  h2oResultPlusExtras <- as.h2o(h2oServer, cbind(c(rep(1, nrow(results)), rep(0, nrow(ExtraDrivers))), 
                                                 signif(rbind(results, ExtraDrivers), digits = 4)))
  h2oResultsNthDriver <- as.h2o(h2oServer, cbind(rep(1, nrow(resultsFull)),
                                                 signif(resultsFull[, -ncol(resultsFull)], digits = 4)))
  print(h2o.ls(h2oServer))
  
  #Cross Validation 
  driverDeepNNModelCV <- h2o.deeplearning(x = seq(2, ncol(h2oResultPlusExtras)), y = 1,
                                          data = h2oResultPlusExtras[randIdxs, ],   
                                          nfolds = 5,
                                          classification = TRUE,
                                          #checkpoint = ifelse("deepNetPath" %in% ls(), checkpointModel, ""),
                                          activation = c("Tanh", "TanhWithDropout"),
                                          input_dropout_ratio = c(0, 0.2),
                                          hidden_dropout_ratio = list(c(0, 0, 0, 0, 0), c(0.5, 0.5, 0.5, 0.5, 0.5)),
                                          l1 = c(0, 1e-5),
                                          l2 = c(0, 1e-5),
                                          rho = c(0.95, 0.99),
                                          epsilon = c(1e-12, 1e-10, 1e-08),
                                          hidden = c(100, 100, 100, 100, 100), 
                                          epochs = 60)  
  #Log Info  
  aucNN <- driverDeepNNModelCV@model[[1]]@model$auc
  activationNN <- driverDeepNNModelCV@model[[1]]@model$params$activation
  l1NN <- driverDeepNNModelCV@model[[1]]@model$params$l1
  l2NN <- driverDeepNNModelCV@model[[1]]@model$params$l2
  rhoNN <- driverDeepNNModelCV@model[[1]]@model$params$rho
  epsilonNN <- driverDeepNNModelCV@model[[1]]@model$params$epsilon
  
  #Best Model
  bestCVNN <- driverDeepNNModelCV@model[[1]]
  
  #Select the best 2 models and keep training them  
  driverDeepNNModel <- h2o.deeplearning(x = seq(2, ncol(h2oResultPlusExtras)), y = 1,
                                        data = h2oResultPlusExtras[randIdxs, ],   
                                        classification = TRUE,
                                        #checkpoint = ifelse("deepNetPath" %in% ls(), checkpointModel, ""),
                                        activation = driverDeepNNModelCV@model[[1]]@model$params$activation,
                                        input_dropout_ratio = driverDeepNNModelCV@model[[1]]@model$params$input_dropout_ratio,
                                        hidden_dropout_ratio = driverDeepNNModelCV@model[[1]]@model$params$hidden_dropout_ratio,                                        
                                        l1 = driverDeepNNModelCV@model[[1]]@model$params$l1,
                                        l2 = driverDeepNNModelCV@model[[1]]@model$params$l2,
                                        rho = driverDeepNNModelCV@model[[1]]@model$params$rho,
                                        epsilon = driverDeepNNModelCV@model[[1]]@model$params$epsilon,
                                        hidden = c(100, 100, 100, 100, 100), 
                                        epochs = 250)
  
  driverDeepNNModel2 <- h2o.deeplearning(x = seq(2, ncol(h2oResultPlusExtras)), y = 1,
                                         data = h2oResultPlusExtras[randIdxs, ],   
                                         classification = TRUE,
                                         #checkpoint = ifelse("deepNetPath" %in% ls(), checkpointModel, ""),
                                         activation = driverDeepNNModelCV@model[[2]]@model$params$activation,
                                         input_dropout_ratio = driverDeepNNModelCV@model[[2]]@model$params$input_dropout_ratio,
                                         hidden_dropout_ratio = driverDeepNNModelCV@model[[2]]@model$params$hidden_dropout_ratio,                                                                                 
                                         l1 = driverDeepNNModelCV@model[[2]]@model$params$l1,
                                         l2 = driverDeepNNModelCV@model[[2]]@model$params$l2,
                                         rho = driverDeepNNModelCV@model[[2]]@model$params$rho,
                                         epsilon = driverDeepNNModelCV@model[[2]]@model$params$epsilon,
                                         hidden = c(100, 100, 100, 100, 100), 
                                         epochs = 250)
  
  print(driverDeepNNModel)
  print(driverDeepNNModel2)

  #probability Predictions on all trips in Nth driver 
  predictionNN <- signif(as.data.frame(h2o.predict(driverDeepNNModel, newdata = h2oResultsNthDriver)[, 3]), digits = 8)
  #Bad Data rounding up to one
  predictionNN[idxBadData, 1] <- 1
  predictionNNRank <- rank(predictionNN[, 1])
  
  predictionNN2 <- signif(as.data.frame(h2o.predict(driverDeepNNModel2, newdata = h2oResultsNthDriver)[, 3]), digits = 8)
  #Bad Data rounding up to one
  predictionNN2[idxBadData, 1] <- 1
  predictionNNRank2 <- rank(predictionNN2[, 1])
  
  predictionNNCV <- signif(as.data.frame(h2o.predict(bestCVNN, newdata = h2oResultsNthDriver)[, 3]), digits = 8)
  #Bad Data rounding up to one
  predictionNNCV[idxBadData, 1] <- 1
  predictionNNCVRank <- rank(predictionNNCV[, 1])
  
  print(h2o.ls(h2oServer))  
  h2o.rm(object = h2oServer, keys = h2o.ls(h2oServer)[, 1])   
  print(paste0("Driver number ", driver, " processed with Deep NNs"))  
  
  #Logging hyperparameters
  write.csv(cbind(aucRF, ntreeRF, depthRF,
                  aucGBM, interaction.depthGBM, shrinkageGBM,
                  aucNN, activationNN, l1NN, l2NN, rhoNN, epsilonNN), 
            file = file.path(logsDirectory, paste0(driver, ".csv")), row.names = FALSE)
  
  #h2oObjects2Remove <- which(!h2o.ls(h2oServer)[, 1] %in% checkpointModelKey)
  #h2o.rm(object = h2oServer, keys = h2o.ls(h2oServer)[h2oObjects2Remove, 1]) 
  print(paste0(which(list.files(driversDirectory) == driver), "/", length(list.files(driversDirectory))))
  
  if (SpotInstance == TRUE){
    write.csv(cbind(lofDriver, lofDriverRanking,
                    predictionRF[, 1], predictionRFRank, 
                    predictionGBM[, 1], predictionGBMRank, 
                    predictionNN[, 1], predictionNNRank, 
                    predictionNN2[, 1], predictionNNRank2, 
                    predictionRFCV[, 1], predictionRFCVRank,
                    predictionGBMCV[, 1], predictionGBMCVRank,
                    predictionNNCV[, 1], predictionNNCVRank), 
              file = file.path(outputDirectory, paste0(driver, ".csv")), row.names = FALSE)
    return(TRUE)
    
  }else{
    return(cbind(lofDriver, lofDriverRanking,
                 predictionRF[, 1], predictionRFRank, 
                 predictionGBM[, 1], predictionGBMRank, 
                 predictionNN[, 1], predictionNNRank,
                 predictionNN2[, 1], predictionNNRank2, 
                 predictionNNRF[, 1], predictionRFCVRank,
                 predictionGBMCV[, 1], predictionGBMCVRank,
                 predictionNNCV[, 1], predictionNNCVRank))
  }  
})
#Shutdown h20 instance
h2o.shutdown(h2oServer, prompt = FALSE)

if (SpotInstance == TRUE){
  #List all possible drivers identities
  driversOutput <- list.files(outputDirectory)   
  driversPredictions <- lapply(driversOutput, function(driver){
    predictions <- fread(file.path(outputDirectory, driver), header = TRUE,
                         stringsAsFactors = FALSE)
    print(paste0("Driver number ", driver, " read"))
    return(predictions)
  })
}

#Concatenate lists into a data.frame
driversPredictions <- do.call(rbind, driversPredictions)

#Write .csv files-------------------------
submissionTemplate <- fread(file.path(otherDataDirectory, "sampleSubmission.csv"), header = TRUE,
                            stringsAsFactors = FALSE, colClasses = c("character", "numeric"))

#Probabilities h2o.ai Deep NN
submissionTemplate$prob <- signif(driversPredictions[, 1], digits = 4)
write.csv(submissionTemplate, file = "lofScoreI.csv", row.names = FALSE)
system('zip lofScoreI.zip lofScoreI.csv')

#Rank h2o.ai Deep NN
submissionTemplate$prob <- signif(driversPredictions[, 2], digits = 4)
write.csv(submissionTemplate, file = "lofRankI.csv", row.names = FALSE)
system('zip lofRankI.zip lofRankI.csv')

#Probabilities h2o.ai RF
submissionTemplate$prob <- signif(driversPredictions[, 3], digits = 4)
write.csv(submissionTemplate, file = "RFProbI.csv", row.names = FALSE)
system('zip RFProbI.zip RFProbI.csv')

#Rank h2o.ai RF
submissionTemplate$prob <- signif(driversPredictions[, 4], digits = 4)
write.csv(submissionTemplate, file = "RFRankI.csv", row.names = FALSE)
system('zip RFRankI.zip RFRankI.csv')

#Probabilities h2o.ai GBM
submissionTemplate$prob <- signif(driversPredictions[, 5], digits = 4)
write.csv(submissionTemplate, file = "GBMProbI.csv", row.names = FALSE)
system('zip GBMProbI.zip GBMProbI.csv')

#Rank h2o.ai GBM
submissionTemplate$prob <- signif(driversPredictions[, 6], digits = 4)
write.csv(submissionTemplate, file = "GBMRankI.csv", row.names = FALSE)
system('zip GBMRankI.zip GBMRankI.csv')

#Probabilities h2o.ai Deep NN
submissionTemplate$prob <- signif(driversPredictions[, 7], digits = 4)
write.csv(submissionTemplate, file = "deepNNProbI.csv", row.names = FALSE)
system('zip deepNNProbI.zip deepNNProbI.csv')

#Rank h2o.ai Deep NN
submissionTemplate$prob <- signif(driversPredictions[, 8], digits = 4)
write.csv(submissionTemplate, file = "deepNNRankI.csv", row.names = FALSE)
system('zip deepNNRankI.zip deepNNRankI.csv')

#Probabilities h2o.ai Deep NN second best model
submissionTemplate$prob <- signif(driversPredictions[, 9], digits = 4)
write.csv(submissionTemplate, file = "deepNNProb2I.csv", row.names = FALSE)
system('zip deepNNProb2I.zip deepNNProb2I.csv')

#Rank h2o.ai Deep NN second best model
submissionTemplate$prob <- signif(driversPredictions[, 10], digits = 4)
write.csv(submissionTemplate, file = "deepNNRank2I.csv", row.names = FALSE)
system('zip deepNNRank2I.zip deepNNRank2I.csv')

#Probabilities h2o.ai best RF CV Model
submissionTemplate$prob <- signif(driversPredictions[, 11], digits = 4)
write.csv(submissionTemplate, file = "RFProbsCVI.csv", row.names = FALSE)
system('zip RFProbsCVI.zip RFProbsCVI.csv')

#Rank h2o.ai best RF CV Model
submissionTemplate$prob <- signif(driversPredictions[, 12], digits = 4)
write.csv(submissionTemplate, file = "RFRankCVI.csv", row.names = FALSE)
system('zip RFRankCVI.zip RFRankCVI.csv')

#Probabilities best GBM CV Model
submissionTemplate$prob <- signif(driversPredictions[, 13], digits = 4)
write.csv(submissionTemplate, file = "GBMProbsCVI.csv", row.names = FALSE)
system('zip GBMProbsCVI.zip GBMProbsCVI.csv')

#Rank best GBM CV Model
submissionTemplate$prob <- signif(driversPredictions[, 14], digits = 4)
write.csv(submissionTemplate, file = "GBMRankCVI.csv", row.names = FALSE)
system('zip GBMRankCVI.zip GBMRankCVI.csv')

#Probabilities best NN CV Model
submissionTemplate$prob <- signif(driversPredictions[, 15], digits = 4)
write.csv(submissionTemplate, file = "NNProbsCVI.csv", row.names = FALSE)
system('zip NNProbsCVI.zip NNProbsCVI.csv')

#Probabilities best NN CV Model
submissionTemplate$prob <- signif(driversPredictions[, 16], digits = 4)
write.csv(submissionTemplate, file = "NNRankCVI.csv", row.names = FALSE)
system('zip NNRankCVI.zip NNRankCVI.csv')