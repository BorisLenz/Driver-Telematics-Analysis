#AXA Driver Telematics Analysis
#Ver 0.9.6 # Updated code to include trip matching/boosting and tripMatchingFun integrated

#Init-----------------------------------------------
rm(list=ls(all=TRUE))

#Libraries, directories, options and extra functions----------------------
require("data.table")
require("miscTools")
require("parallel")
require('doParallel')
require("h2o")
require("DMwR")
require("prospectr")
require("adehabitatLT")
require("ggplot2")
#require("plotly")

Set Working Directory
workingDirectory <- "/home/wacax/Wacax/Kaggle/AXA-Driver-Telematics-Analysis/"
setwd(workingDirectory)
driversDirectory <- "/home/wacax/Wacax/Kaggle/AXA-Driver-Telematics-Analysis/Data/drivers"
otherDataDirectory <- "/home/wacax/Wacax/Kaggle/AXA-Driver-Telematics-Analysis/Data/"
outputDirectory <- "/home/wacax/Wacax/Kaggle/AXA-Driver-Telematics-Analysis/Data/Output"
logsDirectory <- "/home/wacax/Wacax/Kaggle/AXA-Driver-Telematics-Analysis/Data/Logs"
tripMatchingDir <- "/home/wacax/Wacax/Kaggle/AXA-Driver-Telematics-Analysis/Data/TripMatches"

vw77Dir = "/home/wacax/vowpal_wabbit-7.7/vowpalwabbit/"
#h2o location
h2o.jarLoc <- "/home/wacax/R/x86_64-pc-linux-gnu-library/3.1/h2o/java/h2o.jar"

#List all possible drivers identities
drivers <- list.files(driversDirectory)
#Detect available cores
numCores <- detectCores()

#Extra Functions
source(paste0(workingDirectory, "lofAnomalyDetection.R"))
source(paste0(workingDirectory, "AutoencodersAnomalyDetection.R"))
source(paste0(workingDirectory, "TripMatchingFun.R"))

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

#Define function to be passed as parallel
transform2Percentiles <- function(file, driverID){  
  trip <- fread(file.path(driversDirectory, driverID, paste0(file, ".csv")))
  #Velocities
  rawVelocities <- velocitiesKH(trip)
  #n-sigma removal
  velocityData <- quantSigma(rawVelocities)
  speedDist <- velocityData[[1]]
  speedSd <- sd(velocityData[[2]])
  #speedMad <- mad(velocityData[[2]])
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

#Start h2o directly from R
h2oServer <- h2o.init(ip = "localhost", max_mem_size = "5g", nthreads = -1)
#Parallel processing of te first driver
resultsDriver1 <- unlist(mclapply(seq(1, 200), transform2Percentiles, mc.cores = numCores, driverID = 1))
resultsDriver1 <- signif(matrix(resultsDriver1, nrow = 200, byrow = TRUE), digits = 3)
#Real bad data known for driver #1
badDataDriver1Idxs <-  c(14, 21, 28, 39, 53, 55, 65, 90, 106, 114, 121, 123, 133, 146, 155, 159, 161, 169, 177, 186)
#Rand indexes
randIdxs <- sample(seq(1, nrow(resultsDriver1)), nrow(resultsDriver1))

#R matrix conversion to h2o object and stored in the server
h2oResultsDriver1 <- as.h2o(h2oServer, cbind(rbinom(length(randIdxs), 1, 0.5),
                                             signif(resultsDriver1, digits = 4)))

#Deep encoders outlier/bad data detection
ae_model <- h2o.deeplearning(x = seq(2, ncol(h2oResultsDriver1) - 1), y = 1,
                             data = h2oResultsDriver1,   
                             autoencoder = TRUE, 
                             fast_mode = TRUE,
                             ignore_const_cols = FALSE,
                             activation = "Rectifier",
                             l1 = 0,
                             l2 = 0,
                             rho = 0.99,
                             epsilon = 1e-10,
                             hidden = c(200, 200), 
                             epochs = 300)

#Anomaly Detection
test_rec_error <- as.data.frame(h2o.anomaly(h2oResultsDriver1, ae_model))
#Check Deep Features
test_features_deep <- h2o.deepfeatures(h2oResultsDriver1, ae_model, layer=1)
summary(test_features_deep)

#Shutdown Server
h2o.shutdown(h2oServer, prompt = FALSE)

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
driverGridVal <- sample(drivers, 10)
RFGrid <- expand.grid(.numberOfNegativeDrivers = c(4, 10, 25),
                      .numberOfNegativeRows = c(50, 200, 650, 1500))
RFGrid <-  RFGrid[-10, ] #Remove row where sample is larger than the population when 'replace = FALSE' 

#Start h2o directly from R
h2oServer <- h2o.init(ip = "localhost", max_mem_size = "10g", nthreads = -1)

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
    #Validation Drivers Sampling
    ValDrivers <- ExtraDrivers[sample(1:nrow(ExtraDrivers), 7), ]
    #Extra drivers sampling
    ExtraDrivers <- ExtraDrivers[sample(seq(1, nrow(ExtraDrivers)), as.numeric(driversSplit[2])), ]  
    
    #R matrix conversion to h2o object and stored in the server
    #h2oResultPlusExtras <- as.h2o(h2oServer, cbind(c(rep(1, nrow(results)), rep(0, nrow(ExtraDrivers))), 
    #                                               signif(rbind(results, ExtraDrivers), digits = 4)))
    
    #R matrix conversion to h2o object and stored in the server
    h2oResultPlusExtras <- as.h2o(h2oServer, cbind(c(rep(1, nrow(results[1:100, ])), rep(0, nrow(ExtraDrivers))), 
                                                   signif(rbind(results[1:100, ], ExtraDrivers), digits = 4)))
    h2oValidation <- as.h2o(h2oServer, cbind(c(rep(1, nrow(results[101:nrow(results), ])), rep(0, nrow(ValDrivers))), 
                                             signif(rbind(results[101:nrow(results), ], ValDrivers), digits = 4)))
    
    print(h2o.ls(h2oServer))  
        
    #h2o.ai RF algorithm
    #Shuffle indexes
    #set.seed(1001001)
    #randIdxs <- sample(seq(1, nrow(results) + nrow(ExtraDrivers)), nrow(results) + nrow(ExtraDrivers))
    #Cross Validation + Modelling
    #driverRFModelCV <- h2o.randomForest(x = seq(2, ncol(h2oResultPlusExtras)), y = 1,
    #                                    data = h2oResultPlusExtras[randIdxs, ],
    #                                    nfolds = 5,
    #                                    classification = TRUE,
    #                                    ntree = c(50, 75, 100),
    #                                    depth = c(20, 50, 75), 
    #                                    verbose = FALSE)
    #aucError <- driverRFModelCV@model[[1]]@model$auc 
    
    randIdxs <- sample(seq(1, nrow(results[1:100, ]) + nrow(ExtraDrivers)), nrow(results[1:100, ]) + nrow(ExtraDrivers))
    #Simple Validation train 60% data - test 40%
    driverRFModelCV2 <- h2o.randomForest(x = seq(2, ncol(h2oResultPlusExtras)), y = 1,
                                        data = h2oResultPlusExtras[randIdxs, ],                            
                                        validation = h2oValidation,
                                        classification = TRUE,
                                        ntree = c(50, 75, 100),
                                        depth = c(20, 50, 75), 
                                        verbose = FALSE)
    
    aucError2 <- driverRFModelCV2@model[[1]]@model$auc
    print(aucError2)
    
    h2o.rm(object = h2oServer, keys = h2o.ls(h2oServer)[, 1])       
    #return(aucError)
    return(aucError2)
  })
  return(c(as.numeric(driversSplit), errorNthDriver))
})
#h2o Shutdown
h2o.shutdown(h2oServer, prompt = FALSE)

modelsGenerated <- t(modelsGenerated)

averageAUC <- apply(modelsGenerated[, c(-1, -2)], 1, mean)
bestDriverNumbers <- which.max(averageAUC)
driversParameters <- modelsGenerated[bestDriverNumbers, c(1, 2)]

#Modelling---------------------------
#Interrupting Remote Instance (in case this is running on an AWS spot instance)
SpotInstance <- TRUE
if (SpotInstance == TRUE){
  driversProcessed <- gsub(pattern = ".csv", replacement = "", x = list.files(outputDirectory))
  drivers <- c(driversProcessed[length(driversProcessed)],
               list.files(driversDirectory)[!(list.files(driversDirectory) %in% driversProcessed)])
}

driversPredictions <- lapply(drivers, function(driver){  
  #Parallel processing of each driver data
  resultsFull <- unlist(mclapply(seq(1, 200), transform2Percentiles, mc.cores = numCores, driverID = driver))
  resultsFull <- signif(matrix(resultsFull, nrow = 200, byrow = TRUE), digits = 4)
  print(paste0("Driver number ", driver, " processed"))    
  
  #Number of sampled trips
  sampleTrips <- 300
  #Define a number of repeated models with different driver's data
  nModels <- 2
  #Sample data from other drivers  
  numberOfDrivers <- 2
  initialDrivers <- sample(drivers[!drivers %in% driver], numberOfDrivers * nModels)
  negativeDrivers <- sapply(initialDrivers, function(driver){
    results <- unlist(mclapply(seq(1, 200), transform2Percentiles, mc.cores = numCores, driverID = driver))
    print(paste0("Driver number: ", driver, " processed"))
    return(results)
  })
  allNegativeDrivers <- signif(matrix(unlist(negativeDrivers), nrow = numberOfDrivers * nModels * 200, byrow = TRUE), digits = 4)
  
  #Start h2o from command line
  system(paste0("java -Xmx12G -jar ", h2o.jarLoc, " -port 54333 -name AXA &"))
  #Small pause
  Sys.sleep(3)
  #Connect R to h2o
  h2oServer <- h2o.init(ip = "localhost", port = 54333, nthreads = -1)
  
  #Anomaly Detection Scores with the whole group of trips from several drivers
  deepNNAnomalyScore <- AutoencodersAnomalyDetection(dataMatrix = rbind(resultsFull, allNegativeDrivers), nwServer = h2oServer)  
  
  #Insert Anomaly Scores into matrices
  resultsFull <- insertCol(resultsFull, ncol(resultsFull), v = deepNNAnomalyScore[1:200])
  allNegativeDrivers <- insertCol(allNegativeDrivers, ncol(allNegativeDrivers),
                                  v = deepNNAnomalyScore[201:length(deepNNAnomalyScore)])
  
  #Bad Data Removal in positive drivers
  idxBadData <- resultsFull[, ncol(resultsFull)] == 1
  results <- resultsFull[!resultsFull[, ncol(resultsFull)] == 1, -ncol(resultsFull)]
  
  #LOF Algorithm
  cleanDataLOF <- lofAnomalyDetection(fullDataMatrix = resultsFull, cleanDataMatrix = results)
  lofDriver <- cleanDataLOF[[1]]
  lofDriverRanking <- cleanDataLOF[[2]]
  fullDataLOF <- lofAnomalyDetection(fullDataMatrix = resultsFull, exclude_bad_data = FALSE)
  lofDriverFull <- fullDataLOF[[1]]
  lofDriverRankingFull <- fullDataLOF[[2]]
  
  #Create Indices for each extra drivers group
  groupsIntegers <- unlist(lapply(1:nModels, function(n){
    return(rep(n, 200 * numberOfDrivers))
  }))
  groupIndexes <- split(seq(1, nrow(allNegativeDrivers)), groupsIntegers)  
    
  #Select respective group
  ExtraDriversRandGroup <- allNegativeDrivers[groupIndexes[[1]], ]    
  #Bad Data Removal
  ExtraDriversRandGroup <- ExtraDriversRandGroup[!ExtraDriversRandGroup[, ncol(ExtraDriversRandGroup)] == 1, -ncol(ExtraDriversRandGroup)]
  #Extra drivers sampling
  ExtraDriversRandGroup <- ExtraDriversRandGroup[sample(seq(1, nrow(ExtraDriversRandGroup)), sampleTrips), ] 
  
  #Shuffle indexes
  #set.seed(1001001)
  randIdxs <- sample(seq(1, nrow(results) + nrow(ExtraDriversRandGroup)), nrow(results) + nrow(ExtraDriversRandGroup))  
  #R matrix conversion to h2o object and store in the server + bad data removal in negative cases
  h2oResultPlusExtras <- as.h2o(h2oServer, cbind(c(rep(1, nrow(results)), rep(0, nrow(ExtraDriversRandGroup))), 
                                                 signif(rbind(results, ExtraDriversRandGroup), digits = 4)))
  h2oResultsNthDriver <- as.h2o(h2oServer, cbind(rep(1, nrow(resultsFull)),
                                                 signif(resultsFull[, -ncol(resultsFull)], digits = 4)))
  print(h2o.ls(h2oServer))

  #h2o.ai RF algorithm
  #Cross Validation
  driverRFModelCV <- h2o.randomForest(x = seq(2, ncol(h2oResultPlusExtras)), y = 1,
                                      data = h2oResultPlusExtras[randIdxs, ],
                                      nfolds = 4,
                                      classification = TRUE,
                                      ntree = c(50, 75, 100),
                                      depth = c(20, 50, 75), 
                                      verbose = FALSE)
  
  #Log Info
  aucRF <- driverRFModelCV@model[[1]]@model$auc
  ntreeRF <- driverRFModelCV@model[[1]]@model$params$ntree
  depthRF <- driverRFModelCV@model[[1]]@model$params$depth

  print(h2o.ls(h2oServer))  
  h2o.rm(object = h2oServer, keys = h2o.ls(h2oServer)[, 1])  
  
  #Random Forest Multi-Model Loop
  randomForestPredictionsList <- lapply(1:nModels, function(modelNumber){
    #Select respective group
    ExtraDrivers <- allNegativeDrivers[groupIndexes[[modelNumber]], ]    
    #Bad Data Removal
    ExtraDrivers <- ExtraDrivers[!ExtraDrivers[, ncol(ExtraDrivers)] == 1, -ncol(ExtraDrivers)]
    #Extra drivers sampling
    ExtraDrivers <- ExtraDrivers[sample(seq(1, nrow(ExtraDrivers)), sampleTrips), ] 
    
    #Shuffle indexes
    #set.seed(1001001)
    randIdxs <- sample(seq(1, nrow(results) + nrow(ExtraDrivers)), nrow(results) + nrow(ExtraDrivers))
    
    #R matrix conversion to h2o object and stored in the server
    h2oResultPlusExtras <- as.h2o(h2oServer, cbind(c(rep(1, nrow(results)), rep(0, nrow(ExtraDrivers))), 
                                                   signif(rbind(results, ExtraDrivers), digits = 4)))
    h2oResultsNthDriver <- as.h2o(h2oServer, cbind(rep(1, nrow(resultsFull)),
                                                   signif(resultsFull[, -ncol(resultsFull)], digits = 4)))
    print(h2o.ls(h2oServer))
    
    #h2o.ai RF Modelling    
    driverRFModel <- h2o.randomForest(x = seq(2, ncol(h2oResultPlusExtras)), y = 1,
                                      data = h2oResultPlusExtras[randIdxs, ],
                                      classification = TRUE,
                                      type = "BigData",
                                      ntree = ntreeRF,
                                      depth = depthRF, 
                                      verbose = FALSE)    
    
    #probability Prediction of trips in Nth driver 
    predictionRFSingle <- signif(as.data.frame(h2o.predict(driverRFModel, newdata = h2oResultsNthDriver)[, 3]), digits = 4)
    #Bad Data rounding up to mean
    predictionRFSingle[idxBadData, 1] <- mean(predictionRFSingle[!idxBadData, 1], na.rm = TRUE)
    
    print(h2o.ls(h2oServer))  
    h2o.rm(object = h2oServer, keys = h2o.ls(h2oServer)[, 1])   
    print(paste0("Driver number ", driver, " processed with RFs")) 
    print(paste0("Model number ", modelNumber, " of ", nModels))
    
    return(predictionRFSingle)
  })   
  
  predictionRF <- apply(as.data.frame(randomForestPredictionsList), 1, mean)
    
  #GBM Multi-Model Loop
  #gbmPredictionsList <- lapply(1:nModels, function(modelNumber){
  #  #Select respective group
  #  ExtraDrivers <- allNegativeDrivers[groupIndexes[[modelNumber]], ]    
  #  #Bad Data Removal
  #  ExtraDrivers <- ExtraDrivers[!ExtraDrivers[, ncol(ExtraDrivers)] == 1, -ncol(ExtraDrivers)]
  #  #Extra drivers sampling
  #  ExtraDrivers <- ExtraDrivers[sample(seq(1, nrow(ExtraDrivers)), sampleTrips), ] 
  #  
  #  #Shuffle indexes
  #  #set.seed(1001001)
  #  randIdxs <- sample(seq(1, nrow(results) + nrow(ExtraDrivers)), nrow(results) + nrow(ExtraDrivers))
  #  
  #  #R matrix conversion to h2o object and stored in the server
  #  h2oResultPlusExtras <- as.h2o(h2oServer, cbind(c(rep(1, nrow(results)), rep(0, nrow(ExtraDrivers))), 
  #                                                 signif(rbind(results, ExtraDrivers), digits = 4)))
  #  h2oResultsNthDriver <- as.h2o(h2oServer, cbind(rep(1, nrow(resultsFull)),
  #                                                 signif(resultsFull[, -ncol(resultsFull)], digits = 4)))
  #  print(h2o.ls(h2oServer))
  #  
  #  #h2o.ai GBM Modelling    
  #  driverGBMModel <- h2o.gbm(x = seq(2, ncol(h2oResultPlusExtras)), y = 1,
  #                            data = h2oResultPlusExtras[randIdxs, ],                              
  #                            distribution = "bernoulli",
  #                            n.trees = 3000,
  #                            interaction.depth = 5,
  #                            shrinkage = 0.001)    
  #  
  #  #probability Prediction of trips in Nth driver 
  #  predictionGBMSingle <- signif(as.data.frame(h2o.predict(driverGBMModel, newdata = h2oResultsNthDriver)[, 3]), digits = 4)
  #  #Bad Data rounding up to one
  #  predictionGBMSingle[idxBadData, 1] <- mean(predictionGBMSingle[!idxBadData, 1], na.rm = TRUE) 
  #  
  #  print(h2o.ls(h2oServer))  
  #  h2o.rm(object = h2oServer, keys = h2o.ls(h2oServer)[, 1])   
  #  print(paste0("Driver number ", driver, " processed with RFs")) 
  #  print(paste0("Model number ", modelNumber, " of ", nModels))
  #  
  #  return(predictionGBMSingle)
  #})   
  #
  #predictionGBM <- apply(as.data.frame(gbmPredictionsList), 1, mean)
  
  #Shutdown h20 instance
  h2o.shutdown(h2oServer, prompt = FALSE)
  
  #Remove R Data
  rm(resultsFull, results, ExtraDriversRandGroup, allNegativeDrivers, negativeDrivers)
  
  #Logging hyperparameters
  write.csv(cbind(aucRF, ntreeRF, depthRF), 
            file = file.path(logsDirectory, paste0(driver, ".csv")), row.names = FALSE)
  
  print(paste0(which(list.files(driversDirectory) == driver), "/", length(list.files(driversDirectory))))  
  
  if (SpotInstance == TRUE){
    write.csv(cbind(lofDriver, lofDriverRanking,
                    predictionRF,
                    #predictionGBM, 
                    ), 
              file = file.path(outputDirectory, paste0(driver, ".csv")), row.names = FALSE)
    return(TRUE)
    
  }else{
    return(cbind(lofDriver, lofDriverRanking,
                 predictionRF,
                 #predictionGBM,
                 )
  }  
})

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
driversPredictions <- as.data.frame(do.call(rbind, driversPredictions))

#Trip Matching------------------
#Interrupting Remote Instance (in case this is running on an AWS spot instance)
SpotInstance <- TRUE
if (SpotInstance == TRUE){
  driversProcessed <- gsub(pattern = ".csv", replacement = "", x = list.files(tripMatchingDir))
  drivers <- c(driversProcessed[length(driversProcessed)],
               list.files(driversDirectory)[!(list.files(driversDirectory) %in% driversProcessed)])
}

#Run the trip matching algorithm
matchedTripsAllDrivers <- sapply(drivers, TripMatchingFun)

tripsMatchOutput <- list.files(tripMatchingDir)   
matchingPredictions <- lapply(tripsMatchOutput, function(driver){
  predictions <- fread(file.path(tripMatchingDir, driver), header = TRUE,
                       stringsAsFactors = FALSE)
  return(predictions)
})
#Concatenate prediction lists into a data.frame
matchedTripsAllDrivers <- as.data.frame(do.call(rbind, matchingPredictions))

#Prediction Boosting--------------
#Boost Best Predictions
bestBoostedPredictions <- driversPredictions[, 3]
bestPredictionsIdxs <- which(driversPredictions[matchedTripsAllDrivers, 3] > 0.8)
bestBoostedPredictions[matchedTripsAllDrivers[bestPredictionsIdxs]] <- 1 

#Boost All Predictions
driversPredictions[matchedTripsAllDrivers, ] <- 1

#Write .csv files-------------------------
submissionTemplate <- fread(file.path(otherDataDirectory, "sampleSubmission.csv"), header = TRUE,
                            stringsAsFactors = FALSE, colClasses = c("character", "numeric"))

#LOF Anomaly Score
submissionTemplate$prob <- signif(driversPredictions[, 1], digits = 8)
write.csv(submissionTemplate, file = "lofScoreIV.csv", row.names = FALSE)
system('zip lofScoreIV.zip lofScoreIV.csv')

#LOF Anomaly Rank
submissionTemplate$prob <- signif(driversPredictions[, 2], digits = 8)
write.csv(submissionTemplate, file = "lofRankIV.csv", row.names = FALSE)
system('zip lofRankIV.zip lofRankIV.csv')

#Probabilities h2o.ai RF
submissionTemplate$prob <- signif(driversPredictions[, 3], digits = 8)
write.csv(submissionTemplate, file = "RFProbIV.csv", row.names = FALSE)
system('zip RFProbIV.zip RFProbIV.csv')

#Probabilities h2o.ai GBM
#submissionTemplate$prob <- signif(driversPredictions[, 4], digits = 8)
#write.csv(submissionTemplate, file = "GBMProbIV.csv", row.names = FALSE)
#system('zip GBMProbIV.zip GBMProbIV.csv')

#Probabilities h2o.ai RF + boosted predictions with trajectory matching
submissionTemplate$prob <- signif(bestBoostedPredictions, digits = 8)
write.csv(submissionTemplate, file = "RFProbBoostIV.csv", row.names = FALSE)
system('zip RFProbBoostIV.zip RFProbBoostIV.csv')

