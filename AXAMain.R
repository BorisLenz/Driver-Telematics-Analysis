#AXA Driver Telematics Analysis
#Ver 0.7.10 #code altered for possible use in an AWS spot instance

#Init-----------------------------------------------
rm(list=ls(all=TRUE))

#Libraries, directories, options and extra functions----------------------
require("data.table")
require("parallel")
require("h2o")
require("DMwR")
require("ggplot2")
#require("plotly")

#Set Working Directory
workingDirectory <- "/home/wacax/Wacax/Kaggle/AXA-Driver-Telematics-Analysis/"
setwd(workingDirectory)
driversDirectory <- "/home/wacax/Wacax/Kaggle/AXA-Driver-Telematics-Analysis/Data/drivers"
otherDataDirectory <- "/home/wacax/Wacax/Kaggle/AXA-Driver-Telematics-Analysis/Data/"
outputDirectory <- "/home/wacax/Wacax/Kaggle/AXA-Driver-Telematics-Analysis/Data/Output"
vw77Dir = "/home/wacax/vowpal_wabbit-7.7/vowpalwabbit/"
#h2o location
h2o.jarLoc <- "/home/wacax/R/x86_64-pc-linux-gnu-library/3.1/h2o/java/h2o.jar"

#List all possible drivers identities
drivers <- list.files(driversDirectory)
#Detect available cores
numCores <- detectCores()

#Extra Functions
source(paste0(workingDirectory, "pretrainh2oNN.R"))

#Data Mining (Functions)------------------------
#Outlier Removal and transformation to quantiles
quantSigma <- function(vector, sigma = 4, returnVector = TRUE){
  #n-sigma removal
  vectorWithoutOutliers <- vector[!vector > mean(vector) + sd(vector) * sigma]
  if (returnVector == TRUE){
    return(list(quantile(vectorWithoutOutliers, seq(0.05, 1, by=0.05)), vectorWithoutOutliers))
  }else{
    return(quantile(vectorWithoutOutliers, seq(0.05, 1, by=0.05)))    
  }
}
#Transform data to speed distributions
speedDistribution <- function(coordinates){
  speed <- 3.6 * sqrt(diff(coordinates$x)^2 + diff(coordinates$y)^2)
  #n-sigma removal
  speedList <- quantSigma(speed)
  return(speedList)
}
#Angles Mining with Angle Rotation
angleVector <- function(coordinatesAng){
  coorDT <- as.data.table(cbind(diff(coordinatesAng$x), diff(coordinatesAng$y)))
  #Returns the angles between vectors
  #coordinatesAng is a 2D-array of shape (N,M) representing N vectors in M-dimensional space.
  #It returns a vector which is a 1D-array of values of shape (N-1,), with each value between 0 and pi.
  dir2 <- coorDT[2:nrow(coorDT)]
  dir1 <- coorDT[1:(nrow(coorDT) - 1)]
  angles <- acos(rowSums(dir1 * dir2) / ((sqrt(rowSums(dir1^ 2) * rowSums(dir2 ^2))) + 1e-06)) * (180/pi) #1e-06 is included to avoid a division by zero
  #NAs removal
  angles <- na.omit(angles)
  #set orthogonal values to zero
  angles[angles >= 89] <- 0  
  #Outliers removal
  anglesWoOutliers <- quantSigma(angles)
  return(anglesWoOutliers)
}
#Define function to be passed as parallel
transform2Percentiles <- function(file, driverID){
  trip <- fread(file.path(driversDirectory, driverID, paste0(file, ".csv")))
  #Velocities
  speedData <- speedDistribution(trip)
  speedDist <- speedData[[1]]
  #Accelerations
  accelerations <- diff(speedData[[2]])
  accelerationDist <- quantSigma(accelerations, returnVector = FALSE)
  #Distances
  distances <- sqrt((diff(trip$x)^2) + (diff(trip$y)^2))
  distances[distances > mean(distances) + sd(distances) * 5] <- NA
  distances[is.na(distances)] <- mean(distances, na.rm = TRUE)
  distanceTrip <- sum(distances)
  #Angles
  anglesData <- angleVector(trip)
  turningAnglesDist <- anglesData[[1]]
  #anglesPlusSpeedist <- anglesData[[2]] * speedData[[2]][-1] #this doesn't work on the 1004th driver trip number 11
  #anglesPlusSpeedist <- anglesPlusSpeedist[!anglesPlusSpeedist > sd(anglesPlusSpeedist) * 5]
  #anglesPlusSpeedist <- quantile(anglesPlusSpeedist, seq(0.05, 1, by=0.05))
  
  return(c(speedDist, accelerationDist, turningAnglesDist, 
           sd(accelerations[accelerations > 0]), sd(accelerations[accelerations < 0]), distanceTrip, nrow(trip)))
}

#EDA----------------------------------------
## EDA Pt. 1 Determine the minimal PCAs / number of neurons in the middle layer
#Begin with a randomly selected driver to start the PCA calculation
numberOfDrivers <- 100
initialDrivers <- sample(drivers, numberOfDrivers)
driversProcessed <- sapply(initialDrivers, function(driver){
  results <- unlist(mclapply(seq(1, 200), transform2Percentiles, mc.cores = numCores, driverID = driver))
  print(paste0("Driver number: ", driver, " processed"))
  return(results)
})
driversProcessed <- scale(matrix(unlist(driversProcessed), nrow = numberOfDrivers * 200, byrow = TRUE))

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

##EDA Pt. 2 Visualization of Trajectories
driverViz <- sample(drivers, 1)
results <- mclapply(seq(1, 200), function(file, driverID){
  tripCoordinates <- fread(file.path(driversDirectory, driverID, paste0(file, ".csv")))
  return(cbind(tripCoordinates, rep(file, nrow(tripCoordinates))))
}, mc.cores = numCores, driverID = driverViz)https://github.com/h2oai/h2o/blob/master/R/examples/Kaggle/CTR.R
print(paste0("Driver number: ", driverViz, " processed"))
dataTableReady2Plot <- do.call(rbind, results)
qplot(x, y, data = dataTableReady2Plot, colour = V2, geom = "point")
#py <- plotly()
#py$ggplotly()

##EDA Pt. 3 Visualization of Speeds with and without outlier replacement
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

##EDA Pt. 4 Visualization of Turning Angles with and without outlier replacement
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

##EDA Pt. 5 LOF Algorithm visualization
driverViz <- sample(drivers, 1)
results <- unlist(mclapply(seq(1, 200), transform2Percentiles, mc.cores = numCores, driverID = driverViz))
results <- scale(matrix(results, nrow = 200, byrow = TRUE))
#Plot
par(mfrow=c(1, 2))
plot(density(lofactor(results, k=5)))
biplot(prcomp(results), cex=.8)

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
system(paste0("java -Xmx10G -jar ", h2o.jarLoc, " -port 54333 -name AXA &"))
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
  results <- unlist(mclapply(seq(1, 200), transform2Percentiles, mc.cores = numCores, driverID = driver))
  results <- signif(matrix(results, nrow = 200, byrow = TRUE), digits = 3)
  print(paste0("Driver number ", driver, " processed"))
  
  #Sample data from other drivers  
  numberOfDrivers <- 2
  initialDrivers <- sample(drivers[!drivers %in% driver], numberOfDrivers)
  ExtraDrivers <- sapply(initialDrivers, function(driver){
    results <- unlist(mclapply(seq(1, 200), transform2Percentiles, mc.cores = numCores, driverID = driver))
    print(paste0("Driver number: ", driver, " processed"))
    return(results)
  })
  ExtraDrivers <- signif(matrix(unlist(ExtraDrivers), nrow = numberOfDrivers * 200, byrow = TRUE), digits = 3)
  ExtraDrivers <- ExtraDrivers[sample(seq(1, nrow(ExtraDrivers)), 50), ]
  
  #LOF Algorithm
  lofDriver <- signif(-(lofactor(results, k=10)) + 10, digits = 4)
  lofDriverRanking <- rank(lofDriver)
  print(paste0("Driver number ", driver, " processed with the LOF Algorithm"))  
  
  #Shuffle indexes
  #set.seed(1001001)
  randIdxs <- sample(seq(1, nrow(results) + nrow(ExtraDrivers)), nrow(results) + nrow(ExtraDrivers))
  
  #R matrix conversion to h2o object and stored in the server
  h2oResultPlusExtras <- as.h2o(h2oServer, cbind(c(rep(1, nrow(results)), rep(0, nrow(ExtraDrivers))), 
                                                 signif(scale(rbind(results, ExtraDrivers)), digits = 4)))
  h2oResultsNthDriver <- h2oResultPlusExtras[1:nrow(results), ]
  print(h2o.ls(h2oServer))
  
  #h2o.ai GBM algorithm
  #Cross Validation + Modelling
  driverRFModelCV <- h2o.randomForest(x = seq(2, ncol(h2oResultPlusExtras)), y = 1,
                                      data = h2oResultPlusExtras[randIdxs, ],
                                      classification = TRUE,
                                      ntree = c(50, 75, 100),
                                      depth = c(20, 40, 60), 
                                      verbose = TRUE)
  
  driverRFModel <- driverRFModelCV@model[[1]]
  print(driverRFModel)
  
  #probability Prediction of trips in Nth driver 
  predictionRF <- signif(as.data.frame(h2o.predict(driverRFModel, newdata = h2oResultsNthDriver)[, 3]), digits = 4)
  predictionRFRank <- rank(predictionRF[, 1])
  print(h2o.ls(h2oServer))
  
  h2o.rm(object = h2oServer, keys = h2o.ls(h2oServer)[, 1])   
  print(paste0("Driver number ", driver, " processed with RFs")) 
  
  #h2o.ai GBM algorithm
  #R matrix conversion to h2o object and stored in the server
  h2oResultPlusExtras <- as.h2o(h2oServer, cbind(c(rep(1, nrow(results)), rep(0, nrow(ExtraDrivers))), 
                                                 signif(scale(rbind(results, ExtraDrivers)), digits = 4)))
  h2oResultsNthDriver <- h2oResultPlusExtras[1:nrow(results), ]
  print(h2o.ls(h2oServer))
  
  #Cross Validation
  driverGBMModelCV <-  h2o.gbm(x = seq(2, ncol(h2oResultPlusExtras)), y = 1,
                               data = h2oResultPlusExtras[randIdxs, ],
                               distribution = "bernoulli",
                               interaction.depth = c(2, 4, 7),
                               shrinkage = c(0.001, 0.003), 
                               n.trees = 100)
  
  #Modelling
  driverGBMModel <-  h2o.gbm(x = seq(2, ncol(h2oResultPlusExtras)), y = 1,
                             data = h2oResultPlusExtras[randIdxs, ],
                             distribution = "bernoulli",
                             interaction.depth = driverGBMModelCV@model[[1]]@model$params$interaction.depth,
                             shrinkage = driverGBMModelCV@model[[1]]@model$params$shrinkage, 
                             n.trees = 2000)  
  
  print(driverGBMModel)
  
  #probability Prediction of trips in Nth driver 
  predictionGBM <- signif(as.data.frame(h2o.predict(driverGBMModel, newdata = h2oResultsNthDriver)[, 3]), digits = 4)
  predictionGBMRank <- rank(predictionGBM[, 1])
  print(h2o.ls(h2oServer))
  
  h2o.rm(object = h2oServer, keys = h2o.ls(h2oServer)[, 1])     
  print(paste0("Driver number ", driver, " processed with GBMs"))  
  
  #h20.ai deep learning algorithm 
  #R matrix conversion to h2o object and stored in the server
  h2oResultPlusExtras <- as.h2o(h2oServer, cbind(c(rep(1, nrow(results)), rep(0, nrow(ExtraDrivers))), 
                                                 signif(scale(rbind(results, ExtraDrivers)), digits = 4)))
  h2oResultsNthDriver <- h2oResultPlusExtras[1:nrow(results), ]
  
  #Cross Validation + model
  driverDeepNNModelCV <- h2o.deeplearning(x = seq(2, ncol(h2oResultPlusExtras)), y = 1,
                                          data = h2oResultPlusExtras[randIdxs, ],   
                                          classification = TRUE,
                                          #checkpoint = ifelse("deepNetPath" %in% ls(), checkpointModel, ""),
                                          activation = "Tanh",
                                          input_dropout_ratio = c(0, 0.1),
                                          l1 = c(0, 1e-5),
                                          l2 = c(0, 1e-5),
                                          rho = c(0.95, 0.99),
                                          epsilon = c(1e-12, 1e-10),
                                          hidden = c(60, 60, 40, 60, 60), 
                                          epochs = 60)  
  
  driverDeepNNModel <- driverDeepNNModelCV@model[[1]]
  driverDeepNNModel2 <- driverDeepNNModelCV@model[[2]]
  driverDeepNNModel3 <- driverDeepNNModelCV@model[[3]]
  
  print(driverDeepNNModel)
  print(driverDeepNNModel2)
  print(driverDeepNNModel3)
  
  #probability Predictions on all trips in Nth driver 
  predictionNN <- as.data.frame(h2o.predict(driverDeepNNModel, newdata = h2oResultsNthDriver)[, 3])
  predictionNNRank <- rank(predictionNN[, 1])
  
  predictionNN2 <- as.data.frame(h2o.predict(driverDeepNNModel2, newdata = h2oResultsNthDriver)[, 3])
  predictionNNRank2 <- rank(predictionNN2[, 1])
  
  predictionNN3 <- as.data.frame(h2o.predict(driverDeepNNModel3, newdata = h2oResultsNthDriver)[, 3])
  predictionNNRank3 <- rank(predictionNN3[, 1])
  print(h2o.ls(h2oServer))
  
  h2o.rm(object = h2oServer, keys = h2o.ls(h2oServer)[, 1])   
  print(paste0("Driver number ", driver, " processed with Deep NNs"))  
  
  #h2oObjects2Remove <- which(!h2o.ls(h2oServer)[, 1] %in% checkpointModelKey)
  #h2o.rm(object = h2oServer, keys = h2o.ls(h2oServer)[h2oObjects2Remove, 1]) 
  print(paste0(which(drivers == driver), "/", length(list.files(driversDirectory))))
  
  if (SpotInstance == TRUE){
    write.csv(cbind(lofDriver, lofDriverRanking,
                    predictionRF[, 1], predictionRFRank, 
                    predictionGBM[, 1], predictionGBMRank, 
                    predictionNN[, 1], predictionNNRank, 
                    predictionNN2[, 1], predictionNNRank2, 
                    predictionNN3[, 1], predictionNNRank3), 
              file = file.path(outputDirectory, paste0(driver, ".csv")), row.names = FALSE)
    return(TRUE)
    
  }else{
    return(cbind(lofDriver, lofDriverRanking,
                 predictionRF[, 1], predictionRFRank, 
                 predictionGBM[, 1], predictionGBMRank, 
                 predictionNN[, 1], predictionNNRank,
                 predictionNN2[, 1], predictionNNRank2, 
                 predictionNN3[, 1], predictionNNRank3))
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
    return(predictions)
  })
}

#Concatenate lists into a data.frame
driversPredictions <- do.call(rbind, driversPredictions)

lofScore <- driversPredictions[, 1]
lofRank <- driversPredictions[, 2]
RFProb <- driversPredictions[, 3]
RFRank <- driversPredictions[, 4]
GBMProb <- driversPredictions[, 5]
GBMRank <- driversPredictions[, 6]
deepNNProb <- driversPredictions[, 7]
deepNNRank <- driversPredictions[, 8]

#Write .csv files-------------------------
submissionTemplate <- fread(file.path(otherDataDirectory, "sampleSubmission.csv"), header = TRUE,
                            stringsAsFactors = FALSE, colClasses = c("character", "numeric"))

#Probabilities h2o.ai RF
submissionTemplate$prob <- signif(RFProb, digits = 4)
write.csv(submissionTemplate, file = "RFProbI.csv", row.names = FALSE)
system('zip RFProbI.zip RFProbI.csv')

#Probabilities h2o.ai RF
submissionTemplate$prob <- signif(RFRank, digits = 4)
write.csv(submissionTemplate, file = "RFRankI.csv", row.names = FALSE)
system('zip RFRankI.zip RFRankI.csv')

#Probabilities h2o.ai GBM
submissionTemplate$prob <- signif(GBMProb, digits = 4)
write.csv(submissionTemplate, file = "GBMProbI.csv", row.names = FALSE)
system('zip GBMProbI.zip GBMProbI.csv')

#Probabilities h2o.ai GBM
submissionTemplate$prob <- signif(GBMRank, digits = 4)
write.csv(submissionTemplate, file = "GBMRankI.csv", row.names = FALSE)
system('zip GBMRankI.zip GBMRankI.csv')

#Probabilities h2o.ai Deep NN
submissionTemplate$prob <- signif(deepNNProb, digits = 4)
write.csv(submissionTemplate, file = "deepNNProbI.csv", row.names = FALSE)
system('zip deepNNProbI.zip deepNNProbI.csv')

#Probabilities h2o.ai Deep NN
submissionTemplate$prob <- signif(deepNNRank, digits = 4)
write.csv(submissionTemplate, file = "deepNNRankI.csv", row.names = FALSE)
system('zip deepNNRankI.zip deepNNRankI.csv')

#Probabilities h2o.ai Deep NN
submissionTemplate$prob <- signif(lofScore, digits = 4)
write.csv(submissionTemplate, file = "lofScoreI.csv", row.names = FALSE)
system('zip lofScoreI.zip lofScoreI.csv')

#Probabilities h2o.ai Deep NN
submissionTemplate$prob <- signif(lofRank, digits = 4)
write.csv(submissionTemplate, file = "lofRankI.csv", row.names = FALSE)
system('zip lofRankI.zip lofRankI.csv')
