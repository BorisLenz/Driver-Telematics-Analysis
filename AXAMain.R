#AXA Driver Telematics Analysis
#Ver 0.61 #Added: Angles second revision and outlier detection second revision

#Init-----------------------------------------------
rm(list=ls(all=TRUE))

#Libraries, directories, options and extra functions----------------------
require("data.table")
require("parallel")
require("h2o")
require("DMwR")
require("ggplot2")
require("plotly")

#Set Working Directory
workingDirectory <- "/home/wacax/Wacax/Kaggle/AXA Driver Telematics Analysis/"
setwd(workingDirectory)
driversDirectory <- "/home/wacax/Wacax/Kaggle/AXA Driver Telematics Analysis/Data/drivers"
otherDataDirectory <- "/home/wacax/Wacax/Kaggle/AXA Driver Telematics Analysis/Data/"
vw77Dir = "/home/wacax/vowpal_wabbit-7.7/vowpalwabbit/"
#List all possible drivers identities
drivers <- list.files(driversDirectory)
#Detect available cores
numCores <- detectCores()

#Data Mining (Functions)------------------------
#Outlier Removal and transformation to quantiles
quantSigma <- function(vector, sigma = 5, returnVector = TRUE){
  #n-sigma removal
  vectorWithoutOutliers <- vector[!vector > mean(vector) + sd(vector) * sigma]
  if (returnVector == TRUE){
    return(list(quantile(vectorWithoutOutliers, seq(0.05, 1, by=0.05)), vectorWithoutOutliers))
  }else{
    return(quantile(vectorWithoutOutliers, seq(0.05, 1, by=0.05)))    
  }
}
#Transform data to speed distributions
speedDistribution <- function(trip){
  speed <- 3.6 * sqrt(diff(trip$x)^2 + diff(trip$y)^2)
  #n-sigma removal
  speedList <- quantSigma(speed)
  return(speedList)
}
#Angles Mining with Angle Rotation
angleVector <- function(dir, sigma = 4){
  dir <- as.data.table(cbind(diff(dir$x), diff(dir$y)))
  setnames(dir, old = c("V1", "V2"), new = c("x", "y"))
  #Returns the angles between vectors
  #dir is a 2D-array of shape (N,M) representing N vectors in M-dimensional space.
  #It returns a vector which is a 1D-array of values of shape (N-1,), with each value between 0 and pi.
  dir2 <- dir[2:nrow(dir)]
  dir1 <- dir[1:(nrow(dir) - 1)]
  angles <- acos(rowSums(dir1 * dir2) / ((sqrt(rowSums(dir1^ 2) * rowSums(dir2 ^2))) + 1e-06)) * (180/pi) #1e-06 is included to avoid a division by zero
  #NAs removal
  angles <- na.omit(angles)
  #set orthogonal values to zero
  angles[angles >= 90] <- 0  
  #Outliers removal
  anglesWoOutliers <- quantSigma(angles, sigma = 4)
  return(anglesWoOutliers)
}
#Define function to be passed as parallel
transform2Percentiles <- function(file, driverID){
  trip <- fread(file.path(driversDirectory, driverID, paste0(file, ".csv")))
  #Velocities
  speedData <- speedDistribution(trip)
  speedDist <- speedData[[1]]
  accelerations <- diff(speedData[[2]])
  #Accelerations
  accelerationDist <- quantSigma(accelerations, sigma = 4,  returnVector = FALSE)
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
  
  return(c(speedDist, accelerationDist, distanceTrip, turningAnglesDist))
  #return(c(speedDist, accelerationDist, distanceTrip, turningAnglesDist, anglesPlusSpeedist))
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
#If there is need to start h2o from command line:
#system(paste0("java -Xmx55G -jar ", h2o.jarLoc, " -port 54333 -name CTR -data_max_factor_levels 100000000 &"))
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
}, mc.cores = numCores, driverID = driverViz)
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

#Unsupervised Learning and Hyperparameter Tuning--------------
#Neural Network pre-training
#Begin with a randomly selected driver to start the unsupervised learning
numberOfDrivers <- 25
initialDrivers <- sample(drivers, numberOfDrivers)
driversProcessed <- sapply(initialDrivers, function(driver){
  results <- unlist(mclapply(seq(1, 200), transform2Percentiles, mc.cores = numCores, driverID = driver))
  print(paste0("Driver number: ", driver, " processed"))
  return(results)
})
driversProcessed <- scale(matrix(unlist(driversProcessed), nrow = numberOfDrivers * 200, byrow = TRUE))

#Init h2o Server
#If there is need to start h2o from command line:
#system(paste0("java -Xmx55G -jar ", h2o.jarLoc, " -port 54333 -name CTR -data_max_factor_levels 100000000 &"))
#Start h2o directly from R
h2oServer <- h2o.init(ip = "localhost", max_mem_size = "5g", nthreads = -1)
h2oResult <- as.h2o(h2oServer, cbind(rep(c(0, 1, 1, 0), 50 * numberOfDrivers), driversProcessed))
print(h2o.ls(h2oServer))
rm(driversProcessed)
activations <- c("RectifierWithDropout", "TanhWithDropout", "Rectifier", "Tanh")
cvNNModel <- h2o.deeplearning(x = seq(2, ncol(h2oResult)), y = 1,
                              data = h2oResult,
                              autoencoder = TRUE, fast_mode = TRUE,
                              activation = activations,
                              input_dropout_ratio = c(0, 0.1),
                              l1 = c(0, 1e-5),
                              l2 = c(0, 1e-5),
                              rho = c(0.95, 0.99),
                              epsilon = c(1e-12, 1e-10, 1e-08),
                              hidden = c(50, 40, 50), epochs = 50)
checkpointModel <- cvNNModel@model[[1]] #Best NN cv model
deepNetPath <- h2o.saveModel(object = checkpointModel, dir = file.path(workingDirectory), force = TRUE)
h2o.shutdown(h2oServer, prompt = FALSE)

#Modelling---------------------------
#Init h2o Server
#If there is need to start h2o from command line:
#system(paste0("java -Xmx55G -jar ", h2o.jarLoc, " -port 54333 -name CTR -data_max_factor_levels 100000000 &"))
#Start h2o directly from R
h2oServer <- h2o.init(ip = "localhost", max_mem_size = "5g", nthreads = -1)
#Load pretrained model
checkpointModel <- h2o.loadModel(h2oServer, deepNetPath)
print(h2o.ls(h2oServer))
checkpointModelKey <- h2o.ls(h2oServer)[, 1]

driversProb <- sapply(drivers, function(driver){
  #Parallel processing of each driver data
  results <- unlist(mclapply(seq(1, 200), transform2Percentiles, mc.cores = numCores, driverID = driver))
  results <- scale(matrix(results, nrow = 200, byrow = TRUE))
  print(paste0("Driver number ", driver, " processed"))
  
  #Sample data from other drivers  
  numberOfDrivers <- 10
  initialDrivers <- sample(drivers[!drivers %in% driver], numberOfDrivers)
  ExtraDrivers <- sapply(initialDrivers, function(driver){
    results <- unlist(mclapply(seq(1, 200), transform2Percentiles, mc.cores = numCores, driverID = driver))
    print(paste0("Driver number: ", driver, " processed"))
    return(results)
  })
  ExtraDrivers <- scale(matrix(unlist(ExtraDrivers), nrow = numberOfDrivers * 200, byrow = TRUE))  
  
  #LOF Algorithm
  lofDriver <- lofactor(results, k=5)
  lofDriverRanking <- rank(-lofDriver)
  
  #h20.ai deep learning
  #R matrix conversion to h2o object and stored in the server
  h2oResult <- as.h2o(h2oServer, cbind(c(rep(1, nrow(results)), rep(0, nrow(ExtraDrivers))), 
                                       rbind(results, ExtraDrivers)))
  h2oOnlyDriverTrips <- h2oResult[1:200, ]
  h2oResult <- h2o.splitFrame(h2oResult, ratios = 1, shuffle = TRUE)[[1]]
  
  print(h2o.ls(h2oServer))
  driverDeepNNModel <- h2o.deeplearning(x = seq(2, ncol(h2oResult)), y = 1,
                                        activation = optimalActivation,
                                        data = h2oResult, 
                                        checkpoint = checkpointModel,
                                        input_dropout_ratio = c(0, 1),
                                        l1 = c(0, 1e-5),
                                        l2 = c(0, 1e-5),
                                        rho = c(0.95, 0.99),
                                        epsilon = c(1e-12, 1e-10, 1e-08),
                                        hidden = c(50, 40, 50), epochs = 500)
  
  #anomalousTrips <- as.data.frame(h2o.anomaly(h2oResult, driverDeepNNModel))
  
  #MSE error transformation into pseudo-probabilities / chi-squared probability calculation
  #anomalousTrips[, 1] <- pchisq(anomalousTrips[, 1], df = 1)
  
  #probability Prediction of trips in driver #N
  predictionNN <- signif(as.data.frame(h2o.predict(DeepNNModel, newdata = h2oOnlyDriverTrips)), digits = 4)
  predictionNNRank <- rank(predictionNN[, 1])
    
  print(h2o.ls(h2oServer))
  h2oObjects2Remove <- which(!h2o.ls(h2oServer)[, 1] %in% checkpointModelKey)
  h2o.rm(object = h2oServer, keys = h2o.ls(h2oServer)[h2oObjects2Remove, 1])  
  print(paste0(which(drivers == driver), "/", length(drivers)))
  
  return(cbind(predictionNN[, 1], predictionNNRank, lofDriver, lofDriverRanking))
})
h2o.shutdown(h2oServer, prompt = FALSE)

deepNNProb <- 
deepNNRank <- 
lofScore <- 
lofRank <- 

#Write .csv files-------------------------
submissionTemplate <- fread(file.path(otherDataDirectory, "sampleSubmission.csv"), header = TRUE,
                            stringsAsFactors = FALSE, colClasses = c("character", "numeric"))
#Probabilities h2o.ai Deep NN
submissionTemplate$prob <- signif(deepNNProb, digits = 5)
write.csv(submissionTemplate, file = "deepNNProbI.csv", row.names = FALSE)
system('zip deepNNProbI.zip deepNNProbI.csv')

#Probabilities h2o.ai Deep NN
submissionTemplate$prob <- signif(deepNNRank, digits = 5)
write.csv(submissionTemplate, file = "deepNNRankI.csv", row.names = FALSE)
system('zip deepNNRankI.zip SdeepNNRankI.csv')

#Probabilities h2o.ai Deep NN
submissionTemplate$prob <- signif(lofScore, digits = 5)
write.csv(submissionTemplate, file = "lofScoreI.csv", row.names = FALSE)
system('zip lofScoreI.zip lofScoreI.csv')

#Probabilities h2o.ai Deep NN
submissionTemplate$prob <- signif(lofRank, digits = 5)
write.csv(submissionTemplate, file = "lofRankI.csv", row.names = FALSE)
system('zip lofRankI.zip lofRankI.csv')
