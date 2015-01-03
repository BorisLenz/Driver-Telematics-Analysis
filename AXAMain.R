#AXA Driver Telematics Analysis
#Ver 0.3  #Added: Turns, PCA, cross validation and pre-training

#Init-----------------------------------------------
rm(list=ls(all=TRUE))

#Libraries, Options and extra functions----------------------
require("data.table")
require("parallel")
require("h2o")
require("ggplot2")

#Set Working Directory------------------------------
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
#Transform data to distributions
speedDistribution <- function(trip, sigma = 5){
  speed <-  3.6 * sqrt(diff(trip$x)^2 + diff(trip$y)^2) 
  #Six sigma removal
  speedWoOutliers <- speed[!speed > sd(speed) * sigma]  
  return(list(quantile(speedWoOutliers, seq(0.05, 1, by=0.05)), speed))
}

#Define function to be passed as parallel
transform2Percentiles <- function(file, driverID){
  trip <- fread(file.path(driversDirectory, driverID, paste0(file, ".csv")))
  speedData <- speedDistribution(trip)
  speedDist <- speedData[[1]]
  accelerationDist <- quantile(diff(speedData[[2]]), seq(0.05, 1, by=0.05))
  distanceTrip <- sum(sqrt((diff(trip$x)^2) + (diff(trip$y)^2)))
  turningAngles <- diff(atan2(diff(trip$y), diff(trip$x)) * (180/pi))
  turningAngles <- quantile(turningAngles, seq(0.05, 1, by=0.05))
  
  return(c(speedDist, accelerationDist, distanceTrip, turningAngles))
}

#EDA----------------------------------------
## EDA Pt. 1 Determine the minimal PCAs / number of neurons in the middle layer
#Begin with a randomly selected driver to start the PCA calculation
numberOfDrivers <- 300
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

##EDA Pt. 3 Visualization of Speeds with and without outlier replacement
driverViz <- sample(drivers, 1)
fileViz <- sample(1:200, 1)

tripCoordinates <- fread(file.path(driversDirectory, driverViz, paste0(fileViz, ".csv")))
print(paste0("Driver number: ", driverViz, " trip number ", fileViz, " processed"))

speed2Plot <-  3.6 * sqrt(diff(tripCoordinates$x)^2 + diff(tripCoordinates$y)^2)
ggplot(as.data.frame(speed2Plot), aes(x = speed2Plot)) + geom_density()
#Remove data above 5 sigmas (5 standard deviations)
speed2Plot2 <- speed2Plot[!speed2Plot > sd(speed2Plot) * 5]
ggplot(as.data.frame(speed2Plot2), aes(x = speed2Plot2)) + geom_density()

#Unsupervised Learning and Hyperparameter Tuning--------------
#Neural Network pre-training
#Begin with a randomly selected driver to start the unsupervised learning
numberOfDrivers <- 75
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

h2oResult <- as.h2o(h2oServer, cbind(rep(c(0, 1, 1, 0), 50), driversProcessed)) 
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
                              hidden = c(55, 35, 55), epochs = 250)

checkpointModel <- cvNNModel@model[[1]] #Best NN cv model 
deepNetPath <- h2o.saveModel(object = checkpointModel, dir = file.path(workingDirectory), force = TRUE)
optimalActivation <- cvNNModel@model[[1]]@model$params$activation
optimalIDR <- cvNNModel@model[[1]]@model$params$input_dropout_ratio
optimall1 <- cvNNModel@model[[1]]@model$params$l1
optimall2 <- cvNNModel@model[[1]]@model$params$l2
optimalRho <- cvNNModel@model[[1]]@model$params$rho
optimalEpsilon <- cvNNModel@model[[1]]@model$params$epsilon

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
  #KSVMmodel <- ksvm(results, type = "one-svc", kernel="rbfdot", kpar="automatic", prob.model = TRUE)
  #KSVMPrediction <- predict(KSVMmodel, results, type = "response")
    
  #Plotting
  #par(mfrow=c(1, 2))
  #plot(density(lofactor(results, k=5)))
  #biplot(prcomp(results), cex=.8) 
  
  #R matrix conversion to h2o object and stored in the server
  h2oResult <- as.h2o(h2oServer, cbind(rep(c(0, 1, 1, 0), 50), results))
  print(h2o.ls(h2oServer))
  driverDeepNNModel <- h2o.deeplearning(x = seq(2, ncol(h2oResult)), y = 1, 
                                        activation = optimalActivation,
                                        data = h2oResult, autoencoder = TRUE,
                                        checkpoint = checkpointModel,
                                        input_dropout_ratio = optimalIDR,
                                        l1 = optimall1,
                                        l2 = optimall2,
                                        rho = optimalRho,
                                        epsilon = optimalEpsilon,
                                        hidden = c(55, 35, 55), epochs = 500)
  
  anomalousTrips <- as.data.frame(h2o.anomaly(h2oResult, driverDeepNNModel))
  
  #MSE error transformation into pseudo-probabilities / chi-squared probability calculation
  #anomalousTrips[, 1] <- pchisq(anomalousTrips[, 1], df = 1)  
  
  print(h2o.ls(h2oServer))
  h2oObjects2Remove <- which(!h2o.ls(h2oServer)[, 1] %in% checkpointModelKey)
  h2o.rm(object = h2oServer, keys = h2o.ls(h2oServer)[h2oObjects2Remove, 1])  
  print(paste0("Driver number ", driver, " processed"))
  return(anomalousTrips[, 1])
})

h2o.shutdown(h2oServer, prompt = FALSE)  

#MSE error transformation into positive pseudo-probabilities-------------------
#chi-squared probability calculation
driversProb <- pchisq(driversProb, df = 1)
driversProb <- 1 - as.vector(driversProb)

#Write .csv------------------------- 
submissionTemplate <- fread(file.path(otherDataDirectory, "sampleSubmission.csv"), header = TRUE, 
                            stringsAsFactors = FALSE, colClasses = c("character", "numeric"))

submissionTemplate$prob <- signif(driversProb, digits = 5)
write.csv(submissionTemplate, file = "SpeedNNMSEPredictionV.csv", row.names = FALSE)
system('zip SpeedNNMSEPredictionV.zip SpeedNNMSEPredictionV.csv')
