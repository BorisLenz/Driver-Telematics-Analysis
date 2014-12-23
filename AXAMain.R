#AXA Driver Telematics Analysis
#Ver 0.2  #Added extra features

#Init-----------------------------------------------
rm(list=ls(all=TRUE))

#Libraries, Options and extra functions----------------------
require("data.table")
require("parallel")
require("h2o")
require("caret")

#Set Working Directory------------------------------
workingDirectory <- "/home/wacax/Wacax/Kaggle/AXA Driver Telematics Analysis/"
setwd(workingDirectory)
driversDirectory <- "/home/wacax/Wacax/Kaggle/AXA Driver Telematics Analysis/Data/drivers"
otherDataDirectory <- "/home/wacax/Wacax/Kaggle/AXA Driver Telematics Analysis/Data/"

vw77Dir = "/home/wacax/vowpal_wabbit-7.7/vowpalwabbit/"

#Modelling------------------------
#Transform data to distributions
speedDistribution <- function(trip){
  speed <- 3.6 * sqrt(diff(trip$x, 20, 1)^2 + diff(trip$y, 20, 1)^2) / 20
  return(list(quantile(speed, seq(0.05, 1, by=0.05)), speed))
}

#Define function to be passed as parallel
transform2Percentiles <- function(file, driverID){
  trip <- fread(file.path(driversDirectory, driverID, paste0(file, ".csv")))
  speedData <- speedDistribution(trip)
  speedDist <- speedData[[1]]
  accelerationDist <- quantile(diff(speedData[[2]]), seq(0.05, 1, by=0.05))
  distanceTrip <- log(sum(sqrt((diff(trip$x)^2) + (diff(trip$y)^2))))
  
  return(c(speedDist, accelerationDist, distanceTrip))
}

#List all possible drivers identities
drivers <- list.files(driversDirectory)
numCores <- detectCores() 

#Init h2o Server
#If there is need to start h2o from command line:
#system(paste0("java -Xmx55G -jar ", h2o.jarLoc, " -port 54333 -name CTR -data_max_factor_levels 100000000 &"))
#Start h2o directly from R
h2oServer <- h2o.init(ip = "localhost", max_mem_size = "5g", nthreads = -1) 

driversProcessed <- sapply(drivers, function(driver){
  #Parallel processing of each driver data
  results <- unlist(mclapply(seq(1, 200), transform2Percentiles, mc.cores = numCores, driverID = driver))
  results <- matrix(results, nrow = 200, byrow = TRUE)
  #KSVMmodel <- ksvm(results, type = "one-svc", kernel="rbfdot", kpar="automatic", prob.model = TRUE)
  #KSVMPrediction <- predict(KSVMmodel, results, type = "response")
    
  #Plotting
  #par(mfrow=c(1, 2))
  #plot(density(lofactor(results, k=5)))
  #biplot(prcomp(results), cex=.8) 
  
  #R matrix conversion to h2o object and stored in the server
  h2oResult <- as.h2o(h2oServer, cbind(rep(c(0, 1, 1, 0), 50), results))
  print(h2o.ls(h2oServer))
  driverDeepNNModel <- h2o.deeplearning(x = seq(2, ncol(h2oResult)), y = 1, activation = "Tanh",
                                        data = h2oResult, autoencoder = TRUE,
                                        hidden = c(11, 11), epochs = 750)
  anomalousTrips <- as.data.frame(h2o.anomaly(h2oResult, driverDeepNNModel))
  print(h2o.ls(h2oServer))
  h2o.rm(object = h2oServer, keys = h2o.ls(h2oServer)[, 1])  
  print(paste0("Driver number ", driver, " processed"))
  return(anomalousTrips[, 1])
})

h2o.shutdown(h2oServer, prompt = FALSE)  

#Predictions on the MSE scale turned into pseudo-probabilities
driversProcessed <- 1 - as.vector(driversProcessed)

#Write .csv------------------------- 
submissionTemplate <- fread(file.path(otherDataDirectory, "sampleSubmission.csv"), header = TRUE, 
                            stringsAsFactors = FALSE, colClasses = c("character", "numeric"))

submissionTemplate$prob <- signif(driversProcessed, digits = 4)
write.csv(submissionTemplate, file = gzfile("SpeedNNMSEPrediction.csv.gz"), row.names = FALSE)
