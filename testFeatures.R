#small test of the code *do not include in the git*
#This piece of code is useful to test feature creation in every trip.
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
system(paste0("java -Xmx5G -jar ", h2o.jarLoc, " -port 54353 -name AXA &"))
#Connect R to h2o
h2oServer <- h2o.init(ip = "localhost", port = 54353, nthreads = numCores - 1)  
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
  h2o.rm(object = h2oServer, keys = h2o.ls(h2oServer)[, 1])   
  print(paste0(which(list.files(driversDirectory) == driver), "/", length(list.files(driversDirectory))))
  
  if (SpotInstance == TRUE){
    write.csv(cbind(lofDriver, lofDriverRanking), 
              file = file.path(outputDirectory, paste0(driver, ".csv")), row.names = FALSE)
    return(TRUE)
    
  }else{
    return(cbind(lofDriver, lofDriverRanking))
  }  
})

h2o.shutdown(h2oServer, prompt = FALSE)
