pretrainh2oNN <- function(numberOfDrivers){
  
  require("data.table")
  require("parallel")
  require("h2o")
  
  #Start h2o directly from R
  h2oServer <- h2o.init(ip = "localhost", max_mem_size = "5g", nthreads = -1)
  
  #Parallel processing of each driver data
  pretrainingDrivers <- sample(drivers, 1)
  
  results <- unlist(mclapply(seq(1, 200), transform2Percentiles, mc.cores = numCores, driverID = driver))
  results <- scale(matrix(results, nrow = 200, byrow = TRUE))
  print(paste0("Driver number ", driver, " processed"))
  
  numberOfExtraDrivers <- 2
  initialDrivers <- sample(drivers[!drivers %in% driver], numberOfExtraDrivers)
  ExtraDrivers <- sapply(initialDrivers, function(driver){
    results <- unlist(mclapply(seq(1, 200), transform2Percentiles, mc.cores = numCores, driverID = driver))
    print(paste0("Driver number: ", driver, " processed"))
    return(results)
  })
  ExtraDrivers <- scale(matrix(unlist(ExtraDrivers), nrow = numberOfExtraDrivers * 200, byrow = TRUE))  
  ExtraDrivers <- ExtraDrivers[sample(seq(1, nrow(ExtraDrivers)), 50), ]
  
  h2oResultPlusExtras <- as.h2o(h2oServer, cbind(c(rep(1, nrow(results)), rep(0, nrow(ExtraDrivers))), 
                                                 rbind(results, ExtraDrivers)))
  
  driverDeepNNModel <- h2o.deeplearning(x = seq(2, ncol(h2oResultPlusExtras)), y = 1,
                                        activation = "Tanh",
                                        data = h2oResultPlusExtras, fast_mode = TRUE,
                                        classification = TRUE,
                                        input_dropout_ratio = c(0, 0.1),
                                        hidden_dropout_ratios = list(c(0, 0, 0), c(0.2, 0.2, 0.2)),                                      
                                        l1 = c(0, 1e-5),
                                        l2 = c(0, 1e-5),
                                        rho = c(0.95, 0.99),
                                        epsilon = c(1e-10, 1e-8),
                                        hidden = c(50, 40, 50), 
                                        epochs = 5)
  
  checkpointModel <- driverDeepNNModel@model[[1]] #Best NN cv model
  deepNetPath <- h2o.saveModel(object = checkpointModel, dir = file.path(workingDirectory), force = TRUE)
  h2o.shutdown(h2oServer, prompt = FALSE)
  #save(deepNetPath, file = "deepNetPath.RData")
  
  #Start h2o directly from R
  h2oServer <- h2o.init(ip = "localhost", max_mem_size = "5g", nthreads = -1)
  
  #load(file.path(workingDirectory, "deepNetPath.RData"))
  checkpointModel <- h2o.loadModel(h2oServer, deepNetPath)
  checkpointModelKey <- h2o.ls(h2oServer)[, 1]
  
  #Begin with a randomly selected list of drivers to pre-train the model
  numberOfDrivers <- numberOfDrivers
  pretrainingDrivers <- sample(drivers, numberOfDrivers)
  
  driversPredictions <- lapply(pretrainingDrivers, function(driver){
    #Parallel processing of each driver data
    results <- unlist(mclapply(seq(1, 200), transform2Percentiles, mc.cores = numCores, driverID = driver))
    results <- scale(matrix(results, nrow = 200, byrow = TRUE))
    print(paste0("Driver number ", driver, " processed"))
    
    #Sample data from other drivers  
    numberOfExtraDrivers <- 2
    initialDrivers <- sample(drivers[!drivers %in% driver], numberOfExtraDrivers)
    ExtraDrivers <- sapply(initialDrivers, function(driver){
      results <- unlist(mclapply(seq(1, 200), transform2Percentiles, mc.cores = numCores, driverID = driver))
      print(paste0("Driver number: ", driver, " processed"))
      return(results)
    })
    ExtraDrivers <- scale(matrix(unlist(ExtraDrivers), nrow = numberOfExtraDrivers * 200, byrow = TRUE))  
    ExtraDrivers <- ExtraDrivers[sample(seq(1, nrow(ExtraDrivers)), 50), ]
    
    #h20.ai deep learning algorithm
    #R matrix conversion to h2o object and stored in the server
    h2oResultPlusExtras <- as.h2o(h2oServer, cbind(c(rep(1, nrow(results)), rep(0, nrow(ExtraDrivers))), 
                                                   rbind(results, ExtraDrivers)))
    h2oResultPlusExtras <- h2o.splitFrame(h2oResultPlusExtras, ratios = 0.99, shuffle = TRUE)[[1]]
    
    print(h2o.ls(h2oServer))
    driverDeepNNModel <- h2o.deeplearning(x = seq(2, ncol(h2oResultPlusExtras)), y = 1,
                                          activation = "Tanh",
                                          data = h2oResultPlusExtras, fast_mode = TRUE,
                                          classification = TRUE,
                                          checkpoint = checkpointModel,
                                          input_dropout_ratio = c(0, 0.1),
                                          hidden_dropout_ratios = list(c(0, 0, 0), c(0.2, 0.2, 0.2)),                                        
                                          l1 = c(0, 1e-5),
                                          l2 = c(0, 1e-5),
                                          rho = c(0.95, 0.99),
                                          epsilon = c(1e-10, 1e-8),
                                          hidden = c(50, 40, 50), 
                                          epochs = 5)
    
    checkpointModel <- driverDeepNNModel@model[[1]] #Best NN cv model
    print(h2o.ls(h2oServer))
    h2oObjects2Remove <- which(!h2o.ls(h2oServer)[, 1] %in% "NewCheckpointModel")
    h2o.rm(object = h2oServer, keys = h2o.ls(h2oServer)[h2oObjects2Remove, 1]) 
    print(paste0(which(drivers == driver), "/", numberOfDrivers))
    
    return(TRUE)
  })
  
  deepNetPath <- h2o.saveModel(object = checkpointModel, dir = file.path(workingDirectory), force = TRUE)
  #save(deepNetPath, file = "deepNetPath.RData")
  
  h2o.shutdown(h2oServer, prompt = FALSE)
  
}
