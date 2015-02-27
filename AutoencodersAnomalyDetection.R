AutoencodersAnomalyDetection <- function(dataMatrix, nwServer){
  #Required libraries
  require("h2o")
  
  #Bad Data Indices 
  goodDataIdx <- which(!dataMatrix[, ncol(dataMatrix)] == 1)
  #Generate Random indexes
  randIdxs <- sample(goodDataIdx, length(goodDataIdx))
  
  #Good Data Shuffled
  goodData <- cbind(rbinom(length(randIdxs), 1, 0.5), signif(dataMatrix[randIdxs, -ncol(dataMatrix)], digits = 4))
  fullData <- cbind(rbinom(nrow(dataMatrix), 1, 0.5), signif(dataMatrix[, -ncol(dataMatrix)], digits = 4))
  h2oMatrix <- as.h2o(nwServer, goodData)
  allDataH2o <- as.h2o(nwServer, fullData)
  
  #Make a deep learning autoencoder model
  ae_model <- h2o.deeplearning(x = seq(2, ncol(h2oMatrix)), y = 1,
                               data = h2oMatrix,   
                               autoencoder = TRUE, 
                               ignore_const_cols = FALSE,
                               activation = "Tanh",
                               hidden = c(60, 60), 
                               epochs = 100)
  
  #Anomaly Detection
  test_rec_error <- as.data.frame(h2o.anomaly(allDataH2o, ae_model))
  h2o.rm(object = nwServer, keys = h2o.ls(nwServer)[, 1])    
    
  return(test_rec_error[, 1])  
}