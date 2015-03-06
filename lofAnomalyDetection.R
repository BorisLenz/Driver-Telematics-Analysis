lofAnomalyDetection <- function(fullDataMatrix, cleanDataMatrix, kUser = 8, exclude_bad_data = TRUE){
  require("DMwR")
  
  if (exclude_bad_data == TRUE){
    goodDataIdx <- which(fullDataMatrix[, ncol(fullDataMatrix)] == 0)
    resultsLOF <- signif(-(lofactor(cleanDataMatrix[, -ncol(cleanDataMatrix)], k=kUser)) + 12, digits = 4)
    resultsLOFRank <- rank(resultsLOF)
    lofDriver <- rep(max(resultsLOF), 200)
    lofDriver[goodDataIdx] <- resultsLOF
    lofDriverRanking <- rep(max(resultsLOFRank), 200)
    lofDriverRanking[goodDataIdx] <- resultsLOFRank  
    return(list(lofDriver, lofDriverRanking))
  }
  if (exclude_bad_data == FALSE){
    resultsLOF <- signif(-(lofactor(fullDataMatrix[, -ncol(fullDataMatrix)], k=kUser)) + 12, digits = 4)
    resultsLOFRank <- rank(resultsLOF)
    return(list(resultsLOF, resultsLOFRank))
  }
  print(paste0("Driver number ", driver, " processed with the LOF Algorithm"))  
}