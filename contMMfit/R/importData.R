


export("AdaptDataFrame")
AdaptDataFrame <- function (dataFrame, resolution = 1, executionTimes=TRUE)  {
  colnames(dataFrame) <- c("executionTime")
  chainLength <- length(dataFrame$executionTime)
  dataFrame$indices <- seq(1, chainLength)
  dataFrame$outputs <- dataFrame$executionTime
  return(dataFrame)
}





