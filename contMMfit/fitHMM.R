#!/usr/bin/env Rscript
lib <- modules::use("R")

library(depmixS4)
library(ggplot2)
library(data.tree)

EstimateModel <- function(dataFrame, outputDir, maxNStates, nPartitions, prevLikelihood){
  likelihoodsAndTree  <- lib$evalLikelihood$CrossValidationLikelihoodsAndClusteredTree(dataFrame, maxNStates, nPartitions)
  likelihoodsNCluster <- likelihoodsAndTree[[1]]
  plot(1:length(likelihoodsNCluster), likelihoodsNCluster)
  tree <- likelihoodsAndTree[[2]]
  if(likelihoodsNCluster[1] == 0){
    return(0)
  }
  # save to file if first or better
  lastLikelihood <- likelihoodsNCluster[length(likelihoodsNCluster)]
  if ((prevLikelihood > -1e-5) | (lastLikelihood > prevLikelihood)){
    fittedMod <- lib$evalLikelihood$FitMarkovChainFromClusteredTree(dataFrame, tree)
    nStates <- tree$leafCount
    transitionMatrix <- matrix(nrow=nStates, ncol=nStates)
    # get the transition matrix
    for (r in 1:nStates){
      transitionRow = fittedMod@transition[[r]]@parameters$coefficients
      for (c in 1:nStates){
        transitionMatrix[r, c]<- transitionRow[c]
      }
    }
    filename <- paste(outputDir, "/transitionMatrix.txt", sep="")
    write.table(transitionMatrix, file=filename, row.names=FALSE, col.names=FALSE)
    #get the normal distribution mean and stddev
    normalParams <- matrix(nrow=nStates, ncol=2)
    for (r in 1:nStates){
      normalParams[r, 1] <- fittedMod@response[[r]][[1]]@parameters$coefficients[[1]]
      normalParams[r, 2] <- fittedMod@response[[r]][[1]]@parameters$sd
    }
    filename <- paste(outputDir, "/normalParams.txt", sep="")
    write.table(normalParams, file=filename, row.names=FALSE, col.names=FALSE)
    # get the stationary distribution
    # Get the eigenvectors of P, note: R returns right eigenvectors
    r=eigen(transitionMatrix)
    rvec=r$vectors
    # left eigenvectors are the inverse of the right eigenvectors
    lvec=ginv(r$vectors)
    # normalized is the stationary distribution
    pi_eig<-lvec[1,]/sum(lvec[1,])
  
    filename <- paste(outputDir, "/stationaryDistr.txt", sep="")
    write.table(Re(pi_eig), file=filename, row.names=FALSE, col.names=FALSE)
    return(lastLikelihood)
  } 
  return(prevLikelihood) 
}

EstimateValidateModel <- function(dataFrame, outputDir, maxNStates, nPartitions, prevLikelihood){
  likelihoodsAndTree  <- lib$evalLikelihood$CrossValidationLikelihoodsAndClusteredTree(dataFrame, maxNStates, nPartitions)
  likelihoodsNCluster <- likelihoodsAndTree[[1]]
  plot(1:length(likelihoodsNCluster), likelihoodsNCluster)
  tree <- likelihoodsAndTree[[2]]
  if(likelihoodsNCluster[1] == 0){
    return(0)
  }
  # save to file if first or better
  lastLikelihood <- likelihoodsNCluster[length(likelihoodsNCluster)]
  if ((prevLikelihood > -1e-5) | (lastLikelihood > prevLikelihood)){
    fittedMod <- lib$evalLikelihood$FitMarkovChainFromClusteredTree(dataFrame, tree)
    nStates <- tree$leafCount
    transitionMatrix <- matrix(nrow=nStates, ncol=nStates)
    # get the transition matrix
    for (r in 1:nStates){
      transitionRow = fittedMod@transition[[r]]@parameters$coefficients
      for (c in 1:nStates){
        transitionMatrix[r, c]<- transitionRow[c]
      }
    }
    filename <- paste(outputDir, "/transitionMatrix.txt", sep="")
    write.table(transitionMatrix, file=filename, row.names=FALSE, col.names=FALSE)
    #get the normal distribution mean and stddev
    normalParams <- matrix(nrow=nStates, ncol=2)
    for (r in 1:nStates){
      normalParams[r, 1] <- fittedMod@response[[r]][[1]]@parameters$coefficients[[1]]
      normalParams[r, 2] <- fittedMod@response[[r]][[1]]@parameters$sd
    }
    filename <- paste(outputDir, "/normalParams.txt", sep="")
    write.table(normalParams, file=filename, row.names=FALSE, col.names=FALSE)
    # get the stationary distribution
    # Get the eigenvectors of P, note: R returns right eigenvectors
    r=eigen(transitionMatrix)
    rvec=r$vectors
    # left eigenvectors are the inverse of the right eigenvectors
    lvec=ginv(r$vectors)
    # normalized is the stationary distribution
    pi_eig<-lvec[1,]/sum(lvec[1,])
    
    filename <- paste(outputDir, "/stationaryDistr.txt", sep="")
    write.table(Re(pi_eig), file=filename, row.names=FALSE, col.names=FALSE)
    dfsim <- lib$dataConsistencyModelValidation$simulateTrajectory(fittedMod)
    dfsim$index = seq(1, length(dfsim$outputs))
    p <- ggplot(dfsim, aes(x=index, y=outputs)) +
      geom_point()
    print(p)
    filename <- paste(outputDir, "/simExecutionTimeSeq.eps", sep = "")
    ggsave(filename)
    filename <- paste(outputDir, "/simExecutionTimeSeq.png", sep = "")
    ggsave(filename)
    p <- ggplot(dfsim, aes(x=outputs)) +
      geom_histogram(binwidth = 100) 
    print(p)
    filename <- paste(outputDir, "/simExecutionTimeHis.eps", sep = "")
    ggsave(filename)
    filename <- paste(outputDir, "/simExecutionTimeHist.png", sep = "")
    ggsave(filename)
    # simulate data from the fitted model
    # first for estimation of mean and variance
    MP = 100
    M = 100
    zSim1 <- lib$dataConsistencyModelValidation$simZStatesTraj(fittedMod, MP, transitionMatrix, normalParams)
    Ez <- lapply(zSim1, apply, 1, mean)
    Vz <- lapply(zSim1, apply, 1, var)
    # then for consistency comparison
    zSim2 <- lib$dataConsistencyModelValidation$simZStatesTraj(fittedMod, M, transitionMatrix, normalParams)
    TSim <- list()
    for (i in 1:(nStates+1)){
      tmp <- (zSim2[[i]] - Ez[[i]])^2/Vz[[i]]
      TSim[[i]] <- apply(tmp, 2, mean)
    }
    print("get observation likelihoods")
    zObs <- lib$dataConsistencyModelValidation$logLikelihoodStatesTraj(dfsim, fittedMod, transitionMatrix, normalParams)
    TObs <- list()
    for (i in 1:(nStates+1)){
      TObs[[i]] <- mean((zObs[,i] - Ez[[i]])^2/Vz[[i]])
    }
    betaSim <- numeric(nStates + 1)
    for (j in 1:(nStates+1)){
      betaSim[j] = length(which(TSim[[j]] > TObs[[j]]))/ M
    }
    betas <- matrix(1, ncol=nStates+1)
    zObs <- lib$dataConsistencyModelValidation$logLikelihoodStatesTraj(dataFrame, fittedMod, transitionMatrix, normalParams)
    TObs <- list()
    for(j in 1:(nStates+1)){
      TObs[[j]] <- mean((zObs[,j] - Ez[[j]])^2/Vz[[j]])
    }
    for (j in 1:(nStates+1)){
      betas[1, j] = length(which(TSim[[j]] > TObs[[j]]))/ M
    }
    zObs <- lib$dataConsistencyModelValidation$logLikelihoodStatesTraj(dataFrame, fittedMod, transitionMatrix, normalParams)
    TObs <- list()
    for(j in 1:(nStates+1)){
      TObs[[j]] <- mean((zObs[,j] - Ez[[j]])^2/Vz[[j]])
    }
    filename <- paste(outputDir, "/PFAuTable.txt", sep = "")
    write.table(betas, file=filename, row.names=FALSE, col.names=FALSE, sep=" & ")
    
    return(lastLikelihood)
    
  } 
  return(prevLikelihood) 
}


#args = commandArgs(trailingOnly=TRUE)
#if (length(args) < 2){
#  stop("At least two input argument must be given (observations file and output directory)", call.=FALSE)
#}
#if(length(args) < 3){
  # default starting number of states
#  args[3] = 10
#}
#if(length(args) < 4){
  # default number of partitions
#  args[4] = 4
#}
#if(length(args) < 5){
  # default number of iterations to try
#  args[5] = 3
#}
observationsFile = "../data/traces_reports_csv/control_time.csv"
outputDir = "../data/models/cont"
maxNStates = 10
nPartitions = 4
nIters = 3
dataFrame1 = read.csv(observationsFile)
dataFrame1 = lib$importData$AdaptDataFrame(dataFrame1, 1, TRUE)


set.seed(42)

#likelihood = EstimateValidateModel(dataFrame1, outputDir, maxNStates, nPartitions, likelihood)

likelihood = 0
for (i in 1:nIters){
  likelihood = EstimateModel(dataFrame1, outputDir, maxNStates, nPartitions, likelihood)
  print(likelihood)
}




