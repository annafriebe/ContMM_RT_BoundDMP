library(MASS)
library(ggplot2)

set.seed(42)

calcStationaryProbDistr <- function(transitionMatrix){
  r=eigen(transitionMatrix)
  rvec=r$vectors
  # left eigenvectors are the inverse of the right eigenvectors
  lvec=ginv(r$vectors)
  # normalized is the stationary distribution
  pi_eig<-lvec[1,]/sum(lvec[1,])
  print("stationary prob")
  print(pi_eig)
}

countDMSim <- function(startingState, n, transitionMatrix, means, stddevs, deadline, nQs, nStates){
  currentState = startingState
  workload = 0
  stateCount = rep(0, nStates)
  lengthsStatesCount = matrix(data = rep(0, nStates*200), nrow = 200, ncol = nStates)
  lengthsCount = rep(0, 200)
  deadlineMissCountState = rep(0, nStates)
  nWLMax = 0
  nWL = 1
  zeroWLCount = 0
  zeroWLCountPerState = rep(0, nStates)
  for(i in 1:n){
    stateProbs = transitionMatrix[currentState,]
    currentState = sample(1:nStates, size=1, prob = stateProbs)
    stateCount[currentState] = stateCount[currentState]+1
    lengthsCount[nWL] = lengthsCount[nWL]+1
    lengthsStatesCount[nWL,currentState] = lengthsStatesCount[nWL,currentState] + 1

    et = rnorm(1, mean=means[currentState], sd=stddevs[currentState])
    workload = workload + et
    if (workload > deadline){
      deadlineMissCountState[currentState] = deadlineMissCountState[currentState] + 1
    }
    if (workload > nQs){
      nWL = nWL+1
      if (nWL > nWLMax){
        nWLMax = nWL
      }
    }
    else{
      nWL = 1
      zeroWLCount = zeroWLCount + 1
      zeroWLCountPerState[currentState] = zeroWLCountPerState[currentState] +1
    }
    workload = max(0, workload-nQs)
  }
  lengthsStatesProportions = lengthsStatesCount/sum(lengthsStatesCount)
  lengthsProportions = lengthsCount/sum(lengthsCount)
  statesProportions = colSums(lengthsStatesProportions)
  print("states proportions")
  print(statesProportions)
  indLast = 50
  while(lengthsProportions[indLast] < 1e-10){
    indLast = indLast - 1
  }
  beta_est = rep(0, indLast)
  beta_states_est = matrix(data = rep(0, (nStates+1)*indLast), nrow = indLast, ncol = (nStates+1))
  beta_last = 0
  beta_states_last = rep(0, nStates)
  while(indLast > 0){
    beta_last = beta_last + lengthsProportions[indLast]
    beta_est[indLast] = beta_last
    beta_states_last = beta_states_last + lengthsStatesProportions[indLast, ]
    beta_states_est[indLast,1] = indLast - 1
    beta_states_est[indLast,2:(nStates+1)] = beta_states_last
    indLast = indLast - 1
  }
  dmp = deadlineMissCountState/stateCount
  dmpAll = sum(deadlineMissCountState)/n
  return(list(dmp, dmpAll, beta_states_est))
}
  


nStates = 8
trMat = read.csv("../../data/models/cont/transitionMatrix.txt", 
                 header=FALSE, sep = " ")
transitionMatrix = data.matrix(trMat)
normalParams = read.csv("../../data/models/cont/normalParams.txt", 
                        header=FALSE, sep = " ")


means = normalParams[,1]
stddevs = normalParams[,2]


QsList = c(60000, 70000, 80000)
nsList = c(5, 4, 4)
ksList = list(c(8, 10), c(6, 8), c(6,8))

outputQ = c()
outputN = c()
outputK = c()
outputDMP3 = c()
outputDMPAll = c()


for(j in 1:3){
  Qs = QsList[j]
  n = nsList[j]
  nQs = n*Qs
  for(kiter in 1:2){
    k = ksList [[j]][kiter]
    deadline = k*Qs
    nPeriods = 1000000

    calcStationaryProbDistr(transitionMatrix)
    print(Qs)
    print(n)
    print(k)
    
    resDMCountSimpleEx = countDMSim(1, nPeriods, transitionMatrix, means, stddevs, deadline, nQs, nStates)
  
    print(resDMCountSimpleEx[[1]])
    print(resDMCountSimpleEx[[2]])
    print(resDMCountSimpleEx[[3]])
    filename = paste("../../data/sim_cont/betasSim", Qs, "_", n, "_", k, ".csv", sep="")
    write.table(resDMCountSimpleEx[3], row.names = FALSE, col.names = FALSE, sep=",", filename)
    outputQ = append(outputQ, Qs)
    outputN = append(outputN, n)
    outputK = append(outputK, k)
    outputDMP3 = append(outputDMP3, resDMCountSimpleEx[[1]][[3]])
    outputDMPAll = append(outputDMPAll, resDMCountSimpleEx[[2]])
  }
}

output = data.frame(Q = outputQ, n = outputN, k = outputK, dmp3 = outputDMP3, dmpAll = outputDMPAll)
outputFileName="../../data/sim_cont/dmp_result.csv"
write.csv(output, outputFileName)

