library(MASS)
library(ggplot2)

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

countDMSim <- function(startingState, n, transitionMatrix, means, stddevs, deadline, nQs){
  currentState = startingState
  workload = 0
  stateCount = rep(0, 2)
  lengthsStatesCount = matrix(data = rep(0, 60), nrow = 30, ncol = 2)
  lengthsCount = rep(0, 30)
  deadlineMissCountState = rep(0,2)
  nWLMax = 0
  nWL = 1
  zeroWLCount = 0
  zeroWLCountPerState = rep(0, 2)
  for(i in 1:n){
    stateProbs = transitionMatrix[currentState,]
    currentState = sample(1:2, size=1, prob = stateProbs)
    stateCount[currentState] = stateCount[currentState]+1
    lengthsCount[nWL] = lengthsCount[nWL]+1
    lengthsStatesCount[nWL,currentState] = lengthsStatesCount[nWL,currentState] + 1
    
    et = max(0, rnorm(1, mean= means[currentState], sd = stddevs[currentState]))
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
  indLast = 30
  while(lengthsProportions[indLast] < 1e-10){
    indLast = indLast - 1
  }
  beta_est = rep(0, indLast)
  beta_states_est = matrix(data = rep(0, 3*indLast), nrow = indLast, ncol = 3)
  beta_last = 0
  beta_states_last = rep(0, 2)
  while(indLast > 0){
    beta_last = beta_last + lengthsProportions[indLast]
    beta_est[indLast] = beta_last
    beta_states_last = beta_states_last + lengthsStatesProportions[indLast, ]
    beta_states_est[indLast,1] = indLast - 1
    beta_states_est[indLast,2:3] = beta_states_last
    indLast = indLast - 1
  }
  dmp = deadlineMissCountState/stateCount
  dmpAll = sum(deadlineMissCountState)/n
  zwlProbAll = zeroWLCount/n
  zwlProb = zeroWLCountPerState/stateCount
  return(list(dmp, dmpAll, beta_states_est, zwlProb, zwlProbAll))
}
  


nStates = 2
transitionMatrix = matrix(c(0.9, 0.7, 0.1, 0.3),
                          nrow=nStates, ncol=nStates)
means = c(1, 2)
stddevs = c(0.25, 0.5)
Qs = 1
n = 2
nQs = n* Qs
k = 4
deadline = k * Qs
nPeriods = 1000000

dmp_1_idle = 1 - pnorm(4, mean=1, sd=0.5)
dmp_2_idle = 1 - pnorm(4, mean=2, sd=1)
print(dmp_1_idle)
print(dmp_2_idle)

set.seed(3)
calcStationaryProbDistr(transitionMatrix)

resDMCountSimpleEx = countDMSim(2, nPeriods, transitionMatrix, means, stddevs, deadline, nQs)

print(resDMCountSimpleEx[[1]])
print(resDMCountSimpleEx[[2]])
print(resDMCountSimpleEx[[3]])
print(resDMCountSimpleEx[[4]])
print(resDMCountSimpleEx[[5]])
write.table(resDMCountSimpleEx[3], row.names = FALSE, col.names = FALSE, sep=",", "../../data/sim_cont/simpleExampleBetasSim.csv")

