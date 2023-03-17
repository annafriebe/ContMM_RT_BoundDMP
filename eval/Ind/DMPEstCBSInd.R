#!/usr/bin/env Rscript

countDMTest <- function(executionTimes, deadline, nQs){
  workload = 0
  maxLength = 500
  lengthsCount = rep(0, maxLength)
  deadlineMissCount = 0
  nWLMax = 0
  nWL = 1
  zeroWLCount = 0
  deadlineMissAt = c()
  for(i in 1:length(executionTimes)){
    et = executionTimes[i]
    workload = workload + et
    lengthsCount[nWL] = lengthsCount[nWL] + 1
    if (workload > deadline){
      deadlineMissCount = deadlineMissCount + 1
      deadlineMissAt = c(deadlineMissAt, i)
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
    }
    workload = max(0, workload-nQs)
  }
  zWLRatio = zeroWLCount/length(executionTimes)
  dmpAll = deadlineMissCount/length(executionTimes)
  return(list(zWLRatio, nWLMax, dmpAll))
}


testDMPAll <- function(QsList, nsList, ksList, filename, nFiles, nSamplesInd){
  outputQ = c()
  outputN = c()
  outputK = c()
  outputDMP = c()
  for(j in 1:3){
    Qs = QsList[j]
    n = nsList[j]
    nQs = n*Qs
    for(kiter in 1:2){
      k = ksList [[j]][kiter]
      deadline = k*Qs
      print(Qs)
      print(n)
      print(k)
      dataFrame <- read.csv(filename, header=TRUE)
      independentOrder = sample(dataFrame$executionTime, nSamplesInd, replace=TRUE)
      resDMCountSimpleEx = countDMTest(independentOrder, deadline, nQs)
      print(resDMCountSimpleEx[[3]])
      outputQ = append(outputQ, Qs)
      outputN = append(outputN, n)
      outputK = append(outputK, k)
      outputDMP = append(outputDMP, resDMCountSimpleEx[[3]])
    }
  }
  output = data.frame(Q = outputQ, n = outputN, k = outputK, dmp = outputDMP)
  return(output)
}

QsList = c(60000, 70000, 80000)
nsList = c(5, 4, 4)
ksList = list(c(8, 10), c(6, 8), c(6,8))

filename = "../../data/traces_reports_csv/control_time.csv"
set.seed(23)
print("results")
output = testDMPAll(QsList, nsList, ksList, filename, nSequences, 1000000)
outputFileName="../../data/ind/result.csv"
write.csv(output, outputFileName)
