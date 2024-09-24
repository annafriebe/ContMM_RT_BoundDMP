# -*- coding: utf-8 -*-
"""

@author: Anna Friebe
"""
import numpy as np
import csv
import EstZeroWLProbAccum
import time

def createOutputLine(i, result_betas, result_zwl_hi, result_zwl_lo, 
                     result_dmp_hi, nStates):
    outputList = []
    outputList.append(i+1)
    for s in range(nStates):
        outputList.append(result_betas[i][s])
    for s in range(nStates):
        outputList.append(result_zwl_hi[i][s])
    for s in range(nStates):
        outputList.append(result_zwl_lo[i][s])
    for s in range(nStates):
        outputList.append(result_dmp_hi[i][0][s])
    outputList.append(result_dmp_hi[i][1])
    return outputList
    

def boundAllConfigs(outputBaseFilename, QsList, nList, ksList, means, \
                       vars, transitionMatrix, statProbs, beta_first_bound_list, nStates, epsilon):
    nPeriods = [5, 10]
    for p in nPeriods:
        times = []
        for i in range(len(QsList)):
            Qs = QsList[i]
            n = nList[i]
            nQs = n * Qs
            ks = ksList[i]
            beta_first_bound = beta_first_bound_list[i]
            for k in ks:
                deadline = k*Qs
                min_dmp_overall = 1
                min_dmp_3 = 1
                t_before_bound = time.perf_counter_ns()
                (cvProbsPeriods, zwl_LE, zwl_UE, lp_UE, result_betas, result_zwl_lo, result_zwl_hi, result_dmp_hi) = \
                    EstZeroWLProbAccum.createAccVectorsZWLUBLinEqBound(nStates, means, vars, transitionMatrix, nQs, statProbs, beta_first_bound, deadline, epsilon, p)
                t_after_bound = time.perf_counter_ns()
                print("time", t_after_bound -t_before_bound)
                times.append(t_after_bound - t_before_bound)
                if p == nPeriods[-1]:
                    filename = outputBaseFilename + "_bound_" + str(Qs) + "_" + str(n) + "_" + str(k) + ".csv"
                    with open(filename, 'w', newline='') as csvfile:
                        betasZwlWriter = csv.writer(csvfile)
                        for i in range(len(result_betas)):
                            if result_dmp_hi[i][1] < min_dmp_overall:
                                min_dmp_overall = result_dmp_hi[i][1]
                            if result_dmp_hi[i][0][2] < min_dmp_3:
                                min_dmp_3 = result_dmp_hi[i][0][2]
                            outputLine = createOutputLine(i, result_betas, result_zwl_hi, result_zwl_lo, result_dmp_hi, nStates)
                            betasZwlWriter.writerow(outputLine)
                        betasZwlWriter.writerow((min_dmp_3, min_dmp_overall))
        filename = outputBaseFilename + str(p) + "_time.csv"
        time_arr = np.asarray(times)
        with open(filename, 'w', newline='') as csvfile:
            timesWriter = csv.writer(csvfile)
            for t in times:
                timesWriter.writerow((t,))
            timesWriter.writerow((np.mean(time_arr)*1e-9,))
            timesWriter.writerow((np.std(time_arr)*1e-9,))



nStates = 8
transitionMatrix = np.zeros((nStates,nStates))
with open('../../data/models/cont/transitionMatrix.txt', 'r') as csvfile:
    reader = csv.reader(csvfile)
    r = 0
    for row in reader:
        rowElems = row[0].split(sep = ' ')
        for c in range(nStates):
            transitionMatrix[r][c] = float(rowElems[c])
        r += 1

statProbs = np.zeros(nStates)
with open('../../data/models/cont/stationaryDistr.txt', 'r') as csvfile:
    reader = csv.reader(csvfile)
    r = 0
    for row in reader:
        rowElems = row[0].split(sep = ' ')
        statProbs[r] = float(rowElems[0])
        r+= 1

means = np.zeros(nStates)
variances = np.zeros(nStates)
with open('../../data/models/cont/normalParams.txt', 'r') as csvfile:
    reader = csv.reader(csvfile)
    r = 0
    for row in reader:
        rowElems = row[0].split(sep = ' ')
        means[r] = float(rowElems[0])
        variances[r] = pow(float(rowElems[1]), 2)
        r+= 1
        
QsList = [60000, 70000, 80000]
nList = [5, 4, 4]
ksList = [[8, 10], [6,8], [6,8]]

beta_first_bound_list = [np.array([0.000103, 0.001973, 0.003312, 0.000106, 0.000631, 0.000258, 0.000141, 0.000030]),
                         np.array([0.000157, 0.002259, 0.003648, 0.000185, 0.001354, 0.000303, 0.000197, 0.000066]),
                         np.array([0.000041, 0.001596, 0.002748, 0.000057, 0.000301, 0.000201, 0.000076, 0.000005])]

epsilon = 1e-9

outputBaseFilename = "../../data/dmp_bound/pzwl_dmp_control"

boundAllConfigs(outputBaseFilename, QsList, nList, ksList, means, variances, transitionMatrix, statProbs, beta_first_bound_list, nStates, epsilon)


