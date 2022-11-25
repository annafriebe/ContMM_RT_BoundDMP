# -*- coding: utf-8 -*-
"""
Created on Wed Sep  1 13:37:22 2021

@author: afe02
"""
import numpy as np
import csv
import EstZeroWLProbAccum


nStates = 2
transitionMatrix = np.zeros((nStates,nStates))
transitionMatrix[0, 0] = 0.9
transitionMatrix[0, 1] = 0.1
transitionMatrix[1, 0] = 0.7
transitionMatrix[1, 1] = 0.3

means = np.array([20, 40])
vars = np.array([9, 16]) 
statProbs = np.array([0.875, 0.125])
Qs = 8
n = 4
nQs = n* Qs
k = 8
deadline = k * Qs

beta_first_bound = np.array([0.1278, 0.0442])

epsilon = 1e-9

(cvProbsPeriods, zwl_LE, zwl_UE, lp_UE, result_betas, result_zwl_lo, result_zwl_hi, result_dmp_hi) = \
    EstZeroWLProbAccum.createAccVectorsZWLUBLinEqBound(nStates, means, vars, transitionMatrix, nQs, statProbs, beta_first_bound, deadline, epsilon)

with open('../../data/dmp_bound/betas_pzwl_bound_example.csv', 'w', newline='') as csvfile:
    betasZwlWriter = csv.writer(csvfile)
    for i in range(len(result_betas)):
        betasZwlWriter.writerow([i+1, result_betas[i][0], result_betas[i][1], \
                                result_zwl_hi[i][0], result_zwl_hi[i][1], \
                                result_zwl_lo[i][0], result_zwl_lo[i][1], \
                                result_dmp_hi[i][0][0], result_dmp_hi[i][0][1], \
                                result_dmp_hi[i][0][0] - result_betas[i][0]/statProbs[0], 
                                result_dmp_hi[i][0][1] - result_betas[i][1]/statProbs[1], \
                                result_dmp_hi[i][1]])

