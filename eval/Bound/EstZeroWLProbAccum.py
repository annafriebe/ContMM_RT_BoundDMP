# -*- coding: utf-8 -*-
"""

@author: Anna Friebe
"""
import numpy as np
import math
from scipy.stats import norm
from StateCVProbsLimitsAlpha import StateCVProbs

    
def initCombinationVectorProbs(nStates, transitionMatrix, means, vars, nQs, statProbs):
    cvProbsFirstPer = []
    for s in range(nStates):
        cvProbsFirstPer.append({})
    for s in range(nStates):
        combinationVector = [0]*nStates
        combinationVector[s] = 1
        cvTuple = tuple(combinationVector)
        sCVProbs = StateCVProbs(s, cvTuple, nStates)
        for sPrev in range(nStates):
            prevOutputProbs = [0]*nStates
            prevOutputProbs[sPrev] = statProbs[sPrev]
            sCVProbs.addInputProbs(transitionMatrix, sPrev, prevOutputProbs, prevOutputProbs, 0, nQs)
        sCVProbs.calcMeanVarOutputProbs(means, vars, nQs)
        cvProbsFirstPer[s][cvTuple] = sCVProbs
    return cvProbsFirstPer

def initCombinationVectorProbsAlphas(nStates, transitionMatrix, means, vars, alphas, nQs, statProbs):
    cvProbsFirstPer = []
    for s in range(nStates):
        cvProbsFirstPer.append({})
    for s in range(nStates):
        combinationVector = [0]*nStates
        combinationVector[s] = 1
        cvTuple = tuple(combinationVector)
        sCVProbs = StateCVProbs(s, cvTuple, nStates, alphas[s]-nQs)
        for sPrev in range(nStates):
            prevOutputProbs = [0]*nStates
            prevOutputProbs[sPrev] = statProbs[sPrev]
            sCVProbs.addInputProbs(transitionMatrix, sPrev, prevOutputProbs, prevOutputProbs, 0, nQs)
        sCVProbs.calcMeanVarOutputProbs(means, vars, nQs)
        cvProbsFirstPer[s][cvTuple] = sCVProbs
    return cvProbsFirstPer

def addCombinationVectorProbs(nStates, transitionMatrix, means, vars, nQs, 
                              prevCombVectorProbs):
    cvProbsNextPer = []
    for s in range(nStates):
        cvProbsNextPer.append({})
    for sPrev in range(nStates):
        for cvPrev, sCVProbsPrev in prevCombVectorProbs[sPrev].items():
            cvList = list(cvPrev)
            for s in range(nStates):
                cvListNext = cvList.copy()
                cvListNext[s] += 1
                cvTuple = tuple(cvListNext)
                if cvTuple not in cvProbsNextPer[s]:
                    cvProbsNextPer[s][cvTuple] = \
                        StateCVProbs(s, cvTuple, nStates)
                sCVProbs = cvProbsNextPer[s][cvTuple]
                sCVProbs.addInputProbs(transitionMatrix, sPrev, \
                                       sCVProbsPrev.getOutputProbsLo(), \
                                       sCVProbsPrev.getOutputProbsHi(),
                                       sCVProbsPrev.getAlpha(), nQs)
    for s in range(nStates):
        for cv, sCVProbs in cvProbsNextPer[s].items():
            sCVProbs.calcMeanVarOutputProbs(means, vars, nQs)
    return cvProbsNextPer

def calcDMPHi(cvProbs, nPeriods, nStates, zwlUB, lp_UB, statProbs, nQs, deadline):
    dmpOE = np.zeros(nStates)
    sumProbsHi = np.zeros(nStates)
    probsPartHi = np.zeros(nStates)
    for s in range(nStates):
        stateSumProbsHi = np.zeros(nStates)
        dmpSumHi = 0
        for n in range(nPeriods):
            for cv, sCVProbs in cvProbs[n][s].items():
                inProbHi = sCVProbs.getInputProbsHi()
                factor = sum(inProbHi*zwlUB) 
                alpha = sCVProbs.getAlpha() + nQs
                dmpHi = 1
                if deadline > alpha:
                    dmpHi = norm.sf(deadline, sCVProbs.getMean() + nQs, \
                                    math.sqrt(sCVProbs.getVar())) * sCVProbs.getK()
                dmpSumHi += factor * dmpHi
                stateSumProbsHi += inProbHi
        sumProbsHi[s] = sum(stateSumProbsHi*zwlUB)
        probsPartHi[s] = sumProbsHi[s] / statProbs[s] 
        dmpOE[s] = dmpSumHi/statProbs[s] + lp_UB[s]/statProbs[s] 
#    print('probs part: ', probsPartHi)
#    print('dmp: ', dmpOE)
    dmpAll = sum(dmpOE*statProbs)
#    print('dmp overall: ', dmpAll)
    return (dmpOE, dmpAll)


def linEqSolveZWL_UB_maxAllDim(cvProbs, nPeriods, nStates, statProbs):
    estProbLoMatrix = np.zeros((nStates, nStates))
    for s in range(nStates):
        for n in range(nPeriods):
            for cv, sCVProbs in cvProbs[n][s].items():
                estProbLoMatrix[s, :] += sCVProbs.getInputProbsLo()
    zwlEst = np.linalg.solve(estProbLoMatrix, statProbs)
    return zwlEst

def linEqSolveZWL_UB_1DLimit(cvProbs, nPeriods, nStates, statProbs, longerProbs, dim):
    estProbMatrix = np.zeros((nStates, nStates))
    b = statProbs.copy()
    for s in range(nStates):
        if s == dim:
            for n in range(nPeriods):
                for cv, sCVProbs in cvProbs[n][s].items():
                    estProbMatrix[s, :] += sCVProbs.getInputProbsHi()
        else:
            for n in range(nPeriods):
                for cv, sCVProbs in cvProbs[n][s].items():
                    estProbMatrix[s, :] += sCVProbs.getInputProbsLo()
    b[dim] -= longerProbs[dim]
    zwlEst = np.linalg.solve(estProbMatrix, b)
    return zwlEst
        

def linEqSolveZWL_LE_minAllDim(cvProbs, nPeriods, nStates, statProbs, longerProbs):
    estProbLoMatrix = np.zeros((nStates, nStates))
    for s in range(nStates):
        for n in range(nPeriods):
            for cv, sCVProbs in cvProbs[n][s].items():
                estProbLoMatrix[s, :] += sCVProbs.getInputProbsHi()
    zwlEst = np.linalg.solve(estProbLoMatrix, statProbs - longerProbs)
    return zwlEst
        

def linEqSolveZWL_LE_1DLimit(cvProbs, nPeriods, nStates, statProbs, longerProbs, dim):
    estProbMatrix = np.zeros((nStates, nStates))
    b = statProbs - longerProbs
    for s in range(nStates):
        if s == dim:
            for n in range(nPeriods):
                for cv, sCVProbs in cvProbs[n][s].items():
                    estProbMatrix[s, :] += sCVProbs.getInputProbsLo()
        else:
            for n in range(nPeriods):
                for cv, sCVProbs in cvProbs[n][s].items():
                    estProbMatrix[s, :] += sCVProbs.getInputProbsHi()
    b[dim] = statProbs[dim]
    zwlEst = np.linalg.solve(estProbMatrix, b)
    return zwlEst
        


def accountedForProbsLo(cvProbs, nPeriods, nStates, statProbs, zwl):
    accountedForProbs_LE = np.zeros(nStates)
    for s in range(nStates):
        for n in range(nPeriods):
            for cv, sCVProbs in cvProbs[n][s].items():
                accountedForProbs_LE[s] += sum(sCVProbs.getInputProbsLo()*zwl)
    return accountedForProbs_LE

def accountedForProbsHi(cvProbs, nPeriods, nStates, statProbs, zwl):
    accountedForProbs_HI = np.zeros(nStates)
    for s in range(nStates):
        for n in range(nPeriods):
            for cv, sCVProbs in cvProbs[n][s].items():
                accountedForProbs_HI[s] += sum(sCVProbs.getInputProbsHi()*zwl)
    return accountedForProbs_HI

def longerProbs_UE(cvProbs, nPeriods, nStates, statProbs, zwl_LE):
    lpe = statProbs - accountedForProbsLo(cvProbs, nPeriods, nStates, statProbs, zwl_LE)
    for s in range(nStates):
        if lpe[s] < 0:
            lpe[s] = 0
    return lpe



def lastPeriodProbs(cvProbs, nPeriods, nStates, zwl):
    lpe = np.zeros(nStates)
    for s in range(nStates):
        probInLastScalar = 0
        for cv, sCVProbs in cvProbs[nPeriods-1][s].items():
            probInLastScalar += sum(sCVProbs.getInputProbsLo()*zwl)
        lpe[s] = probInLastScalar
    return lpe

def calcBeta(cvProbs, nPeriods, nStates, zwl, prevBeta, statProbs):
    lpe = lastPeriodProbs(cvProbs, nPeriods, nStates, zwl)
    beta = prevBeta - lpe
    betaTest = longerProbs_UE(cvProbs, nPeriods, nStates, statProbs, zwl)
    for s in range(nStates):
        if betaTest[s] < beta[s]:
            beta[s] = betaTest[s]
    return beta

# precondition start valid
def limitFactor(start, stop):
    if stop <= 1 and stop >= 0:
        return 1
    if stop > 1:
        return (1-start)/(stop-start)
    return start/(start-stop)

# precondition one valid
def limitFactors(start, stop):
    lf = []
    if start <= 1 and start >= 0:
        lf.append(0)
    else:
        f = limitFactor(stop, start)
        lf.append(1-f)
    if stop <= 1 and stop >= 0:
        lf.append(1)
    else:
        f = limitFactor(start, stop)
        lf.append(f)
    return lf


def lowerInvalid(point, nStates):
    for s in range(nStates):
        if point[s] < 0:
            return True
    return False

def higherInvalid(point, nStates):
    for s in range(nStates):
        if point[s] > 1:
            return True
    return False

def invalid(point, nStates):
    return lowerInvalid(point, nStates) or higherInvalid(point, nStates)

def invalidVector(start, stop, nStates):
    for s in range(nStates):
        if start[s] < 0 and stop[s] < 0:
            return True
        if start[s] > 1 and stop[s] > 1:
            return True
    return False

# precondition startVec invalid
def findLimitFactors(startVec, stopVec, nStates):
    lf = []
    if invalid(stopVec, nStates):
        if invalidVector(startVec, stopVec, nStates):
            return lf
    lf = [0, 1]
    for s in range(nStates):
        stateFactors = limitFactors(startVec[s], stopVec[s])
        if stateFactors[0] > lf[0]:
            lf[0] = stateFactors[0]
        if stateFactors[1] < lf[1]:
            lf[1] = stateFactors[1]
    return lf


def findMinLimitFactor(startVec, stopVec, nStates):
    fact = 1
    for s in range(nStates):
        stateFact = limitFactor(startVec[s], stopVec[s])
        if stateFact < fact:
            fact = stateFact
    return fact

def validEndpoints(points, nStates):
    vep = []
    if invalid(points[0], nStates):
        for i in range(1, len(points)):
            lf = findLimitFactors(points[0], points[i], nStates)
            for f in lf:
                vep.append(points[0] + f*(points[i] - points[0]))
    else:
        vep.append(points[0])
        for i in range(1, len(points)):
            fact = findMinLimitFactor(points[0], points[i], nStates)
            vep.append(points[0] + fact*(points[i] - points[0]))
    return vep


def maxProbs(probsList, nStates):
    mp = probsList[0].copy()
    for i in range(1, len(probsList)):
        for s in range(nStates):
            if probsList[i][s] > mp[s]:
                mp[s] = probsList[i][s]
    return mp


def findMaxProbs(cvProbs, nStates, statProbs, limits):
    pointsList = []
    zwl_UE = linEqSolveZWL_UB_maxAllDim(cvProbs, len(cvProbs), nStates, statProbs)
    pointsList.append(zwl_UE)
    for s in range(nStates):
        zwl_UE_test = linEqSolveZWL_UB_1DLimit(cvProbs, len(cvProbs), nStates, statProbs, limits, s)
        pointsList.append(zwl_UE_test)
    probsList = validEndpoints(pointsList, nStates)
    if not probsList:
#        print("No valid upper point")
        return np.ones(nStates)
    zwl_UB = maxProbs(probsList, nStates)
    return zwl_UB


def minProbs(probsList, nStates):
    mp = probsList[0].copy()
    for i in range(1, len(probsList)):
        for s in range(nStates):
            if probsList[i][s] < mp[s]:
                mp[s] = probsList[i][s]
    return mp


def findMinProbs(cvProbs, nStates, statProbs, limits):
    pointsList = []
    zwl_LE = linEqSolveZWL_LE_minAllDim(cvProbs, len(cvProbs), nStates, statProbs, limits)
    pointsList.append(zwl_LE)
    for s in range(nStates):
        zwl_LE_test = linEqSolveZWL_LE_1DLimit(cvProbs, len(cvProbs), nStates, statProbs, limits, s)
        pointsList.append(zwl_LE_test)
    probsList = validEndpoints(pointsList, nStates)
    if not probsList:
#        print("no valid lower point")
        return np.zeros(nStates)
    zwl_LB = minProbs(probsList, nStates)
    return zwl_LB


def createAccVectorsZWLUBLinEqBound(nStates, means, vars, transitionMatrix, nQs, \
                                  statProbs, beta_bound, deadline, epsilon, max_accum_periods = 20, period_bound = 1):
    cvProbsPeriods = []
    cvProbs = initCombinationVectorProbs(nStates, transitionMatrix, \
                                         means, vars, nQs, statProbs)
    cvProbsPeriods.append(cvProbs)
    for i in range(period_bound-1):
        cvProbsNext = addCombinationVectorProbs(nStates, transitionMatrix, \
                                                means, vars, nQs, cvProbs)
        cvProbsPeriods.append(cvProbsNext)
        cvProbs = cvProbsNext
    ubPLonger = beta_bound 
    result_betas = []
    result_zwl_lo = []
    result_zwl_hi = []
    result_dmp_hi = []
    result_betas.append(ubPLonger.copy())
    zwl_LE = findMinProbs(cvProbsPeriods, nStates, statProbs, ubPLonger)
    lp_UE = longerProbs_UE(cvProbsPeriods, 1, nStates, statProbs, zwl_LE)

    zwl_UE = findMaxProbs(cvProbsPeriods, nStates, statProbs, ubPLonger)
    dmp = calcDMPHi(cvProbsPeriods, len(cvProbsPeriods), nStates, zwl_UE, ubPLonger, statProbs, nQs, deadline)
    result_zwl_lo.append(zwl_LE.copy())
    result_zwl_hi.append(zwl_UE.copy())
    result_dmp_hi.append(dmp)

    ue_started_decrease = np.full((nStates), False)
    ue_stopped_decrease = np.full((nStates), False)
    le_started_increase = np.full((nStates), False)
    le_stopped_increase = np.full((nStates), False)
    eps_vec = np.full((nStates), 1e-9)

    while not (np.all(ue_stopped_decrease) or np.all(le_stopped_increase)) and (len(cvProbsPeriods) < max_accum_periods):
        cvProbsNext = addCombinationVectorProbs(nStates, transitionMatrix, \
                                                means, vars, nQs, cvProbs)
        cvProbsPeriods.append(cvProbsNext)
#        print(len(cvProbsPeriods))
        cvProbs = cvProbsNext
        prev_zwl_LE = zwl_LE.copy()
        prev_zwl_UE = zwl_UE.copy()
        prevBeta = ubPLonger.copy()
        ubPLonger = calcBeta(cvProbsPeriods, len(cvProbsPeriods), nStates, zwl_LE, prevBeta, statProbs)
        result_betas.append(ubPLonger.copy())
        zwl_LE = findMinProbs(cvProbsPeriods, nStates, statProbs, ubPLonger)
        lp_UE = longerProbs_UE(cvProbsPeriods, len(cvProbsPeriods), nStates, statProbs, zwl_LE)
        zwl_UE = findMaxProbs(cvProbsPeriods, nStates, statProbs, ubPLonger)
        dmp = calcDMPHi(cvProbsPeriods, len(cvProbsPeriods), nStates, zwl_UE, ubPLonger, statProbs, nQs, deadline)
        result_zwl_lo.append(zwl_LE.copy())
        result_zwl_hi.append(zwl_UE.copy())
        result_dmp_hi.append(dmp)
        ue_started_decrease = np.logical_or(zwl_UE < prev_zwl_UE - eps_vec, ue_started_decrease)
        ue_stopped_decrease = np.logical_and(np.logical_or(zwl_UE > prev_zwl_UE - eps_vec, ue_stopped_decrease), ue_started_decrease)
        le_started_increase = np.logical_or(zwl_LE > prev_zwl_LE + eps_vec, le_started_increase)
        le_stopped_increase = np.logical_and(np.logical_or(zwl_LE < prev_zwl_LE + eps_vec, le_stopped_increase), le_started_increase)
    #    print(ue_started_decrease)
    #    print(ue_stopped_decrease)
    #    print(le_started_increase)
    #    print(le_stopped_increase)
        
    return (cvProbsPeriods, zwl_LE, zwl_UE, ubPLonger, result_betas, result_zwl_lo, result_zwl_hi, result_dmp_hi)




def createAccVectorsZWLUBLinEqBoundMerge(nStates, means, vars, alphas, 
                                         transitionMatrix, nQs, statProbs, 
                                         beta_bound, deadline, epsilon, 
                                         max_accum_periods = 20, period_bound = 1):
    cvProbsPeriods = []
    cvProbs = initCombinationVectorProbsAlphas(nStates, transitionMatrix, 
                                               means, vars, alphas, nQs, statProbs)
    cvProbsPeriods.append(cvProbs)
    for i in range(period_bound-1):
        cvProbsNext = addCombinationVectorProbs(nStates, transitionMatrix, \
                                                means, vars, nQs, cvProbs)
        cvProbsPeriods.append(cvProbsNext)
        cvProbs = cvProbsNext
    ubPLonger = beta_bound 
    result_betas = []
    result_zwl_lo = []
    result_zwl_hi = []
    result_dmp_hi = []
    result_betas.append(ubPLonger.copy())
    zwl_LE = findMinProbs(cvProbsPeriods, nStates, statProbs, ubPLonger)
    lp_UE = longerProbs_UE(cvProbsPeriods, 1, nStates, statProbs, zwl_LE)

    zwl_UE = findMaxProbs(cvProbsPeriods, nStates, statProbs, ubPLonger)
    dmp = calcDMPHi(cvProbsPeriods, len(cvProbsPeriods), nStates, zwl_UE, ubPLonger, statProbs, nQs, deadline)
    result_zwl_lo.append(zwl_LE.copy())
    result_zwl_hi.append(zwl_UE.copy())
    result_dmp_hi.append(dmp)

    ue_started_decrease = np.full((nStates), False)
    ue_stopped_decrease = np.full((nStates), False)
    le_started_increase = np.full((nStates), False)
    le_stopped_increase = np.full((nStates), False)
    eps_vec = np.full((nStates), 1e-9)

    while not (np.all(ue_stopped_decrease) or np.all(le_stopped_increase)) and (len(cvProbsPeriods) < max_accum_periods):
        cvProbsNext = addCombinationVectorProbs(nStates, transitionMatrix, \
                                                means, vars, nQs, cvProbs)
        cvProbsPeriods.append(cvProbsNext)
#        print(len(cvProbsPeriods))
        cvProbs = cvProbsNext
        prev_zwl_LE = zwl_LE.copy()
        prev_zwl_UE = zwl_UE.copy()
        prevBeta = ubPLonger.copy()
        ubPLonger = calcBeta(cvProbsPeriods, len(cvProbsPeriods), nStates, zwl_LE, prevBeta, statProbs)
        result_betas.append(ubPLonger.copy())
        zwl_LE = findMinProbs(cvProbsPeriods, nStates, statProbs, ubPLonger)
        lp_UE = longerProbs_UE(cvProbsPeriods, len(cvProbsPeriods), nStates, statProbs, zwl_LE)
        zwl_UE = findMaxProbs(cvProbsPeriods, nStates, statProbs, ubPLonger)
        dmp = calcDMPHi(cvProbsPeriods, len(cvProbsPeriods), nStates, zwl_UE, ubPLonger, statProbs, nQs, deadline)
        result_zwl_lo.append(zwl_LE.copy())
        result_zwl_hi.append(zwl_UE.copy())
        result_dmp_hi.append(dmp)
        ue_started_decrease = np.logical_or(zwl_UE < prev_zwl_UE - eps_vec, ue_started_decrease)
        ue_stopped_decrease = np.logical_and(np.logical_or(zwl_UE > prev_zwl_UE - eps_vec, ue_stopped_decrease), ue_started_decrease)
        le_started_increase = np.logical_or(zwl_LE > prev_zwl_LE + eps_vec, le_started_increase)
        le_stopped_increase = np.logical_and(np.logical_or(zwl_LE < prev_zwl_LE + eps_vec, le_stopped_increase), le_started_increase)
        
    return (cvProbsPeriods, zwl_LE, zwl_UE, ubPLonger, result_betas, result_zwl_lo, result_zwl_hi, result_dmp_hi)




