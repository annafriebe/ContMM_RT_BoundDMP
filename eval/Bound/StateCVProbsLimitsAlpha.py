# -*- coding: utf-8 -*-
"""

@author: Anna Friebe
"""
import numpy as np
import math
from scipy.stats import norm

def tableStringSN(probs, nStates):
    ts = ""
    for s in range(nStates):
        ts += " & "
        ts += "{:.2}".format(probs[s])
    return ts

class StateCVProbs:
    def __init__(self, state, cv, nStates):
        self._state = state
        self._cv = cv
        self._probInLo = np.zeros(nStates)
        self._probInHi = np.zeros(nStates)
        self._probOutLo = np.zeros(nStates)
        self._probOutHi = np.zeros(nStates)
        self._probOutLoDistr = 0
        self._probOutHiDistr = 0
        self._nStates = nStates
        self._mean = 0
        self._var = 0
        self._alpha = 0
        self._prevAlphaMax = 0
        self._kInv = 1
    def addInputProbs(self, transitionMatrix, prevState, prevOutputProbsLo, \
                      prevOutputProbsHi, prevAlpha, nQs):
        self._prevAlphaMax = max(self._prevAlphaMax, prevAlpha - nQs)
        for s in range(self._nStates):
            self._probInLo[s] += transitionMatrix[prevState, self._state]*\
                prevOutputProbsLo[s]
            self._probInHi[s] += transitionMatrix[prevState, self._state]*\
                prevOutputProbsHi[s]
        return
    def getInputProbsLo(self):
        return self._probInLo
    def getInputProbsHi(self):
        return self._probInHi
    def calcMeanVarOutputProbs(self, means, vars, nQs):
        mean = 0
        var = 0
        n = 0
        for s in range(self._nStates):
            mean += means[s]*self._cv[s]
            var += vars[s]*self._cv[s]
            n += self._cv[s]
        mean -= n*nQs
        if n > 1:
            meanPrevWL = mean - means[self._state] + nQs
            varPrevWL = var - vars[self._state]
            self._kInv = norm.cdf((meanPrevWL - self._prevAlphaMax)/math.sqrt(varPrevWL))
            self._alpha = norm.isf(self._kInv, loc=mean, scale=math.sqrt(var))
        self._probOutDistrHi = 1
        if self._alpha < 1e-10:
            self._probOutDistrHi =  norm.sf(0, mean, math.sqrt(var))/self._kInv
        self._probOutDistrLo = norm.sf(0, mean, math.sqrt(var))
        self._mean = mean
        self._var = var
        self._probOutHi = self._probInHi.copy()
        self._probOutHi *= self._probOutDistrHi
        self._probOutLo = self._probInLo.copy()
        self._probOutLo *= self._probOutDistrLo
        return
    def getOutputProbsHi(self):
        return self._probOutHi.copy()
    def getOutputProbsLo(self):
        return self._probOutLo.copy()
    def getProbDistrHi(self):
        return self._probOutDistrHi.copy()
    def getProbDistrLo(self):
        return self._probOutDistrLo.copy()
    def getMean(self):
        return self._mean
    def getVar(self):
        return self._var            
    def getAlpha(self):
        return self._alpha
    def getK(self):
        return 1/self._kInv
    def outputForTable(self):
        print(str(self._cv) + tableStringSN(self._probInLo, self._nStates) + \
              tableStringSN(self._probInHi, self._nStates) + \
              tableStringSN(self._probOutLo, self._nStates) + \
              tableStringSN(self._probOutHi, self._nStates) + " \\\\")
        print("\hline")
    
