# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 15:43:33 2021

@author: afe02
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
import math
import csv



def calcMean(z, mu_t, mu_wl, var_t, var_wl):
#    return ((z - mu_t)*var_t + mu_wl*var_wl)/(var_t + var_wl)
    return ((z - mu_t)*var_wl + mu_wl*var_t)/(var_t + var_wl)
    
def calcVar(var_t, var_wl):
    return var_t*var_wl/(var_t+var_wl)

def scalingFactor(mu_wl, var_wl):
    return 1/norm.cdf(mu_wl/math.sqrt(var_wl))

def calcFactor(z, mu_t, mu_wl, var_t, var_wl):
    mean = calcMean(z, mu_t, mu_wl, var_t, var_wl)
    var = calcVar(var_t, var_wl)
    res = norm.sf(0, loc=mean, scale = math.sqrt(var))
    return res
    
def scaledNorm(z, factors, mu_t, mu_wl, var_t, var_wl):
    mean = mu_t + mu_wl
    var = var_t + var_wl
    return norm.pdf(z, loc=mean, scale=math.sqrt(var))*factors

def normZ(z, mu_t, mu_wl, var_t, var_wl):
    mean = mu_t + mu_wl
    var = var_t + var_wl
    return norm.pdf(z, loc=mean, scale=math.sqrt(var))

def scaledMeanFactor(z, mu_t, mu_wl, var_t, var_wl):
    mean = mu_t + mu_wl
    factorMean = calcFactor(mean, mu_t, mu_wl, var_t, var_wl)
    print(factorMean)
    var = var_t + var_wl
    return factorMean*norm.pdf(z, loc=mean, scale=math.sqrt(var))

def shiftNorm(z, mu_t, mu_wl, var_t, var_wl):
    mean = mu_t + mu_wl
    factorMean = calcFactor(mean, mu_t, mu_wl, var_t, var_wl)
    var = var_t + var_wl
    shift = math.sqrt(-2*var*math.log(factorMean))
    mean = mean + shift
    return norm.pdf(z, loc=mean, scale=math.sqrt(var))

def findAlpha(mu_t, mu_wl, var_t, var_wl):
    KInv = norm.cdf(mu_wl/math.sqrt(var_wl))
    print(KInv)
    alpha = norm.isf(KInv, loc = mu_t + mu_wl, scale = math.sqrt(var_t + var_wl))
    print(alpha)
    return alpha, 1/KInv

def eval(mu_t, mu_wl, var_t, var_wl, z):
    factor = calcFactor(z, mu_t, mu_wl, var_t, var_wl)
#    sca = scalingFactor(mu_wl, var_wl)
    convResult = scaledNorm(z, factor, mu_t, mu_wl, var_t, var_wl)
    nz = normZ(z, mu_t, mu_wl, var_t, var_wl)
    sn = shiftNorm(z, mu_t, mu_wl, var_t, var_wl)
    alpha, K = findAlpha(mu_t, mu_wl, var_t, var_wl)
    return [factor, convResult, nz, sn, alpha, K]

def findHalf(factor, z):
    for i in range(len(factor)-1):
        if (factor[i] < 0.5) and (factor[i+1]>=0.5):
            return (z[i]+z[i+1])/2
    return 0

z = np.linspace(0, 10, num=400)

res = eval(1, 0, 0.25, 1, z)
plt.plot(z, res[1]*res[5])
plt.plot(z, res[2]*res[5])
plt.plot(z, res[2])
#plt.plot(z, res[0]/10)
#half = findHalf(res[0], z)
#print(half)
#plt.axvline(half, color="blue")
print(res[4])
plt.axvline(res[4], color="red")

#plt.axvline(21.5, color="red")
plt.show()
filename = 'ub_conv_ill_norm_simple_example.csv'
with open(filename, 'w', newline='') as f:
    writer = csv.writer(f, delimiter=',')
    for i in range(len(z)):
        partGaussVal =  res[2][i]*res[5]
        if z[i] < res[4]:
            partGaussVal = 0
        writer.writerow((z[i], res[0][i], res[1][i]*res[5], partGaussVal, res[2][i]))


        
        
