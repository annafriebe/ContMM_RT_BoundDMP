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





def normZ(z, mu_t, mu_wl, var_t, var_wl):
    mean = mu_t + mu_wl
    var = var_t + var_wl
    return norm.pdf(z, loc=mean, scale=math.sqrt(var))

def calcKInvAlphaDelta(mu_h, var_h, nQ, alpha_delta, mu_s, var_s):
    KInv = norm.cdf((mu_h - nQ - alpha_delta)/math.sqrt(var_h))
    alpha = norm.isf(KInv, loc = mu_h + mu_s - nQ, scale = math.sqrt(var_h + var_s))
    return KInv, alpha

def eval_alphas(mu_h, var_h, nQ, alpha_delta_1, alpha_delta_2, mu_s, var_s, z):
    K1Inv, alpha1 = calcKInvAlphaDelta(mu_h, var_h, nQ, alpha_delta_1, mu_s, var_s)
    K2Inv, alpha2 = calcKInvAlphaDelta(mu_h, var_h, nQ, alpha_delta_2, mu_s, var_s)
    nz = normZ(z, mu_s, mu_h-nQ, var_s, var_h)
    return (nz, K1Inv, alpha1, K2Inv, alpha2)

def tail_Ks(mu_h, var_h, nQ, alpha_delta_1, alpha_delta_2, z):
    KInv1 = norm.cdf((mu_h - nQ - alpha_delta_1)/math.sqrt(var_h))
    KInv2 = norm.cdf((mu_h - nQ - alpha_delta_2)/math.sqrt(var_h))
    n_tail_1 = norm.pdf(z, loc=mu_h-nQ, scale=math.sqrt(var_h))/KInv1
    n_tail_2 = norm.pdf(z, loc=mu_h-nQ, scale=math.sqrt(var_h))/KInv2
    for i in range(len(z)):
        if z[i] < alpha_delta_1:
            n_tail_1[i] = 0
        if z[i] < alpha_delta_2:
            n_tail_2[i] = 0
    return (n_tail_1, n_tail_2)



# comparing
print(0.5/norm.cdf(2))

z = np.linspace(0, 4, num=400)

mu_h = 1  # 1 + 2 - 2
var_h = 1.25
mu_s = 2
var_s = 1
nQ = 2
alpha_delta_1 = 0
alpha_delta_2 = 1.236

norm_arr, K1Inv, alpha1, K2Inv, alpha2 = \
    eval_alphas(mu_h, var_h, nQ, alpha_delta_1, alpha_delta_2, mu_s, var_s, z)

tail_1, tail_2 = tail_Ks(mu_h, var_h, nQ, alpha_delta_1, alpha_delta_2, z)

plt.plot(z, tail_1)
plt.plot(z, tail_2)

plt.axvline(alpha_delta_1, color="blue")
plt.axvline(alpha_delta_2, color="red")


#plt.plot(z, norm/K1Inv)
#plt.plot(z, norm/K2Inv)
#plt.plot(z, res[0]/10)
#half = findHalf(res[0], z)
#print(alpha1)
#plt.axvline(alpha1, color="blue")
#print(alpha2)
#plt.axvline(alpha2, color="red")

#plt.axvline(21.5, color="red")
plt.show()
filename = 'simple_example_vector_combine.csv'
with open(filename, 'w', newline='') as f:
    writer = csv.writer(f, delimiter=',')
    for i in range(len(z)):
        partGaussVal1 =  norm_arr[i]/K1Inv
        partGaussVal2 =  norm_arr[i]/K2Inv
        if z[i] < alpha1:
            partGaussVal1 = 0
        if z[i] < alpha2:
            partGaussVal2 = 0
        writer.writerow((z[i], tail_1[i], tail_2[i], partGaussVal1, partGaussVal2))


        
        
