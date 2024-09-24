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



def scaledCDF(z, mu, sigma, K, alpha, cdf_alpha):
    cdf = norm.cdf(z, loc=mu, scale=sigma)
    for i in range(len(z)):
        if z[i] < alpha:
            cdf[i] = 0
        else:
            cdf[i] -= cdf_alpha
            cdf[i] *= K
    return cdf
    
z = np.linspace(0, 6, num=200)

mu = 1
sigma_1 = 1
sigma_2 = 2
alpha = 1
cdf_alpha1 = norm.cdf((alpha-mu)/sigma_1)
cdf_alpha2 = norm.cdf((alpha-mu)/sigma_2)
scaleK1 = 1/norm.cdf((mu-alpha)/sigma_1)
scaleK2 = 1/norm.cdf((mu-alpha)/sigma_2)

res1 = scaledCDF(z, mu, sigma_1, scaleK1, alpha, cdf_alpha1)
res2 = scaledCDF(z, mu, sigma_2, scaleK2, alpha, cdf_alpha2)
plt.plot(z, res1)
plt.plot(z, res2)
#plt.axvline(res[4], color="red")
plt.show()
filename = 'lemma_var_incr.csv'
with open(filename, 'w', newline='') as f:
    writer = csv.writer(f, delimiter=',')
    for i in range(len(z)):
        writer.writerow((z[i], res1[i], res2[i]))


        
        
