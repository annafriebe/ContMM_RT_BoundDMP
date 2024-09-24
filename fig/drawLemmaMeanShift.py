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



def scaledCDF(z, mu, sigma, K, delta_alpha, cdf_alpha):
    cdf = norm.cdf(z, loc=mu, scale=sigma)
    for i in range(len(z)):
        if z[i] < mu + delta_alpha:
            cdf[i] = 0
        else:
            cdf[i] -= cdf_alpha
            cdf[i] *= K
    return cdf
    
z = np.linspace(-1, 5, num=200)

mu_1 = 1
mu_2 = 2
sigma = 1
delta_alpha = -1
cdf_alpha = norm.cdf(delta_alpha/sigma)
scaleK = 1/norm.cdf(-delta_alpha/sigma)

res1 = scaledCDF(z, mu_1, sigma, scaleK, delta_alpha, cdf_alpha)
res2 = scaledCDF(z, mu_2, sigma, scaleK, delta_alpha, cdf_alpha)
plt.plot(z, res1)
plt.plot(z, res2)
#plt.axvline(res[4], color="red")
plt.show()
filename = 'lemma_merge_shift.csv'
with open(filename, 'w', newline='') as f:
    writer = csv.writer(f, delimiter=',')
    for i in range(len(z)):
        writer.writerow((z[i], res1[i], res2[i]))


        
        
