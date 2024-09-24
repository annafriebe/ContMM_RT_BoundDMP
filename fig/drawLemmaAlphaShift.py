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
    
z = np.linspace(-1, 5, num=200)

mu = 1
sigma = 1
alpha_1 = 0
alpha_2 = 1
cdf_alpha1 = norm.cdf((alpha_1-mu)/sigma)
cdf_alpha2 = norm.cdf((alpha_2-mu)/sigma)
scaleK1 = 1/norm.cdf((mu-alpha_1)/sigma)
scaleK2 = 1/norm.cdf((mu-alpha_2)/sigma)

res1 = scaledCDF(z, mu, sigma, scaleK1, alpha_1, cdf_alpha1)
res2 = scaledCDF(z, mu, sigma, scaleK2, alpha_2, cdf_alpha2)
plt.plot(z, res1)
plt.plot(z, res2)
#plt.axvline(res[4], color="red")
plt.show()
filename = 'lemma_alpha_shift.csv'
with open(filename, 'w', newline='') as f:
    writer = csv.writer(f, delimiter=',')
    for i in range(len(z)):
        writer.writerow((z[i], res1[i], res2[i]))


        
        
