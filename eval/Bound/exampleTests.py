import numpy as np
from scipy.stats import norm

# Section 4.5
# 
p_in_2_zwl_1_mult = 0.875*0.1
p_in_2_zwl_2_mult = 0.125*0.3
p_co_2_zwl_1_mult_lo = p_in_2_zwl_1_mult * 0.5 
p_co_2_zwl_2_mult_lo = p_in_2_zwl_2_mult * 0.5

# norm.cdf(2/1)
p_co_2_zwl_1_mult_hi = p_in_2_zwl_1_mult * 0.5/norm.cdf(2) 
p_co_2_zwl_2_mult_hi = p_in_2_zwl_2_mult * 0.5/norm.cdf(2)
print(0.5/norm.cdf(2)) 
print("p co 2 mult 1 lo", p_co_2_zwl_1_mult_lo) 
print("p co 2 mult 2 lo", p_co_2_zwl_2_mult_lo) 
print("p co 2 mult 1 hi", p_co_2_zwl_1_mult_hi) 
print("p co 2 mult 2 hi", p_co_2_zwl_2_mult_hi) 

p_in_1_from_2_zwl1_mult_lo = p_co_2_zwl_1_mult_lo*0.7
p_in_1_from_2_zwl2_mult_lo = p_co_2_zwl_2_mult_lo*0.7
p_in_1_from_2_zwl1_mult_hi = p_co_2_zwl_1_mult_hi*0.7
p_in_1_from_2_zwl2_mult_hi = p_co_2_zwl_2_mult_hi*0.7

print("p in 1 from 2 mult 1 lo", p_in_1_from_2_zwl1_mult_lo) 
print("p in 1 from 2 mult 2 lo", p_in_1_from_2_zwl2_mult_lo) 
print("p in 1 from 2 mult 1 hi", p_in_1_from_2_zwl1_mult_hi) 
print("p in 1 from 2 mult 2 hi", p_in_1_from_2_zwl2_mult_hi) 

# Section 4.6

a = np.array([[0.7875, 0.0875],[0.0875, 0.0375]])

b = np.array([0.875, 0.125])

b1 = np.array([0.875-0.093, 0.125])
b2 = np.array([0.875, 0.125-0.026])


x = np.linalg.solve(a, b)
print(x)

x1 = np.linalg.solve(a, b1)
print(x1)

x2 = np.linalg.solve(a, b2)
print(x2)

b3 = np.array([0.875-0.093, 0.125-0.026])
x3 = np.linalg.solve(a, b3)
print(x3)

diff1 = x1 - x3
p1 = (1 - x3[1])/diff1[1]
#print(p1)
print(x3[0] + p1*diff1[0])

diff2 = x2 - x3
p2 = (1-x3[0])/diff2[0]
print(x3[1] + p2*diff2[1])

#print(0.875 - 0.88*0.7875 - 0.31*0.0875)
p_zwl_1_lo = x3[0] + p1*diff1[0]
p_zwl_2_lo = x3[1] + p2*diff2[1]

# Section 4.7

p_in_1_zwl_1_mult = 0.875*0.9
p_in_1_zwl_2_mult = 0.125*0.7
p_co_1_zwl_1_mult_lo = p_in_1_zwl_1_mult * norm.sf(2, loc=1, scale=0.5) 
p_co_1_zwl_2_mult_lo = p_in_1_zwl_2_mult * norm.sf(2, loc=1, scale=0.5)

# norm.cdf(1/0.5)=norm.cdf(2)
p_co_1_zwl_1_mult_hi = p_in_1_zwl_1_mult * norm.sf(2, loc=1, scale=0.5)/norm.cdf(2) 
p_co_1_zwl_2_mult_hi = p_in_1_zwl_2_mult * norm.sf(2, loc=1, scale=0.5)/norm.cdf(2)

print("p co 1 mult 1 lo", p_co_1_zwl_1_mult_lo) 
print("p co 1 mult 2 lo", p_co_1_zwl_2_mult_lo) 
#print("p co 1 mult 1 hi", p_co_1_zwl_1_mult_hi) 
#print("p co 1 mult 2 hi", p_co_1_zwl_2_mult_hi) 

p_in_1_from_1_zwl1_mult_lo = p_co_1_zwl_1_mult_lo*0.9
p_in_1_from_1_zwl2_mult_lo = p_co_1_zwl_2_mult_lo*0.9
p_in_1_from_1_zwl1_mult_hi = p_co_1_zwl_1_mult_hi*0.9
p_in_1_from_1_zwl2_mult_hi = p_co_1_zwl_2_mult_hi*0.9

print("p in 1 from 1 mult 1 lo", p_in_1_from_1_zwl1_mult_lo) 
print("p in 1 from 1 mult 2 lo", p_in_1_from_1_zwl2_mult_lo) 
#print("p in 1 from 1 mult 1 hi", p_in_1_from_1_zwl1_mult_hi) 
#print("p in 1 from 1 mult 2 hi", p_in_1_from_1_zwl2_mult_hi) 

p_in_sum_1_acc_1_lo = p_in_1_zwl_1_mult * p_zwl_1_lo + p_in_1_zwl_2_mult*p_zwl_2_lo
print("p_in lo sum state 1 acc period 1", p_in_sum_1_acc_1_lo)

p_in_sum_1_acc_2_lo = (p_in_1_from_1_zwl1_mult_lo + p_in_1_from_2_zwl1_mult_lo) * p_zwl_1_lo + \
(p_in_1_from_1_zwl2_mult_lo + p_in_1_from_2_zwl2_mult_lo) * p_zwl_2_lo
print("p_in lo sum state 1 acc period 2", p_in_sum_1_acc_2_lo)

beta_2_state_1 = 0.093 - p_in_sum_1_acc_2_lo
print("beta state 1 period 2 first eq", beta_2_state_1)

beta_2_state_1_sp_calc = 0.875 - p_in_sum_1_acc_1_lo - p_in_sum_1_acc_2_lo
print("beta state 1 period 2 second eq", beta_2_state_1_sp_calc)

# section 4.8
dmp_1_hi = norm.sf(4, loc=1, scale=0.5)/norm.cdf(2)
print("dmp 1 first", dmp_1_hi)
dmp_2_hi = norm.sf(4, loc=2, scale=1)/norm.cdf(2)
print("dmp 2 first", dmp_2_hi)

# zwl bound 1 for both
p_in_1_hi = p_in_1_zwl_1_mult + p_in_1_zwl_2_mult
print("p_in_1_hi", p_in_1_hi)
# zwl bound 1 for both
p_in_2_hi = p_in_2_zwl_1_mult + p_in_2_zwl_2_mult
print("p_in_2_hi", p_in_2_hi)

# /norm.cdf(0) = /0.5 = *2
K_inv_from_1 = norm.cdf(-2) # (1-2)/0.5
alpha_1_from_1 = norm.isf(K_inv_from_1, loc=(1+1-2), scale = np.sqrt(0.5))
print("alpha 1 from 1", alpha_1_from_1)
dmp_1_from_1_hi = norm.sf(4, loc = (1+1-2), scale = np.sqrt(0.25 + 0.25))/norm.cdf((1+1-2-alpha_1_from_1)/np.sqrt(0.5))
print("dmp 1 from 1", dmp_1_from_1_hi)
p_in_1_from_1_hi = p_in_1_from_1_zwl1_mult_hi + p_in_1_from_1_zwl2_mult_hi
print("p_in_1_from_1_hi", p_in_1_from_1_hi)

K_inv_from_2 = 0.5 # norm.cdf(0) # (2-2)/1
alpha_1_from_2 = norm.isf(K_inv_from_2, loc=(1+2-2), scale = np.sqrt(1.25))
print("alpha 1 from 2", alpha_1_from_2)
dmp_1_from_2_hi = norm.sf(4, loc = (1+2-2), scale = np.sqrt(0.25 + 1))/norm.cdf((1+2-2-alpha_1_from_2)/np.sqrt(1.25))
print("dmp 1 from 2", dmp_1_from_2_hi)
p_in_1_from_2_hi = p_in_1_from_2_zwl1_mult_hi + p_in_1_from_2_zwl2_mult_hi
print("p_in_1_from_2_hi", p_in_1_from_2_hi)

alpha_2_from_1 = norm.isf(K_inv_from_1, loc=(2+1-2), scale = np.sqrt(1.25))
print("alpha 2 from 1", alpha_2_from_1)
dmp_2_from_1 = norm.sf(4, loc=(1+2-2), scale = np.sqrt(0.25 + 1))/norm.cdf((1+2-2-alpha_2_from_1)/np.sqrt(1.25))
print("dmp 2 from 1", dmp_2_from_1)
p_in_2_from_1_zwl1_mult_lo = p_co_1_zwl_1_mult_lo*0.1
p_in_2_from_1_zwl2_mult_lo = p_co_1_zwl_2_mult_lo*0.1
p_in_2_from_1_zwl1_mult_hi = p_co_1_zwl_1_mult_hi*0.1
p_in_2_from_1_zwl2_mult_hi = p_co_1_zwl_2_mult_hi*0.1
p_in_2_from_1_hi = p_in_2_from_1_zwl1_mult_hi + p_in_2_from_1_zwl2_mult_hi
print("p_in_2_from_1_hi", p_in_1_from_2_hi)

alpha_2_from_2 = norm.isf(K_inv_from_2, loc=(2+2-2), scale = np.sqrt(2))
print("alpha 2 from 2", alpha_2_from_2)
dmp_2_from_2 = norm.sf(4, loc=(2+2-2), scale = np.sqrt(2))/norm.cdf((2+2-2-alpha_2_from_2)/np.sqrt(2))
print("dmp 2 from 2", dmp_2_from_2)
p_in_2_from_2_zwl1_mult_lo = p_co_2_zwl_1_mult_lo*0.3
p_in_2_from_2_zwl2_mult_lo = p_co_2_zwl_2_mult_lo*0.3
p_in_2_from_2_zwl1_mult_hi = p_co_2_zwl_1_mult_hi*0.3
p_in_2_from_2_zwl2_mult_hi = p_co_2_zwl_2_mult_hi*0.3
p_in_2_from_2_hi = p_in_2_from_2_zwl1_mult_hi + p_in_2_from_2_zwl2_mult_hi
print("p_in_2_from_2_hi", p_in_2_from_2_hi)

p_in_sum_2_acc_1_lo = p_in_2_zwl_1_mult * p_zwl_2_lo + p_in_2_zwl_2_mult*p_zwl_2_lo
print("p_in lo sum state 2 acc period 1", p_in_sum_2_acc_1_lo)

p_in_sum_2_acc_2_lo = (p_in_2_from_1_zwl1_mult_lo + p_in_2_from_2_zwl1_mult_lo) * p_zwl_1_lo + \
(p_in_2_from_1_zwl2_mult_lo + p_in_2_from_2_zwl2_mult_lo) * p_zwl_2_lo
print("p_in lo sum state 2 acc period 2", p_in_sum_2_acc_2_lo)

beta_2_state_2 = 0.026 - p_in_sum_2_acc_2_lo
print("beta state 2 period 2 first eq", beta_2_state_2)

beta_2_state_2_sp_calc = 0.125 - p_in_sum_2_acc_1_lo - p_in_sum_2_acc_2_lo
print("beta state 2 period 2 second eq", beta_2_state_2_sp_calc)

dmp_1_hi_beta_part = beta_2_state_1/0.875
dmp_1_hi_sum_part = (p_in_1_hi*dmp_1_hi + p_in_1_from_1_hi*dmp_1_from_1_hi + p_in_1_from_2_hi*dmp_1_from_2_hi)/0.875
print("dmp 1 beta part", dmp_1_hi_beta_part)
print("dmp 1 sum part", dmp_1_hi_sum_part)
print("dmp 1", dmp_1_hi_beta_part + dmp_1_hi_sum_part)

dmp_2_hi_beta_part = beta_2_state_2/0.125
dmp_2_hi_sum_part = (p_in_2_hi*dmp_2_hi + p_in_2_from_1_hi*dmp_2_from_1 + p_in_2_from_2_hi*dmp_2_from_2)/0.125
print("dmp 2 beta part", dmp_2_hi_beta_part)
print("dmp 2 sum part", dmp_2_hi_sum_part)


dmp_beta_part = beta_2_state_1 + beta_2_state_2
dmp_sum_part = p_in_1_hi*dmp_1_hi + p_in_1_from_1_hi*dmp_1_from_1_hi + p_in_1_from_2_hi*dmp_1_from_2_hi + \
    p_in_2_hi*dmp_2_hi + p_in_2_from_1_hi*dmp_2_from_1 + p_in_2_from_2_hi*dmp_2_from_2
print("dmp beta part", dmp_beta_part)
print("dmp sum part", dmp_sum_part)
print("dmp", dmp_sum_part + dmp_beta_part)

