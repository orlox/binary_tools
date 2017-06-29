#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 13:09:52 2017

@author: kaliroe
"""

import random as rd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import maxwell
rd.seed(9)
from functions import rand_theta(), rand_phi(), rand_velocity(c)


phi_array = np.zeros(10000)
for i in range(0,len(phi_array)):
    phi_array[i] = rand_phi()
    
theta_array = np.zeros(10000)
for j in range(0,len(theta_array)):
    theta_array[j] = rand_theta()

test_c = 10

velocity_array = np.zeros(10000)
for k in range(0,len(velocity_array)):
    velocity_array[k] = rand_velocity(test_c)

'''histogram'''

vals_phi, bins_phi, patches_phi = plt.hist(phi_array, bins=np.linspace(0,np.pi,20))
plt.title("phi distibution")
plt.xlabel("phi value")
plt.ylabel("distribution")
plt.show()

for k in range(0,len(vals_phi)):
    tolerance = max(vals_phi)/10**3
    prob_hist = vals_phi[k]/(sum(vals_phi))
    prob_test = (1/(2*np.pi))*(bins_phi[k+1]-bins_phi[k])
    assert((prob_test-tolerance) < prob_hist < (prob_test +tolerance)), "phi histogram is not a good fit"

vals_theta, bins_theta, patches_theta = plt.hist(theta_array, bins=np.linspace(0,np.pi,20))
plt.title("theta distibution")
plt.xlabel("theta value")
plt.ylabel("distribution")
plt.show()

for k in range(0,len(vals_theta)):
    tolerance = max(vals_theta)/10**3
    prob_hist = vals_theta[k]/(sum(vals_theta))
    prob_test = -(np.cos(bins_theta[k+1])-np.cos(bins_theta[k]))/2
    assert((prob_test-tolerance) < prob_hist < (prob_test +tolerance)), "theta histogram is not a good fit"

vals_velocity, bins_velocity, patches_velocity = plt.hist(velocity_array, bins=np.linspace(0,5*test_c,20))
plt.title("velocity distibution")
plt.xlabel("velocity value")
plt.ylabel("distribution")
plt.show()

for k in range(0,len(vals_velocity)):
    tolerance = max(vals_velocity)/10**3
    prob_hist = vals_velocity[k]/(sum(vals_velocity))
    prob_test = maxwell.cdf(bins_velocity[k+1],0,test_c) - maxwell.cdf(bins_velocity[k+1],0,test_c)
    assert((prob_test-tolerance) < prob_hist < (prob_test +tolerance)), "velcity histogram is not a good fit"




