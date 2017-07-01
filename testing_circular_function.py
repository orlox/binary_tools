#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 20:47:47 2017

@author: kaliroe
"""

import random as rd
import numpy as np
import matplotlib.pyplot as plt
import constants as cn
from scipy.stats import maxwell
rd.seed(10)
from functions import rand_phi, rand_theta, rand_velocity
from circular_function import circular_conditions

test_c = 100
test_M1 = 5.5
test_M2 = 55
test_Ai = 133
test_Mns = 1.4

"""creating an array of eccentricity and semi-major axis values""" 
  
testing_function = np.zeros([100,2])

for i in range(len(testing_function)):
    testing_function[i][0] = circular_conditions(test_Ai, test_M1, test_M2, test_Mns, rand_theta(),rand_phi(),rand_velocity(test_c))[1]
    testing_function[i][1] = circular_conditions(test_Ai, test_M1, test_M2, test_Mns, rand_theta(),rand_phi(),rand_velocity(test_c))[2]

"""changing the semi-major axis to period values in days"""

for k in range(len(testing_function)):
    testing_function[k][1] = (np.sqrt(((testing_function[k][1])**3)*4*(np.pi**2)/(cn.cgrav*(test_M1+test_Mns)*cn.Msun)))/(3600*24)
    
plt.plot(testing_function[:,0], testing_function[:,1], "o")
plt.xlim([0.0726,0.0728])
    