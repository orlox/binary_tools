#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 20:47:47 2017

@author: kaliroe
"""

import random as rd
import numpy as np
import matplotlib.pyplot as plt
from binary_tools.utils.constants import *
from binary_tools.utils import kicks
from binary_tools.utils.Keplers_laws import *
from scipy.stats import maxwell

__author__ = "Kaliroe Pappas"
__credits__ = ["Kaliroe Pappas"]
__license__ = "GPL"
__version__ = "3.0"
__maintainer__ = "Kaliroe Pappas"
__email__ = "kmpappa2@illinois.edu"


"""A function that tests the post_explosion_params_circular function""" 

def testing_circular_function(test_sigma = 15, test_M1 = 5.5, test_M2 = 55, test_Ai = 133, test_Mns = 1.4, seed="Flay",sample_velocity = 30):
    """Test that the graph of the eccentricity vs the period looks correct
    Arguments:
        - test_sigma: a sample sigma for the rand_velocity function
        - test_M1 = solar mass of the first mass pre-explosion
        - test_M2 = solar mass of the second mass
        - test_Ai = the initial maximum separation of the two masses
            pre-explosion
        - test_Mns = solar mass of the first mass post-explosion
        - seed: the seed used for the random number generator
        - sample_velocity = a constant velocity over which a line is
            drawn on the graph
    Returns: 'Compare graph to paper'
    """  
    rd.seed(seed)
    testing_function = np.zeros([1000,2])
    constant_velocity = np.zeros([1000,2])

    for i in range(len(testing_function)):  
        testing_function[i][0] = kicks.post_explosion_params_circular(test_Ai, test_M1, test_M2, test_Mns, kicks.rand_theta(),kicks.rand_phi(),kicks.rand_velocity(test_sigma))[0]
        testing_function[i][1] = kicks.post_explosion_params_circular(test_Ai, test_M1, test_M2, test_Mns, kicks.rand_theta(),kicks.rand_phi(),kicks.rand_velocity(test_sigma))[1]
    
    for j in range(len(constant_velocity)):
        constant_velocity[j][0] = kicks.post_explosion_params_circular(test_Ai, test_M1, test_M2, test_Mns, kicks.rand_theta(),kicks.rand_phi(),sample_velocity)[0]
        constant_velocity[j][1] = kicks.post_explosion_params_circular(test_Ai, test_M1, test_M2, test_Mns, kicks.rand_theta(),kicks.rand_phi(),sample_velocity)[1]
    
    """changing the semi-major axis to period values in days"""
    
    for k in range(len(testing_function)):
        testing_function[k][0] = keplers_third_law(testing_function[k][0],test_M2,test_Mns)
        constant_velocity[k][0] = keplers_third_law(constant_velocity[k][0],test_M2,test_Mns)
    
    plt.plot(testing_function[:,0], testing_function[:,1], "o")
    plt.plot(constant_velocity[:,0], constant_velocity[:,1], 'k-')
    plt.show()
    plt.close()
    
    rd.seed()
    
    return "Compare graph to paper"
    