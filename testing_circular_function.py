#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 20:47:47 2017

@author: kaliroe
"""

import random as rd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import maxwell
rd.seed(9)

from functions import rand_phi, rand_theta, rand_velocity

G = 6.67408*10**(-11)  #m**3 kg**-1 s**-2 

from circular_function import circular_conditions

Mo = 1.989*10**30 #kg
