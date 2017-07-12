#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 12:00:19 2017

@author: kaliroe
"""

import numpy as np
from binary_tools.utils.constants import *

__author__ = "Kaliroe Pappas"
__credits__ = ["Kaliroe Pappas"]
__license__ = "GPL"
__version__ = "3.0"
__maintainer__ = "Kaliroe Pappas"
__email__ = "kmpappa2@illinois.edu"

"""Kepler's laws in function form"""

def keplers_third_law(semi_major_axis, mass_1, mass_2):
    """Takes accepted arguments and returns the period
    of the two-mass system
    Arguments:
        semi_major_axis: largest seperation of the masses in units of solar radii
        mass_1: mass of the first mass in solar masses
        mass_2: mass of the second mass in solar masses
    Returns the period of the system in units of days"""
    
    mass_1 = mass_1*Msun
    mass_2 = mass_2*Msun
    semi_major_axis = semi_major_axis*Rsun
    period = np.sqrt(4*np.pi**2*(semi_major_axis**3)/(cgrav*(mass_1+mass_2)))/(3600*24)
    return period
    