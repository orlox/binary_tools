#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 11:13:47 2017

@author: kaliroe
"""
import random as rd
import numpy as np
import matplotlib.pyplot as plt

G = 6.67408*10**(-11)  #m**3 kg**-1 s**-2

def circular_conditions(Ai, M1, M2, Mns, theta, phi, Vk):
    Vr = (G*(M1+M2)/Ai)**(1/2)
    Vkx = np.cos(theta)*Vk
    Vky = np.sin(theta)*np.sin(phi)*Vk
    Vkz = np.sin(theta)*np.cos(phi)*Vk
    Af = G*(Mns + M2)*(2*G*(Mns+M2)/Ai -Vk**2 -Vr**2 - 2*Vky*Vr)**(-1)
    e = (1 - (Vkz**2 + Vky**2 + Vr**2 + 2*Vky*Vr)*Af**2/(G*(Mns+M2)*Af))**(1/2)
    theta_new = np.arccos((Vky + Vr)*((Vky+Vr)**2 +Vkz**2)**(-1/2))
    return theta_new, e, Af



