#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 11:13:47 2017

@author: kaliroe
"""
import random as rd
import numpy as np
import matplotlib.pyplot as plt
import constants as cn


def circular_conditions(Ai, M1, M2, Mns, theta, phi, Vk): #mass in solar masses, length in solar radians
    M1 = M1*cn.Msun
    M2 = M2*cn.Msun
    Mns = Mns*cn.Msun
    Ai = Ai*cn.Rsun
    Vr = np.sqrt(cn.cgrav*(M1+M2)/Ai)
    Vkz = np.cos(theta)*Vk
    Vky = np.sin(theta)*np.sin(phi)*Vk
    Vkx = np.sin(theta)*np.cos(phi)*Vk
    Af = cn.cgrav*(Mns + M2)/(2*cn.cgrav*(Mns+M2)/Ai -Vk**2 -Vr**2 - 2*Vky*Vr)
    e = np.sqrt(1 - (Vkz**2 + Vky**2 + Vr**2 + 2*Vky*Vr)*Ai**2/(cn.cgrav*(Mns+M2)*Af))
    theta_new = np.arccos((Vky + Vr)/np.sqrt((Vky+Vr)**2 +Vkz**2))
    return theta_new, e, Af



