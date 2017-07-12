#!/usr/bin/env python
import random as rd
import numpy as np
from binary_tools.utils.constants import *
from scipy.stats import maxwell

__author__ = "Kaliroe Pappas"
__credits__ = ["Kaliroe Pappas"]
__license__ = "GPL"
__version__ = "3.0"
__maintainer__ = "Kaliroe Pappas"
__email__ = "kmpappa2@illinois.edu"

"""Variety of functions meant for Montecarlo sampling of kicks,
and their effect on orbital parameters.
"""

def rand_kick(sigma):
    """Randomly sample an isotropic kick with a 1D sigma equal
    to the provided one. This means that in cartesian coordinates
    each component of the velocity follows a Gaussian distribution,

    f(v_x) = (2*pi*sigma^2)^(-1/2)*exp(-v_x^2/(2*sigma^2))

    Returns: three floating points corresponding to two angles
        and a velocity
    """
    return rand_theta(), rand_phi(), rand_velocity(sigma)

def rand_theta():
    """Randomly sample the polar angle of the kick direction.
    We define theta=0 as a kick in the same direction of the orbital motion

    Returns: the randomly sampled value of theta
    """
    rand_theta = np.arccos(1 - 2*rd.random())
    return rand_theta

def rand_phi():
    """Randomly sample the azimuthal angle of the kick direction.

    Returns: the randomly sampled value of phi
    """
    rand_phi = 2*np.pi*rd.random()
    return rand_phi

def rand_velocity(sigma):
    """Randomly sample the velocity of the kick direction.
    The kick is assumed to follow a gaussian distribution given
    by sigma in each cartesian component of the kick. This results in
    the full kick velocity following a Maxwellian distribution.

    Returns: the randomly sampled value of v, in the units of sigma
    """
    rand_velocity = maxwell.isf(rd.random(), 0, scale = sigma)
    return rand_velocity

def post_explosion_params_circular(Ai, M1, M2, M1f, theta, phi, Vk):
    """Computes the post explosion orbital parameters of a binary system,
    assuming the pre-explosion system is in a circular orbit.
    Calculation follows Kalogera (1996), ApJ, 471, 352

    Arguments:
        - Ai: initial orbital separation in Rsun
        - M1: initial mass of the exploding star in Msun
        - M2: mass of the companion star in Msun
        - M1f: final mass of the exploding star in Msun
        - theta: polar angle of the kick direction
        - phi: azimuthal angle of the kick direction
        - Vk: kick velocity, in km/s

    Returns: The final separation in Rsun, final eccentricity,
        and the angle between the pre and post explosion orbital
        planes.
    """

    #turn input into CGS
    M1 = M1*Msun
    M2 = M2*Msun
    M1f = M1f*Msun
    Ai = Ai*Rsun
    Vk = Vk*1e5

    #compute the final relative velocity of the system
    #See Kalogera (1996) for a definition of the axes
    Vr = np.sqrt(cgrav*(M1+M2)/Ai)
    Vky = np.cos(theta)*Vk
    Vkz = np.sin(theta)*np.sin(phi)*Vk
    Vkx = np.sin(theta)*np.cos(phi)*Vk

    #compute post explosion orbital parameters
    # Eqs. (3), (4) and (5) of Kalogera (1996)
    Af = cgrav*(M1f + M2)/(2*cgrav*(M1f+M2)/Ai -Vk**2 -Vr**2 - 2*Vky*Vr)
    e = np.sqrt(1 - (Vkz**2 + Vky**2 + Vr**2 + 2*Vky*Vr)*Ai**2/(cgrav*(M1f+M2)*Af))
    theta_new = np.arccos((Vky + Vr)/np.sqrt((Vky+Vr)**2 +Vkz**2))
    
    unbound = False
    if e > 1:
        unbound = True
    

    return Af/Rsun, e, theta_new, unbound

