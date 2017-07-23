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
    Calcuted from the following distribution:
        f(x) = (1/2)sin(theta)
        
    Returns: the randomly sampled value of theta
    """
    rand_theta = np.arccos(1 - 2*rd.random())
    return rand_theta

def rand_phi():
    """Randomly sample the azimuthal angle of the kick direction.
    Calcuted from the following distribution:
        f(x) = 1
    
    Returns: the randomly sampled value of phi
    """
    rand_phi = 2*np.pi*rd.random()
    return rand_phi

def rand_velocity(sigma):
    """Randomly sample the velocity of the kick direction.
    The kick is assumed to follow a gaussian distribution given
    by sigma in each cartesian component of the kick, which makes
    up the full maxwellian using the velocity in three dimensions:
    f(v) = (2*pi*sigma^2)^(-1/2)*exp(-v^2/(2*sigma^2))
    
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
    the angle between the pre and post explosion orbital plane, and 
    a boolean describing if the orbit is unbound.
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
    
    bound = True
    if e > 1:
        bound = False
        
    

    return Af/Rsun, e, theta_new, bound

def rand_true_anomaly(e):
    """Computes a random true anomaly an eliptical orbit using a
     Montecarlo sampling method and drawing a random number from 
     an array of possibilities. Calculation follows Dosopoulou
     (2016), ApJ, 825, 70D
     Calcuted from the following distribution:
        f(x) = np.sqrt(1/(4*np.pi))*(1-e**2)**(1.5)/((1+e*np.cos(u))**2)
        
    Arguments:
        - e: eccentricity of the orbit
        -the maximum value was chosen based on the principle that
            a mass in an eliptical orbit moves the slowest at apogee
            and therefore is most likely to be found there. In this 
            case apogee occurs at pi.
    Returns: a random true_anomaly
    """
    
    max = np.sqrt(1/(4*np.pi))*(1-e**2)**(1.5)/((1+e*np.cos(np.pi))**2)
    
    while True:
        randx = rd.uniform(0,2*np.pi)
        randy = rd.uniform(0, max)
        
        if np.sqrt(1/(4*np.pi))*(1-e**2)**(1.5)/((1+e*np.cos(randx))**2)>randy:
            return randx
            
    



def rand_separation(e, Ai):
    """Computes a random separation using the rand_true_anomaly
    function.
    Calcuted from the following distribution:
        f(x) = Ai*(1 - e**2)/(1+e*np.cos(u))
    Arguments:
        - e: eccentricity of the orbit
        - Ai: initial orbital separation in Rsun
    Returns: a random separation in Rsun
    """
    u = rand_true_anomaly(e)
    separation = Ai*(1 - e**2)/(1+e*np.cos(u))
    
    return separation


def post_explosion_params_general( Ai, M1, M2, Mns, e, theta, phi, Vk, true_anomaly):
    """Computes the post explosion orbital parameters of a binary system,
    assuming the pre-explosion system is in a circular orbit.
    Calculation follows Kalogera (1996), ApJ, 471, 352

    Arguments:
        - Ai: initial orbital separation in Rsun
        - M1: initial mass of the exploding star in Msun
        - M2: mass of the companion star in Msun
        - Mns: final mass of the exploding star in Msun
        - e: eccentricity of the orbit
        - theta: polar angle of the kick direction
        - phi: azimuthal angle of the kick direction
        - Vk: kick velocity, in km/s
        - true_anomaly: the true anomaly in radians

    Returns: The final separation in Rsun, final eccentricity, and 
    a boolean describing if the orbit is unbound.
    """

    #turn input into CGS
    M1 = M1*Msun
    M2 = M2*Msun
    Mns = Mns*Msun
    Ai = Ai*Rsun
    Vk = Vk*1e5
    separation = Ai*(1 - e**2)/(1+e*np.cos(true_anomaly))
    
    V_theta = np.sqrt(cgrav*(M1+M2)*Ai*(1-e**2))/separation
    V_radius = np.sqrt(cgrav*(M1+M2)*(2/separation-1/Ai-Ai*(1-e**2)/(separation**2)))
    
    V_squared = (Vk**2)*np.cos(theta)**2 + 2*Vk*np.cos(theta)*V_theta +\
    V_theta**2 + V_radius**2 + (Vk**2)*np.sin(theta)**2 +\
    2*Vk*np.sin(theta)*np.cos(phi)*V_radius
    
    Af = cgrav*(Mns + M2)*(2*cgrav*(Mns + M2)/separation - V_squared)**(-1)
    
    e_final = np.sqrt(1 - ((Vk*np.cos(theta))**2 + 2*Vk*V_theta*np.cos(theta)\
    + V_theta**2 + (Vk*np.sin(theta)*np.sin(phi))**2)*separation**2/(cgrav*(Mns + M2)*Af))
    
    bound = True
    if e_final > 1:
        bound = False
    

    return Af/Rsun, e_final, bound
        

