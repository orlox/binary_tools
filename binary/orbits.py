from binary_tools.constants import *
import numpy as np

__author__ = "Kaliroe Pappas"
__credits__ = ["Kaliroe Pappas"]
__license__ = "GPL"
__version__ = "3.0"
__maintainer__ = "Kaliroe Pappas"
__email__ = "kmpappa2@illinois.edu"

"""Provides various functions to get information on a binary system"""

def orbital_angular_momentum(a,m1,m2,e):
    """Calculates the orbital angular 
    momentum of a two-body system in units of 
    grams*cm^2/s
    Arguments:    
        - a: the semi-major axis of the system in solar radii
        - m1: mass of one of the components in solar masses
        - m2: mass of the other component in solar masses
        - e: eccentricity of the orbit
    Returns the angular momentum of the system
    """
    L =  m1*Msun*m2*Msun*np.sqrt(cgrav*a*Rsun*(1-e**2)/((m1+m2)*Msun))
    return L


def binary_separation(a,e,true_anomaly):
    """Calculates the separation of a binary system given the
    true anomaly and the semimajor axis
    Arguments:
        - a: the semi-major axis of the system
        - e: eccentricity of the orbit
        - true_anomaly: the current true anomaly
    Returns the separation in the same units as the input semimajor axis.
    """
    separation = a*(1 - e**2)/(1+e*np.cos(true_anomaly))
    return separation

def kepler3_a(P,m1,m2):
    """Calculates the semimajor axis of a binary from its period and masses
    using Kepler's third law
    Arguments:
        - P: the orbital period in days
        - e: eccentricity of the orbit
        - true_anomaly: the current true anomaly
    Returns the semimajor axis in solar radii.
    """
    return float(((P*24.*3600.)**2*cgrav*(m1+m2)*Msun/(4.*np.pi**2))**(1./3.)/Rsun);

def kepler3_P(a,m1,m2):
    """Calculates the period of a binary from its semimajor axis and masses
    using Kepler's third law
    Arguments:
        - a: the semi-major axis of the system in solar radii
        - e: eccentricity of the orbit
        - true_anomaly: the current true anomaly
    Returns the semimajor axis in solar radii.
    """
    return float((4.*np.pi**2*(a*Rsun)**3/(cgrav*(m1+m2)*Msun))**(1./2.)/(24.*3600.));
