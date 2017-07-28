
from binary_tools.utils.constants import *
import numpy as np

__author__ = "Kaliroe Pappas"
__credits__ = ["Kaliroe Pappas"]
__license__ = "GPL"
__version__ = "3.0"
__maintainer__ = "Kaliroe Pappas"
__email__ = "kmpappa2@illinois.edu"

"""various funcitons used in the test_kicks file"""

def angular_momentum(A,M1,M2,e):
    """A function that calculates the  angular 
    momentum of a two-body system in units of 
    grams*cm^2/s
    Arguments:    
        - A: the semi-major axis of the system
        - M1: solar mass of the first mass
        - M2: solar mass of the second mass
        - e: eccentricity of the orbit
    Returns the angular momentum of the system
    """
    L =  M1*Msun*M2*Msun*np.sqrt(cgrav*A*Rsun*(1-e**2)/((M1+M2)*Msun))
    return L


def separation(A,e,true_anomaly):
    """A function that calculates the separation
    of a binary system
    Arguments:
        - A: the semi-major axis of the system
        - e: eccentricity of the orbit
        - true_anomaly: the current true anomaly
    Returns the separation in solar radians.
    """
    separation = A*(1 - e**2)/(1+e*np.cos(true_anomaly))
    return separation