#!/usr/bin/env python
from binary_tools.utils.constants import *
from binary_tools.utils import kicks
import matplotlib.pyplot as plt
from scipy.stats import maxwell
import random as rd
import numpy as np

__author__ = "Kaliroe Pappas"
__credits__ = ["Kaliroe Pappas"]
__license__ = "GPL"
__version__ = "3.0"
__maintainer__ = "Kaliroe Pappas"
__email__ = "kmpappa2@illinois.edu"

"""Test the functions in the kicks.py module
"""

def test_rand_phi(num_sample=100000, nbins=20, tolerance = 1e-3, seed="Jean", plot=False):
    """Test that phi is sampled as flat from zero to pi
    Arguments:
        - num_sample: number of random phi generated
        - nbins: random sampled numbers will be binned to compute
            probabilities and compare to expected values. This variable
            specifies the number of bins used.
        - tolerance: tolerance for the test
        - seed: the seed used for the random number generator
        - plot: if true, plot results
    Returns: True if the test is succesful, False otherwise
    """
    rd.seed(seed)
    phi_array = np.zeros(num_sample)
    for i in range(0,len(phi_array)):
        phi_array[i] = kicks.rand_phi()

    #do a histogram
    #TODO: use numpy histogram, to avoid plotting if its not neccesary
    vals_phi, bins_phi = np.histogram(phi_array, bins=np.linspace(0,np.pi,nbins))
    if plot:
        plt.hist(phi_array, bins=np.linspace(0,np.pi,20))
        plt.title("phi distribution")
        plt.xlabel("phi value")
        plt.ylabel("distribution")
        plt.show()
        plt.close()

    #check if the probability computed for each bin is within the tolerance
    success = True
    tolerance = max(vals_phi)/10**3
    for k in range(0,len(vals_phi)):
        prob_hist = vals_phi[k]/(sum(vals_phi))
        prob_test = (1/(2*np.pi))*(bins_phi[k+1]-bins_phi[k])
        if abs(prob_test-prob_hist)>tolerance:
            success = False
            break

    #re-seed the random number generator
    rd.seed()

    return success

def test_rand_theta(num_sample=100000, nbins=20, tolerance = 1e-3, seed="Jubilee", plot=False):
    """Test that theta is sampled as a sign graph from zero to pi
    Arguments:
        - num_sample: number of random theta generated
        - nbins: random sampled numbers will be binned to compute
            probabilities and compare to expected values. This variable
            specifies the number of bins used.
        - tolerance: tolerance for the test
        - seed: the seed used for the random number generator
        - plot: if true, plot results
    Returns: True if the test is succesful, False otherwise
    """
    rd.seed(seed)
    theta_array = np.zeros(num_sample)
    for i in range(0,len(theta_array)):
        theta_array[i] = kicks.rand_theta()

    #do a histogram
    #TODO: use numpy histogram, to avoid plotting if its not neccesary
    vals_theta, bins_theta = np.histogram(theta_array, bins=np.linspace(0,np.pi,nbins))
    if plot:
        plt.hist(theta_array, bins=np.linspace(0,np.pi,20))
        plt.title("theta distribution")
        plt.xlabel("theta value")
        plt.ylabel("distribution")
        plt.show()
        plt.close()
    #check if the probability computed for each bin is within the tolerance
    success = True
    tolerance = max(vals_theta)/10**3
    for k in range(0,len(vals_theta)):
        prob_hist = vals_theta[k]/(sum(vals_theta))
        prob_test = -(np.cos(bins_theta[k+1])-np.cos(bins_theta[k]))/2
        if abs(prob_test-prob_hist)>tolerance:
            success = False
            break

    #re-seed the random number generator
    rd.seed()

    return success

def test_rand_velocity(sigma, num_sample=10000, nbins=20, tolerance=1e-3, seed="Dimitris", plot=False):
    """Test that the velocity output is sampled as a maxwellian
    Arguments:
        - sigma: argument needed to run rand_velocity
        - num_sample: number of random theta generated
        - nbins: random sampled numbers will be binned to compute
            probabilities and compare to expected values. This variable
            specifies the number of bins used.
        - tolerance: tolerance for the test
        - seed: the seed used for the random number generator
        - plot: if true, plot results
    Returns: True if the test is succesful, False otherwise
    """
    rd.seed(seed)
    velocity_array = np.zeros(num_sample)
    for k in range(0,len(velocity_array)):
        velocity_array[k] = kicks.rand_velocity(sigma)
    
    #do a histogram
    vals_velocity, bins_velocity = np.histogram(velocity_array, bins=np.linspace(0,5*sigma,20))
    if plot:    
        plt.hist(velocity_array, bins=np.linspace(0,5*sigma,20))
        plt.title("velocity distribution")
        plt.xlabel("velocity value")
        plt.ylabel("distribution")
        plt.show()
        plt.close()
    #check if the probability computed from each bin is within the tolerance
    success = True
    tolerance = max(vals_velocity)/10**3
    for k in range(0,len(vals_velocity)):
        prob_hist = vals_velocity[k]/(sum(vals_velocity))
        prob_test = maxwell.cdf(bins_velocity[k+1],0,sigma) - maxwell.cdf(bins_velocity[k+1],0,sigma)
        if abs(prob_test-prob_hist)>tolerance:
            success = False
            break
    
    #re-seed the random number generator
    rd.seed()
    
    return success

        
        