#!/usr/bin/env python
from binary_tools.utils.constants import *
from binary_tools.utils import kicks
import matplotlib.pyplot as plt
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

def test_rand_phi(num_sample=10000, nbins=20, tolerance = 1e-3, seed="Jean", plot=False):
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
    vals_phi, bins_phi, patches_phi = plt.hist(phi_array, bins=np.linspace(0,np.pi,nbins))
    if plot:
        plt.title("phi distibution")
        plt.xlabel("phi value")
        plt.ylabel("distribution")
        plt.show()

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
