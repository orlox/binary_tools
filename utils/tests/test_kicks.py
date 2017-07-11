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

def test_rand_phi(num_sample=10000, nbins=20, tolerance = 1e-3, seed="Jean", plot=False, save=True):
    """Test that phi is sampled as flat from zero to pi
    Arguments:
        - num_sample: number of random phi generated
        - nbins: random sampled numbers will be binned to compute
            probabilities and compare to expected values. This variable
            specifies the number of bins used.
        - tolerance: tolerance for the test
        - seed: the seed used for the random number generator
        - plot: if true, plot results
        - save: saves the plot
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
        plt.hist(phi_array, bins=np.linspace(0,np.pi,nbins))
        plt.title("phi distribution")
        plt.xlabel("phi value")
        plt.ylabel("distribution")
        plt.show()
        plt.close()
    if save:
        plt.hist(phi_array, bins=np.linspace(0,np.pi,nbins))
        plt.title("phi distribution")
        plt.xlabel("phi value")
        plt.ylabel("distribution")
        plt.savefig("phi_distribution.png")
        plt.close()

    #check if the probability computed for each bin is within the tolerance
    success = True
    tolerance = max(vals_phi)*tolerance
    for k in range(0,len(vals_phi)):
        prob_hist = vals_phi[k]/(sum(vals_phi))
        prob_test = (1/(2*np.pi))*(bins_phi[k+1]-bins_phi[k])
        if abs(prob_test-prob_hist)>tolerance:
            success = False
            break

    #re-seed the random number generator
    rd.seed()

    return success

def test_rand_theta(num_sample=10000, nbins=20, tolerance = 1e-3, seed="Jubilee", plot=False, save=True):
    """Test that theta is sampled as a sign graph from zero to pi
    Arguments:
        - num_sample: number of random theta generated
        - nbins: random sampled numbers will be binned to compute
            probabilities and compare to expected values. This variable
            specifies the number of bins used.
        - tolerance: tolerance for the test
        - seed: the seed used for the random number generator
        - plot: if true, plot results
        - save: saves the plot
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
        plt.hist(theta_array, bins=np.linspace(0,np.pi,nbins))
        plt.title("theta distribution")
        plt.xlabel("theta value")
        plt.ylabel("distribution")
        plt.show()
        plt.close()
    if save:
        plt.hist(theta_array, bins=np.linspace(0,np.pi,nbins))
        plt.title("theta distribution")
        plt.xlabel("theta value")
        plt.ylabel("distribution")
        plt.savefig("theta_distribution.png")
        plt.close()
    #check if the probability computed for each bin is within the tolerance
    success = True
    tolerance = max(vals_theta)*tolerance
    for k in range(0,len(vals_theta)):
        prob_hist = vals_theta[k]/(sum(vals_theta))
        prob_test = -(np.cos(bins_theta[k+1])-np.cos(bins_theta[k]))/2
        if abs(prob_test-prob_hist)>tolerance:
            success = False
            break

    #re-seed the random number generator
    rd.seed()

    return success

def test_rand_velocity(sigma, num_sample=10000, nbins=20, tolerance=1e-3, seed="Dimitris", plot=False, save=True):
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
        - save: saves the plot
    Returns: True if the test is succesful, False otherwise
    """
    rd.seed(seed)
    velocity_array = np.zeros(num_sample)
    for k in range(0,len(velocity_array)):
        velocity_array[k] = kicks.rand_velocity(sigma)
    
    #do a histogram
    vals_velocity, bins_velocity = np.histogram(velocity_array, bins=np.linspace(0,3*sigma,nbins))
    if plot:    
        plt.hist(velocity_array, bins=np.linspace(0,3*sigma,nbins))
        plt.title("velocity distribution")
        plt.xlabel("velocity value")
        plt.ylabel("distribution")
        plt.show()
        plt.close()
    if save:
        plt.hist(velocity_array, bins=np.linspace(0,3*sigma,nbins))
        plt.title("velocity distribution")
        plt.xlabel("velocity value")
        plt.ylabel("distribution")
        plt.savefig("velocity_distribution.png")
        plt.close()
    #check if the probability computed from each bin is within the tolerance
    success = True
    tolerance = max(vals_velocity)*tolerance
    for k in range(0,len(vals_velocity)):
        prob_hist = vals_velocity[k]/(sum(vals_velocity))
        prob_test = maxwell.cdf(bins_velocity[k+1],0,sigma) - maxwell.cdf(bins_velocity[k+1],0,sigma)
        if abs(prob_test-prob_hist)>tolerance:
            success = False
            break
    
    #re-seed the random number generator
    rd.seed()
    
    return success

def testing_circular_function_momentum(Ai, M1, M2, Mns, test_sigma, num_sample=1000, seed = "Lela", tolerance=1e-3):
    """Test that the post_explosion_params_circular function produces
    a correct momentum against a calculated momentum
    Arguments:
        - Ai: the initial maximum separation of the two masses
            pre-explosion
        - M1: solar mass of the first mass pre-explosion
        - M2: solar mass of the second mass
        - Mns: solar mass of the first mass post-explosion
        - test_sigma: a sample sigma for the rand_velocity function
        - num_sample: number of points sampled
        - seed: the seed used for the random number generator
        - tolerance: tolerance for the test
    Returns: True or False as to whether the test was successful
    """  
    rd.seed(seed)

    
    for i in range(num_sample):
        #establishing random parameters
        Vk = kicks.rand_velocity(test_sigma)*1e5
        theta = kicks.rand_theta()
        phi = kicks.rand_phi()
        
        #getting values from the post_explosion_params_circular function
        separation, e, angle, boolean = kicks.post_explosion_params_circular(Ai, M1, M2, Mns,theta,phi,Vk)
        
        #calculating the momentum using the results from the function
        Momentum_function = Mns*Msun*M2*Msun*np.sqrt(cgrav*separation*Rsun*(1-e**2)/((Mns+M2)*Msun))
        
        #Calculating the momentum without using the results of the function
        
        #establishing angular velocity and separations from the center of
            #mass in the center of mass frame
        omega = np.sqrt(cgrav*Msun*(M1+M2)/(Ai*Rsun)**3) #rad/second
        R2 = Ai*Rsun*Mns/(Mns+M2) #cm
        R1 = Ai*Rsun - R2 #cm
        
        #velocities of the masses before the kick 
        V1_initial = M2/(M1+M2)*omega*Ai*Rsun #cm/second
        V2 = M1/(M1+M2)*omega*Ai*Rsun #cm/second,-y direction
        
        #velocities after the kick, V2 is unaffected by the kick
        V1x_final = Vk*1e5*np.sin(phi)*np.cos(theta)
        V1y_final = V1_initial + Vk*1e5*np.cos(theta)
        V1z_final = Vk*1e5*np.sin(phi)*np.sin(theta)
        
        #calculating the components of the angular momentum using the cross product
        Momentum_1y = -R1*V1z_final*Mns*Msun
        Momentum_1z = R1*V1y_final*Mns*Msun
        Momentum_2 = R2*V2*M2*Msun #z direction
        
        #absolute value of the components of the angular momentum
        Momentum_calculated = np.sqrt(Momentum_1y**2 + Momentum_1z**2 + Momentum_2**2)
        
        
        #checking that the two momentums are relatively equal
        if abs(Momentum_function - Momentum_calculated)/Momentum_function>tolerance:
            return False
            
    rd.seed()
        
    return True
        






        
        