#!/usr/bin/env python2
import random as rd
import numpy as np
import matplotlib.pyplot as plt
#from binary_tools.utils.constants import *
from scipy.stats import maxwell
from scipy.integrate import quad

__author__ = "Kaliroe Pappas"
__credits__ = ["Kaliroe Pappas"]
__license__ = "GPL"
__version__ = "3.0"
__maintainer__ = "Kaliroe Pappas"
__email__ = "kmpappa2@illinois.edu"

"""Variety of functions that compute probabilities
of binary systems initial setup
"""

def default_function(M,alpha):
    return M**(-alpha)

def default_function_log(M,alpha):
    return 10**(M*(1-alpha))

def prob_mass(M1,M2,Mlow,Mup,alpha=2.3,f_M=default_function):
    """computes the probability for a star
    to be inbetween M1 and M2 solar masses
    Arguments:
        - M1: the lower mass in solar masses
        - M2: the upper mass in solar masses
        - Mlow: the lowest absolute mass in solar masses
        - Mup: the highest absolute mass in solar masses
        - alpha: the exponential value in the function
        - f_M: the function used to calculate the probability
    Returns the probability
    """
    alpha = alpha
    FM = lambda Mtemp: f_M(Mtemp,alpha)
    normalize = quad(FM,Mlow,Mup)[0]
    probability = quad(FM,M1,M2)[0]/normalize
    return probability

def prob_mass_log(log_M1,log_M2,Mlow,Mup,alpha=2.3,f_M=default_function_log):
    """computes the probability for a star to be inbetween 
    M1 and M2 solar masses in a log base 10 scale
    Arguments:
        - log_M1: the lower mass in log base 10 solar masses
        - log_M2: the upper mass in log base 10 solar masses
        - Mlow: the lowest absolute mass in solar masses
        - Mup: the highest absolute mass in solar masses
        - alpha: the exponential value in the function
        - f_M: the function used to calculate the probability
    Returns the probability
    """
    alpha = alpha
    FM = lambda Mtemp: f_M(Mtemp,alpha)
    normalize = quad(FM,Mlow,Mup)[0]
    probability = quad(FM,log_M1,log_M2)[0]/normalize
    return probability


def prob_mass_ratio(R1,R2,Rlow,Rup,gamma=0,f_M=default_function):
    """computes the probability for a the ratio
    of the two masses with the assumption that the 
    larger mass is in the denominator making the
    ratio < 1
    Arguments:
        - R1: the lower ratio
        - R2: the upper ratio
        - Rlow: the lowest absolute ratio
        - Rup: the highest absolute ratio
        - gamma: the exponential value in the function
        - f_M: the function used to calculate the probability
    Returns the probability
    """
    alpha = gamma
    FM = lambda Mtemp: f_M(Mtemp,alpha)
    normalize = quad(FM,Rlow,Rup)[0]
    probability = quad(FM,R1,R2)[0]/normalize
    return probability


def prob_period_log(log_P1,log_P2,Plow,Pup,beta=0,f_M=default_function_log):
    """computes the probability for a binary system
    the have a period between log_P1 and log_P2
    Arguments:
        - log_P1: the lower period in log base 10
        - log_P2: the upper period in log base 10
        - Plow: the lowest absolute period
        - Pup: the highest absolute period
        - beta: the exponential value in the function
        - f_M: the function used to calculate the probability
    Returns the probability
    """
    alpha = beta
    FM = lambda Mtemp: f_M(Mtemp,alpha)
    normalize = quad(FM,Plow,Pup)[0]
    probability = quad(FM,log_P1,log_P2)[0]/normalize
    return probability


"""A series of random generating functions"""


def rand_mass(Mlow,Mup,alpha=2.3,f_M=default_function):
    """Randomly samples the possible range of
    masses using a Montecarlo sampling method 
    and drawing a random number from an array of 
    possibilities using the distribution:
        f(x) = M**(-alpha)
    Arguments:
        - Mlow: the lowest absolute mass in solar masses
        - Mup: the highest absolute mass in solar masses
        - alpha: the exponential value in the function
        - f_M: the function used to calculate the probability
    Returns a random mass
    """
    alpha = alpha
    
    FM = lambda Mtemp: f_M(Mtemp,alpha)
    normalize = quad(FM,Mlow,Mup)[0]
    
    max = FM(Mlow)/normalize
    
    while True:
        randx = rd.uniform(Mlow,Mup)
        randy = rd.uniform(0, max)
        
        if (f_M(randx,alpha)/normalize) >randy:
            return randx

def rand_mass_log(Mlow,Mup,alpha=2.3,f_M=default_function_log):
    """Randomly samples the possible range of
    masses using a Montecarlo sampling method 
    and drawing a random number from an array of 
    possibilities using the distribution:
        f(x) = 10**(M*(1-beta))
    Arguments:
        - Mlow: the lowest absolute mass in solar masses
        - Mup: the highest absolute mass in solar masses
        - alpha: the exponential value in the function
        - f_M: the function used to calculate the probability
    Returns a random mass
    """
    alpha = alpha
    
    FM = lambda Mtemp: f_M(Mtemp,alpha)
    normalize = quad(FM,Mlow,Mup)[0]
    
    max = FM(Mlow)/normalize
    
    while True:
        randx = rd.uniform(Mlow,Mup)
        randy = rd.uniform(0, max)
        
        if (f_M(randx,alpha)/normalize) >randy:
            return 10**randx


def rand_ratio(Rlow,Rup,gamma=0,f_M=default_function):
    """Randomly samples the possible range of
    ratios using a Montecarlo sampling method 
    and drawing a random number from an array of 
    possibilities using the distribution:
        f(x) = M**(-gamma)
    Arguments:
        - Rlow: the lowest absolute ratio
        - Rup: the highest absolute ratio
        - gammma: the exponential value in the function
        - f_M: the function used to calculate the probability
    Returns a random ratio
    """
    alpha = gamma
    
    FM = lambda Mtemp: f_M(Mtemp,alpha)
    normalize = quad(FM,Rlow,Rup)[0]
    
    max = FM(Rlow)/normalize
    
    while True:
        randx = rd.uniform(Rlow,Rup)
        randy = rd.uniform(0, max)
        
        if (f_M(randx,alpha)/normalize) >randy:
            return randx

def rand_period_log(Plow,Pup,beta=0,f_M=default_function_log):
    """Randomly samples the possible range of
    periods using a Montecarlo sampling method 
    and drawing a random number from an array of 
    possibilities using the distribution:
        f(x) = 10**(M*(1-beta))
    Arguments:
        - Plow: the lowest absolute period
        - Pup: the highest absolute period
        - beta: the exponential value in the function
        - f_M: the function used to calculate the probability
    Returns a random period in days
    """
    alpha = beta
    
    FM = lambda Mtemp: f_M(Mtemp,alpha)
    normalize = quad(FM,Plow,Pup)[0]
    
    if (1-alpha)>0:
        M=Pup
    else:
        M=Plow
    
    max = FM(M)/normalize
    
    while True:
        randx = rd.uniform(Plow,Pup)
        randy = rd.uniform(0, max)
        
        if (f_M(randx,alpha)/normalize) >randy:
            return 10**randx



"""Alternate rand functions"""

def rand_mass_cdf(Mlow,Mup,alpha=2.3,f_M=default_function):
    """Randomly samples the possible range of
    masses using the inverse cdf method and 
    the distribution:
        f(x) = M**(-alpha) 
    Arguments:
        - Mlow: the lowest absolute mass in solar masses
        - Mup: the highest absolute mass in solar masses
        - alpha: the exponential value in the function
        - f_M: the function used to calculate the probability
    Returns a random mass
    """
    alpha = alpha
    
    FM = lambda Mtemp: f_M(Mtemp,alpha)
    normalize = quad(FM,Mlow,Mup)[0]
    
    rand_mass = ((1-alpha)*normalize*rd.random()+Mlow**(1-alpha))**(1/(1-alpha))
    return rand_mass

def rand_mass_log_cdf(Mlow,Mup,alpha=2.3,f_M=default_function_log):
    """Randomly samples the possible range of
    masses using the inverse cdf method and 
    the distribution:
        f(x) = 10**(M*(1-beta))
    Arguments:
        - Mlow: the lowest absolute mass in solar masses
        - Mup: the highest absolute mass in solar masses
        - alpha: the exponential value in the function
        - f_M: the function used to calculate the probability
    Returns a random mass
    """
    alpha = alpha
    
    FM = lambda Mtemp: f_M(Mtemp,alpha)
    normalize = quad(FM,Mlow,Mup)[0]
    
    rand_mass_log = np.log10(normalize*(1-alpha)*np.log(10)*rd.random()+10**((1-alpha)*Mlow))/(1-alpha)
    rand_mass_log = 10**(rand_mass_log)
    return rand_mass_log

def rand_ratio_cdf(Rlow,Rup,gamma=0,f_M=default_function):
    """Randomly samples the possible range of
    ratios using the inverse cdf method and 
    the distribution:
        f(x) = M**(-gamma)
    Arguments:
        - Rlow: the lowest absolute ratio
        - Rup: the highest absolute ratio
        - gammma: the exponential value in the function
        - f_M: the function used to calculate the probability
    Returns a random ratio
    """
    alpha = gamma
    
    FM = lambda Mtemp: f_M(Mtemp,alpha)
    normalize = quad(FM,Rlow,Rup)[0]
    
    rand_ratio = ((1-alpha)*normalize*rd.random()+Rlow**(1-alpha))**(1/(1-alpha))
    return rand_ratio
        
def rand_period_log_cdf(Plow,Pup,beta=0,f_M=default_function_log):
    """Randomly samples the possible range of
    periods using the inverse cdf method and 
    the distribution:
        f(x) = 10**(M*(1-beta))
    Arguments:
        - Plow: the lowest absolute period
        - Pup: the highest absolute period
        - beta: the exponential value in the function
        - f_M: the function used to calculate the probability
    Returns a random period
    """
    alpha = beta
    
    FM = lambda Mtemp: f_M(Mtemp,alpha)
    normalize = quad(FM,Plow,Pup)[0]
    
    rand_period_log = np.log10(normalize*(1-alpha)*np.log(10)*rd.random()+10**((1-alpha)*Plow))/(1-alpha)
    rand_period_log = 10**(rand_period_log)
    return rand_period_log
        
        
"""Tests for both the probability functions and the random functions"""


def test_mass_functions(num_sample=10000, tolerance = 1e-3, seed="Kitt"):
    """Test that prob_mass and prob_mass_log return
    equal values for scaled values
    Arguments:
        - num_sample: number of random phi generated
        - tolerance: tolerance for the test
        - seed: the seed used for the random number generator
    Returns: True if the test is succesful, False otherwise
    """
    rd.seed(seed)
    success = True
    
    for mass in range(1,100):
        test_mass = prob_mass(1,mass,1,100)
        test_mass_log = prob_mass_log(0,np.log10(mass),0,2)
        if abs(test_mass_log - test_mass)>tolerance:
            print(mass)
            success = False
    
    return success

def test_rand_mass(Mlow,Mup,num_sample=10000, nbins=20, tolerance = 1e-3, seed="Feyra", cdf=True, plot=False, save=True):
    """Test that rand_mass is sampled as the expected
    from prob_mass
    Arguments:
        - M1: the lower mass in solar masses
        - M2: the upper mass in solar masses
        - num_sample: number of random phi generated
        - nbins: random sampled numbers will be binned to compute
            probabilities and compare to expected values. This variable
            specifies the number of bins used.
        - tolerance: tolerance for the test
        - seed: the seed used for the random number generator
        - cdf: uses the alterate random function that uses the
            cdf Montecarlo method
        - plot: if true, plot results
        - save: if true, saves the plot
    Returns: True if the test is succesful, False otherwise
    """
    rd.seed(seed)
    mass_array = np.zeros(num_sample)
    if cdf:
        for i in range(0,len(mass_array)):
            mass_array[i] = rand_mass_cdf(Mlow,Mup)
    else:
        for i in range(0,len(mass_array)):
            mass_array[i] = rand_mass(Mlow,Mup)
    
    #do a histogram
    #TODO: use numpy histogram, to avoid plotting if its not neccesary
    vals_mass, bins_mass = np.histogram(mass_array, bins=np.linspace(Mlow,Mup,nbins))
    
    prob_test = np.zeros(len(vals_mass))
    
    for k in range(0,len(vals_mass)):
        prob_test[k] = prob_mass(bins_mass[k],bins_mass[k+1],Mlow,Mup)
    
    test_array = []
    
    for j in range(0,len(prob_test)):
        test_array = np.append(test_array, np.ones(int(round(prob_test[j]*num_sample)))*bins_mass[j])
        
    if plot:
        plt.hist(mass_array, bins=np.linspace(Mlow,Mup,nbins), alpha = 0.5, label = "function output")
        plt.hist(test_array, bins=np.linspace(Mlow,Mup,nbins), alpha = 0.5, label = "expected value")
        plt.title("mass distribution")
        plt.xlabel("mass value")
        plt.ylabel("distribution")
        plt.legend(loc='upper right')
        plt.show()
        plt.close()
    if save:
        plt.hist(mass_array, bins=np.linspace(Mlow,Mup,nbins), alpha = 0.5, label = "function output")
        plt.hist(test_array, bins=np.linspace(Mlow,Mup,nbins), alpha = 0.5, label = "expected value")
        plt.title("mass distribution")
        plt.xlabel("mass value")
        plt.ylabel("distribution")
        plt.legend(loc='upper right')
        plt.savefig("mass_distribution.png")
        plt.close()

    #check if the probability computed for each bin is within the tolerance
    success = True
    tolerance = max(vals_mass)*tolerance
    for k in range(0,len(vals_mass)):
        prob_hist = vals_mass[k]/(sum(vals_mass))
        if abs(prob_test[k]-prob_hist)>tolerance:
            success = False
            break

    #re-seed the random number generator
    rd.seed()

    return success



def test_rand_mass_log(Mlow,Mup,num_sample=10000, nbins=20, tolerance = 1e-3, seed="Stewart", cdf=True, plot=False, save=True):
    """Test that rand_mass is sampled as the expected
    from prob_mass
    Arguments:
        - Mlow: the lower mass in solar masses
        - Mup: the upper mass in solar masses
        - num_sample: number of random phi generated
        - nbins: random sampled numbers will be binned to compute
            probabilities and compare to expected values. This variable
            specifies the number of bins used.
        - tolerance: tolerance for the test
        - seed: the seed used for the random number generator
        - cdf: uses the alterate random function that uses the
            cdf Montecarlo method
        - plot: if true, plot results
        - save: if true, saves the plot
    Returns: True if the test is succesful, False otherwise
    """
    rd.seed(seed)
    mass_log_array = np.zeros(num_sample)
    if cdf:
        for i in range(0,len(mass_log_array)):
            mass_log_array[i] = np.log10(rand_mass_log_cdf(Mlow,Mup))
    else:
        for i in range(0,len(mass_log_array)):
            mass_log_array[i] = np.log10(rand_mass_log(Mlow,Mup))
    
    #do a histogram
    #TODO: use numpy histogram, to avoid plotting if its not neccesary
    vals_mass_log, bins_mass_log = np.histogram(mass_log_array, bins=np.linspace(Mlow,Mup,nbins))
    
    prob_test = np.zeros(len(vals_mass_log))
    
    for k in range(0,len(vals_mass_log)):
        prob_test[k] = prob_mass_log(bins_mass_log[k],bins_mass_log[k+1],Mlow,Mup)
    
    test_array = []
    
    for j in range(0,len(prob_test)):
        test_array = np.append(test_array, np.ones(int(round(prob_test[j]*num_sample)))*bins_mass_log[j])
        
    if plot:
        plt.hist(mass_log_array, bins=np.linspace(Mlow,Mup,nbins), alpha = 0.5, label = "function output")
        plt.hist(test_array, bins=np.linspace(Mlow,Mup,nbins), alpha = 0.5, label = "expected value")
        plt.gca().set_xscale("log")
        plt.title("mass log distribution")
        plt.xlabel("mass log value")
        plt.ylabel("distribution")
        plt.legend(loc='upper right')
        plt.show()
        plt.close()
    if save:
        plt.hist(mass_log_array, bins=np.linspace(Mlow,Mup,nbins), alpha = 0.5, label = "function output")
        plt.hist(test_array, bins=np.linspace(Mlow,Mup,nbins), alpha = 0.5, label = "expected value")
        plt.gca().set_xscale("log")
        plt.title("mass log distribution")
        plt.xlabel("mass log value")
        plt.ylabel("distribution")
        plt.legend(loc='upper right')
        plt.savefig("mass_log_distribution.png")
        plt.close()

    #check if the probability computed for each bin is within the tolerance
    success = True
    tolerance = max(vals_mass_log)*tolerance
    for k in range(0,len(vals_mass_log)):
        prob_hist = vals_mass_log[k]/(sum(vals_mass_log))
        if abs(prob_test[k]-prob_hist)>tolerance:
            success = False
            break

    #re-seed the random number generator
    rd.seed()

    return success

def test_rand_ratio(Rlow, Rup, test_alpha=2, num_sample=10000, nbins=20, tolerance = 1e-3, seed="Lee", cdf=True, plot=False, save=True):
    """Test that rand_ratio is sampled as the expected
    from prob_ratio
    Arguments:
        - Rlow: the lowest absolute ratio
        - Rup: the highest absolute ratio
        - test_alpha: a test alpha for rand_ratio and prob_ratio
        - num_sample: number of random phi generated
        - nbins: random sampled numbers will be binned to compute
            probabilities and compare to expected values. This variable
            specifies the number of bins used.
        - tolerance: tolerance for the test
        - seed: the seed used for the random number generator
        - cdf: uses the alterate random function that uses the
            cdf Montecarlo method
        - plot: if true, plot results
        - save: if true, saves the plot
    Returns: True if the test is succesful, False otherwise
    """
    rd.seed(seed)
    ratio_array = np.zeros(num_sample)
    if cdf:
        for i in range(0,len(ratio_array)):
            ratio_array[i] = rand_ratio_cdf(Rlow,Rup,gamma=test_alpha)
    else:
        for i in range(0,len(ratio_array)):
            ratio_array[i] = rand_ratio(Rlow,Rup,gamma=test_alpha)
    #do a histogram
    #TODO: use numpy histogram, to avoid plotting if its not neccesary
    vals_ratio, bins_ratio = np.histogram(ratio_array, bins=np.linspace(Rlow,Rup,nbins))
    
    prob_test = np.zeros(len(vals_ratio))
    
    for k in range(0,len(vals_ratio)):
        prob_test[k] = prob_mass_ratio(bins_ratio[k],bins_ratio[k+1],Rlow,Rup,gamma=test_alpha)
    
    test_array = []
    
    for j in range(0,len(prob_test)):
        test_array = np.append(test_array, np.ones(int(round(prob_test[j]*num_sample)))*bins_ratio[j])
        
    if plot:
        plt.hist(ratio_array, bins=np.linspace(Rlow,Rup,nbins), alpha = 0.5, label = "function output")
        plt.hist(test_array, bins=np.linspace(Rlow,Rup,nbins), alpha = 0.5, label = "expected value")
        plt.title("ratio distribution")
        plt.xlabel("ratio value")
        plt.ylabel("distribution")
        plt.legend(loc='upper right')
        plt.show()
        plt.close()
    if save:
        plt.hist(ratio_array, bins=np.linspace(0,1,nbins), alpha = 0.5, label = "function output")
        plt.hist(test_array, bins=np.linspace(0,1,nbins), alpha = 0.5, label = "expected value")
        plt.title("ratio distribution")
        plt.xlabel("ratio value")
        plt.ylabel("distribution")
        plt.legend(loc='upper right')
        plt.savefig("ratio_distribution.png")
        plt.close()

    #check if the probability computed for each bin is within the tolerance
    success = True
    tolerance = max(vals_ratio)*tolerance
    for k in range(0,len(vals_ratio)):
        prob_hist = vals_ratio[k]/(sum(vals_ratio))
        if abs(prob_test[k]-prob_hist)>tolerance:
            success = False
            break

    #re-seed the random number generator
    rd.seed()

    return success

def test_rand_period_log(Plow,Pup,test_alpha=0,num_sample=10000,nbins=20,tolerance = 1e-3,seed="Hardy",cdf=True,plot=False,save=True):
    """Test that rand_period is sampled as the expected
    from prob_period_log
    Arguments:
        - Plow: the lowest absolute ratio
        - Pup: the highest absolute ratio
        - num_sample: number of random phi generated
        - nbins: random sampled numbers will be binned to compute
            probabilities and compare to expected values. This variable
            specifies the number of bins used.
        - tolerance: tolerance for the test
        - seed: the seed used for the random number generator
        - cdf: uses the alterate random function that uses the
            cdf Montecarlo method
        - plot: if true, plot results
        - save: if true, saves the plot
    Returns: True if the test is succesful, False otherwise
    """
    rd.seed(seed)
    period_array = np.zeros(num_sample)
    if cdf:
        for i in range(0,len(period_array)):
            period_array[i] = np.log10(rand_period_log_cdf(Plow,Pup,beta=test_alpha))
    else:
        for i in range(0,len(period_array)):
            period_array[i] = np.log10(rand_period_log(Plow,Pup,beta=test_alpha))
    
    #do a histogram
    #TODO: use numpy histogram, to avoid plotting if its not neccesary
    vals_period, bins_period = np.histogram(period_array, bins=np.linspace(Plow,Pup,nbins))
    
    prob_test = np.zeros(len(vals_period))
    
    for k in range(0,len(vals_period)):
        prob_test[k] = prob_period_log(bins_period[k],bins_period[k+1],Plow,Pup,beta=test_alpha)
    
    test_array = []
    
    for j in range(0,len(prob_test)):
        test_array = np.append(test_array, np.ones(int(round(prob_test[j]*num_sample)))*bins_period[j])
        
    if plot:
        plt.hist(period_array, bins=np.linspace(Plow,Pup,nbins), alpha = 0.5, label = "function output")
        plt.hist(test_array, bins=np.linspace(Plow,Pup,nbins), alpha = 0.5, label = "expected value")
        plt.gca().set_xscale("log")
        plt.title("period distribution")
        plt.xlabel("period value")
        plt.ylabel("distribution")
        plt.legend(loc='upper right')
        plt.show()
        plt.close()
    if save:
        plt.hist(period_array, bins=np.linspace(Plow,Pup,nbins), alpha = 0.5, label = "function output")
        plt.hist(test_array, bins=np.linspace(Plow,Pup,nbins), alpha = 0.5, label = "expected value")
        plt.gca().set_xscale("log")
        plt.title("period distribution")
        plt.xlabel("period value")
        plt.ylabel("distribution")
        plt.legend(loc='upper right')
        plt.savefig("period_distribution.png")
        plt.close()

    #check if the probability computed for each bin is within the tolerance
    success = True
    tolerance = max(vals_period)*tolerance
    for k in range(0,len(vals_period)):
        prob_hist = vals_period[k]/(sum(vals_period))
        if abs(prob_test[k]-prob_hist)>tolerance:
            success = False
            break

    #re-seed the random number generator
    rd.seed()

    return success

