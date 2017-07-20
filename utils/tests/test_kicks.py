#!/usr/bin/env python
from binary_tools.utils.constants import *
from binary_tools.utils import kicks
from binary_tools.utils.Keplers_laws import *
import matplotlib.pyplot as plt
from scipy.stats import maxwell
from scipy.integrate import quad
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
        - save: if true, saves the plot
    Returns: True if the test is succesful, False otherwise
    """
    rd.seed(seed)
    phi_array = np.zeros(num_sample)
    for i in range(0,len(phi_array)):
        phi_array[i] = kicks.rand_phi()
    
    #do a histogram
    #TODO: use numpy histogram, to avoid plotting if its not neccesary
    vals_phi, bins_phi = np.histogram(phi_array, bins=np.linspace(0,2*np.pi,nbins))
    
    prob_test = np.zeros(len(vals_phi))
    
    for k in range(0,len(vals_phi)):
        prob_test[k] = (1/(2*np.pi))*(bins_phi[k+1]-bins_phi[k])
    
    test_array = []
    
    for j in range(0,len(prob_test)):
        test_array = np.append(test_array, np.ones(int(round(prob_test[j]*num_sample)))*bins_phi[j])
        
    if plot:
        plt.hist(phi_array, bins=np.linspace(0,2*np.pi,nbins), alpha = 0.5, label = "function")
        plt.hist(test_array, bins=np.linspace(0,2*np.pi,nbins), alpha = 0.5, label = "comparison")
        plt.title("phi distribution")
        plt.xlabel("phi value")
        plt.ylabel("distribution")
        plt.legend(loc='upper right')
        plt.show()
        plt.close()
    if save:
        plt.hist(phi_array, bins=np.linspace(0,2*np.pi,nbins), alpha = 0.5, label = "function")
        plt.hist(test_array, bins=np.linspace(0,2*np.pi,nbins), alpha = 0.5, label = "comparison")
        plt.title("phi distribution")
        plt.xlabel("phi value")
        plt.ylabel("distribution")
        plt.legend(loc='upper right')
        plt.savefig("phi_distribution.png")
        plt.close()

    #check if the probability computed for each bin is within the tolerance
    success = True
    tolerance = max(vals_phi)*tolerance
    for k in range(0,len(vals_phi)):
        prob_hist = vals_phi[k]/(sum(vals_phi))
        if abs(prob_test[k]-prob_hist)>tolerance:
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
        - save: if true, saves the plot
    Returns: True if the test is succesful, False otherwise
    """
    rd.seed(seed)
    theta_array = np.zeros(num_sample)
    for i in range(0,len(theta_array)):
        theta_array[i] = kicks.rand_theta()

    #do a histogram
    #TODO: use numpy histogram, to avoid plotting if its not neccesary
    vals_theta, bins_theta = np.histogram(theta_array, bins=np.linspace(0,np.pi,nbins))
    
    prob_test = np.zeros(len(vals_theta))
    
    for k in range(0,len(vals_theta)):
        prob_test[k] = -(np.cos(bins_theta[k+1])-np.cos(bins_theta[k]))/2
    
    test_array = []
    
    for j in range(0,len(prob_test)):
        test_array = np.append(test_array, np.ones(int(round(prob_test[j]*num_sample)))*bins_theta[j])
    
    if plot:
        plt.hist(theta_array, bins=np.linspace(0,np.pi,nbins), alpha = 0.5, label = "function")
        plt.hist(test_array, bins=np.linspace(0,np.pi,nbins), alpha = 0.5, label = "comparison")
        plt.title("theta distribution")
        plt.xlabel("theta value")
        plt.ylabel("distribution")
        plt.legend(loc='upper right')
        plt.show()
        plt.close()
    if save:
        plt.hist(theta_array, bins=np.linspace(0,np.pi,nbins), alpha = 0.5, label = "function")
        plt.hist(test_array, bins=np.linspace(0,np.pi,nbins), alpha = 0.5, label = "comparison")
        plt.title("theta distribution")
        plt.xlabel("theta value")
        plt.ylabel("distribution")
        plt.legend(loc='upper right')
        plt.savefig("theta_distribution.png")
        plt.close()
    #check if the probability computed for each bin is within the tolerance
    success = True
    tolerance = max(vals_theta)*tolerance
    for k in range(0,len(vals_theta)):
        prob_hist = vals_theta[k]/(sum(vals_theta))
        if abs(prob_test[k]-prob_hist)>tolerance:
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
        - save: if true, saves the plot
    Returns: True if the test is succesful, False otherwise
    """
    rd.seed(seed)
    velocity_array = np.zeros(num_sample)
    for k in range(0,len(velocity_array)):
        velocity_array[k] = kicks.rand_velocity(sigma)
    
    #do a histogram
    vals_velocity, bins_velocity = np.histogram(velocity_array, bins=np.linspace(0,3*sigma,nbins))
    
    prob_test = np.zeros(len(vals_velocity))
    
    for k in range(0,len(vals_velocity)):
        prob_test[k] = maxwell.cdf(bins_velocity[k+1],0,sigma) - maxwell.cdf(bins_velocity[k],0,sigma)
    
    test_array = []
    
    for j in range(0,len(prob_test)):
        test_array = np.append(test_array, np.ones(int(round(prob_test[j]*num_sample)))*bins_velocity[j])
    
    if plot:    
        plt.hist(velocity_array, bins=np.linspace(0,3*sigma,nbins), alpha = 0.5, label = "function")
        plt.hist(test_array, bins=np.linspace(0,3*sigma,nbins), alpha = 0.5, label = "comparison")
        plt.title("velocity distribution")
        plt.xlabel("velocity value")
        plt.ylabel("distribution")
        plt.legend(loc='upper right')
        plt.show()
        plt.close()
    if save:
        plt.hist(velocity_array, bins=np.linspace(0,3*sigma,nbins), alpha = 0.5, label = "function")
        plt.hist(test_array, bins=np.linspace(0,3*sigma,nbins), alpha = 0.5, label = "comparison")
        plt.title("velocity distribution")
        plt.xlabel("velocity value")
        plt.ylabel("distribution")
        plt.legend(loc='upper right')
        plt.savefig("velocity_distribution.png")
        plt.close()
    #check if the probability computed from each bin is within the tolerance
    success = True
    tolerance = max(vals_velocity)*tolerance
    for k in range(0,len(vals_velocity)):
        prob_hist = vals_velocity[k]/(sum(vals_velocity))
        if abs(prob_test[k]-prob_hist)>tolerance:
            success = False
            break
    
    #re-seed the random number generator
    rd.seed()
    
    return success



def test_rand_true_anomaly(e,num_sample=10000, nbins=20, tolerance = 1e-3, seed="Rhysand", plot=False, save=True):
    """Test that phi is sampled as flat from zero to pi
    Arguments:
        - e: eccentricity of the orbit
        - num_sample: number of random phi generated
        - nbins: random sampled numbers will be binned to compute
            probabilities and compare to expected values. This variable
            specifies the number of bins used.
        - tolerance: tolerance for the test
        - seed: the seed used for the random number generator
        - plot: if true, plot results
        - save: if true, saves the plot
    Returns: True if the test is succesful, False otherwise
    """
    rd.seed(seed)
    true_anomaly_array = np.zeros(num_sample)
    for i in range(0,len(true_anomaly_array)):
        true_anomaly_array[i] = kicks.rand_true_anomaly(e)

    #do a histogram
    #TODO: use numpy histogram, to avoid plotting if its not neccesary
    vals_true_anomaly, bins_true_anomaly = np.histogram(true_anomaly_array, bins=np.linspace(0,2*np.pi,nbins))

    prob_test = np.zeros(len(vals_true_anomaly))

    def func(x):
        return np.sqrt(1/(4*np.pi))*(1-e**2)**(1.5)/((1+e*np.cos(x))**2)
    
    normalize = quad(func,bins_true_anomaly[0],bins_true_anomaly[nbins-1])[0]
    
    for k in range(0,len(vals_true_anomaly)):
        prob_test[k] = quad(func,bins_true_anomaly[k],bins_true_anomaly[k+1])[0]
    
    prob_test = prob_test/normalize
    
    test_array = []
    
    for j in range(0,len(prob_test)):
        test_array = np.append(test_array, np.ones(int(round(prob_test[j]*num_sample)))*bins_true_anomaly[j])
    

    if plot:
        plt.hist(true_anomaly_array, bins=np.linspace(0,2*np.pi,nbins), alpha = 0.5, label = "function")
        plt.hist(test_array, bins=np.linspace(0,2*np.pi,nbins), alpha = 0.5, label = "comparison")
        plt.title("true anomaly distribution")
        plt.xlabel("true anomaly value")
        plt.ylabel("distribution")
        plt.legend(loc='upper right')
        plt.show()
        plt.close()
    if save:
        plt.hist(true_anomaly_array, bins=np.linspace(0,2*np.pi,nbins), alpha = 0.5, label = "function")
        plt.hist(test_array, bins=np.linspace(0,2*np.pi,nbins), alpha = 0.5, label = "comparison")
        plt.title("true anomaly distribution")
        plt.xlabel("true anomaly value")
        plt.ylabel("distribution")
        plt.legend(loc='upper right')
        plt.savefig("true_anomaly_distribution.png")
        plt.close()

    #check if the probability computed for each bin is within the tolerance
    success = True
    tolerance = max(vals_true_anomaly)*tolerance
    
    for k in range(0,len(vals_true_anomaly)):
        prob_hist = vals_true_anomaly[k]/(sum(vals_true_anomaly))
        if abs(prob_test[k]-prob_hist)>tolerance:
            success = False
            break

    #re-seed the random number generator
    rd.seed()

    return success






def testing_circular_function_momentum(Ai=133, M1=5.5, M2=55, Mns=1.4, test_sigma=100, num_sample=1000, seed = "Lela", tolerance=1e-3):
    """Test that the post_explosion_params_circular function produces
    a correct angular momentum against a calculated angular momentum. This 
    angular momentum is calculated by first finding the anglar velocity, and 
    the individual relative velocities in the center of mass frame pre-
    explosion. The components of the kick velocity were then added to the 
    velocity of mass 1, which is the super nova mass. The separation in the 
    final center of mass frame accounting for the loss of mass of mass 1 is 
    crossed with the velocities to get the components of angular momentum. The 
    total angular momentum is calculated by taking the absolute value of these 
    components. 
    Arguments:
        - Ai: the initial semi-major axis of the system
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
        
        #establishing angular velocity 
        omega = np.sqrt(cgrav*Msun*(M1+M2)/(Ai*Rsun)**3) #rad/second
        
        #velocities of the masses before the kick 
        V1_initial = M2/(M1+M2)*omega*Ai*Rsun #cm/second
        V2 = M1/(M1+M2)*omega*Ai*Rsun #cm/second,-y direction
        
        #velocities after the kick, V2 is unaffected by the kick
        V1x_final = Vk*1e5*np.sin(phi)*np.cos(theta)
        V1y_final = V1_initial + Vk*1e5*np.cos(theta)
        V1z_final = Vk*1e5*np.sin(phi)*np.sin(theta)
        
        #separations from the center of mass in the center of mass frame post-explosion
        R2 = Ai*Rsun*Mns/(Mns+M2) #cm
        R1 = Ai*Rsun - R2 #cm
        
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
        



def testing_circular_function_graph(test_sigma = 100, test_M1 = 5.5, test_M2 = 55, test_Ai = 133, test_Mns = 1.4, seed="Flay",sample_velocity = 100, npoints =10000, plot=False, save =True):
    """Test that the graph of the eccentricity vs the period looks correct
    Arguments:
        - test_sigma: a sample sigma for the rand_velocity function
        - test_M1: solar mass of the first mass pre-explosion
        - test_M2: solar mass of the second mass
        - test_Ai: the initial semi-major axis of the system
            pre-explosion
        - test_Mns: solar mass of the first mass post-explosion
        - seed: the seed used for the random number generator
        - sample_velocity: a constant velocity over which a line is
            drawn on the graph
        - npoints: number of points sampled
        - plot: if true, plot results
        - save: if true, saves the plot
    Returns: 'Compare graph to paper'
    """  
    rd.seed(seed)
    testing_function = np.zeros([npoints,2])
    constant_velocity = np.zeros([npoints,2])

    for i in range(len(testing_function)):  
        separation, e, angle, boolean = kicks.post_explosion_params_circular(test_Ai, test_M1, test_M2, test_Mns, kicks.rand_theta(),kicks.rand_phi(),kicks.rand_velocity(test_sigma))
        testing_function[i][0] = separation
        testing_function[i][1] = e
    
    theta = np.linspace(0,3.14,npoints)
    velocity = np.linspace(0,400,npoints)
    
    for j in range(len(constant_velocity)):
        separation, e, angle, boolean = kicks.post_explosion_params_circular(test_Ai, test_M1, test_M2, test_Mns,theta[j],0,sample_velocity)
        constant_velocity[j][0] = separation 
        constant_velocity[j][1] = e
    
    """changing the semi-major axis to period values in days"""
    
    for k in range(len(testing_function)):
        testing_function[k][0] = keplers_third_law(testing_function[k][0],test_M2,test_Mns)
        constant_velocity[k][0] = keplers_third_law(constant_velocity[k][0],test_M2,test_Mns)
    
    if plot:    
        plt.plot(testing_function[:,0], testing_function[:,1], "o")
        plt.xlim(0,50)
        plt.ylim(0,1)
        plt.plot(constant_velocity[:,0], constant_velocity[:,1], 'k-')
        plt.title("post-explosion results")
        plt.xlabel("Period in days")
        plt.ylabel("Eccentricity")
        plt.show()
        plt.close()
    if save:
        plt.plot(testing_function[:,0], testing_function[:,1], "o")
        plt.xlim(0,50)
        plt.ylim(0,1)
        plt.plot(constant_velocity[:,0], constant_velocity[:,1], 'k-')
        plt.title("post-explosion results")
        plt.xlabel("Period in days")
        plt.ylabel("Eccentricity")
        plt.savefig("post_explosion_circular_graph.png")
        plt.close()
    
    rd.seed()
    
    return "True"


def testing_eccentric_function_graph(test_sigma = 100, test_M1 = 5.5, test_M2 = 55, test_Ai = 133, test_Mns = 1.4, seed="Flay",sample_velocity = 100,npoints=10000, plot=True, save =False):
    """Test that the graph of the eccentricity vs the period looks correct
    Arguments:
        - test_sigma: a sample sigma for the rand_velocity function
        - test_M1: solar mass of the first mass pre-explosion
        - test_M2: solar mass of the second mass
        - test_Ai: the initial semi-major axis of the system
        - test_Mns: solar mass of the first mass post-explosion
        - seed: the seed used for the random number generator
        - sample_velocity: a constant velocity over which a line is
            drawn on the graph
        - npoints: number of points sampled
        - plot: if true, plot results
        - save: if true, saves the plot
    Returns: 'Compare graph to paper'
    """  
    rd.seed(seed)
    testing_function = np.zeros([npoints,2])
    constant_velocity = np.zeros([npoints,2])

    for i in range(len(testing_function)):  
        separation, e, boolean = kicks.post_explosion_params_general(test_Ai,\
        test_M1, test_M2, test_Mns, kicks.rand_phi(), kicks.rand_velocity(test_sigma),kicks.rand_true_anomaly(0),0)
        testing_function[i][0] = separation
        testing_function[i][1] = e
    
    theta = np.linspace(0,3.14,npoints)
    velocity = np.linspace(0,400,npoints)
    
    for j in range(len(constant_velocity)):
        separation, e, boolean = kicks.post_explosion_params_general(test_Ai, test_M1, test_M2, test_Mns,theta[j],kicks.rand_velocity(test_sigma),0,0)
        constant_velocity[j][0] = separation 
        constant_velocity[j][1] = e
    
    """changing the semi-major axis to period values in days"""
    
    for k in range(len(testing_function)):
        testing_function[k][0] = keplers_third_law(testing_function[k][0],test_M2,test_Mns)
        constant_velocity[k][0] = keplers_third_law(constant_velocity[k][0],test_M2,test_Mns)
    
    if plot:    
        plt.plot(testing_function[:,0], testing_function[:,1], "o")
        plt.xlim(0,50)
        plt.ylim(0,1)
        plt.plot(constant_velocity[:,0], constant_velocity[:,1], 'k-')
        plt.title("post-explosion results")
        plt.xlabel("Period in days")
        plt.ylabel("Eccentricity")
        plt.show()
        plt.close()
    if save:
        plt.plot(testing_function[:,0], testing_function[:,1], "o")
        plt.xlim(0,50)
        plt.ylim(0,1)
        plt.plot(constant_velocity[:,0], constant_velocity[:,1], 'k-')
        plt.title("post-explosion results")
        plt.xlabel("Period in days")
        plt.ylabel("Eccentricity")
        plt.savefig("post_explosion_circular_graph.png")
        plt.close()
    
    rd.seed()
    
    return "True"

def testing_eccentric_function_momentum(Ai=133, M1=5.5, M2=55, Mns=1.4, test_sigma=100, num_sample=10000, seed = "Lela", tolerance=1e-3):
    """Test that the post_explosion_params_general function produces
    a correct angular momentum against a calculated angular momentum using an 
    eccentricityof zero. This angular momentum is calculated by first finding 
    the anglar velocity, and the individual relative velocities in the center 
    of mass frame pre-explosion. The components of the kick velocity were then 
    added to the velocity of mass 1, which is the super nova mass. The 
    separation in the final center of mass frame accounting for the loss of 
    mass in mass 1 is crossed with the velocities to get the components of 
    angular momentum. The total angular momentum is calculated by taking the 
    absolute value of these components. 
    Arguments:
        - Ai: the initial semi-major axis of the system
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
        phi = kicks.rand_true_anomaly(0)
        
        #getting values from the post_explosion_params_circular function
        separation, e, boolean = kicks.post_explosion_params_general(Ai, M1, M2, Mns,theta,Vk,phi,0)
        
        #calculating the momentum using the results from the function
        Momentum_function = Mns*Msun*M2*Msun*np.sqrt(cgrav*separation*Rsun*(1-e**2)/((Mns+M2)*Msun))
        
        #Calculating the momentum without using the results of the function
        
        #establishing angular velocity 
        omega = np.sqrt(cgrav*Msun*(M1+M2)/(Ai*Rsun)**3) #rad/second
        
        #velocities of the masses before the kick 
        V1_initial = M2/(M1+M2)*omega*Ai*Rsun #cm/second
        V2 = M1/(M1+M2)*omega*Ai*Rsun #cm/second,-y direction
        
        #velocities after the kick, V2 is unaffected by the kick
        V1x_final = Vk*1e5*np.sin(phi)*np.cos(theta)
        V1y_final = V1_initial + Vk*1e5*np.cos(theta)
        V1z_final = Vk*1e5*np.sin(phi)*np.sin(theta)
        
        #separations from the center of mass in the center of mass frame post-explosion
        R2 = Ai*Rsun*Mns/(Mns+M2) #cm
        R1 = Ai*Rsun - R2 #cm
        
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
        


