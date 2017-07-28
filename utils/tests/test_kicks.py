#!/usr/bin/env python
from binary_tools.utils.constants import *
from binary_tools.utils import kicks
from binary_tools.utils.Keplers_laws import *
from binary_tools.utils.orbits import *
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
        plt.hist(phi_array, bins=np.linspace(0,2*np.pi,nbins), alpha = 0.5, label = "function output")
        plt.hist(test_array, bins=np.linspace(0,2*np.pi,nbins), alpha = 0.5, label = "expected value")
        plt.title("phi distribution")
        plt.xlabel("phi value")
        plt.ylabel("distribution")
        plt.legend(loc='upper right')
        plt.show()
        plt.close()
    if save:
        plt.hist(phi_array, bins=np.linspace(0,2*np.pi,nbins), alpha = 0.5, label = "function output")
        plt.hist(test_array, bins=np.linspace(0,2*np.pi,nbins), alpha = 0.5, label = "expected value")
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
        plt.hist(theta_array, bins=np.linspace(0,np.pi,nbins), alpha = 0.5, label = "function output")
        plt.hist(test_array, bins=np.linspace(0,np.pi,nbins), alpha = 0.5, label = "expected value")
        plt.title("theta distribution")
        plt.xlabel("theta value")
        plt.ylabel("distribution")
        plt.legend(loc='upper right')
        plt.show()
        plt.close()
    if save:
        plt.hist(theta_array, bins=np.linspace(0,np.pi,nbins), alpha = 0.5, label = "function output")
        plt.hist(test_array, bins=np.linspace(0,np.pi,nbins), alpha = 0.5, label = "expected value")
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
        plt.hist(velocity_array, bins=np.linspace(0,3*sigma,nbins), alpha = 0.5, label = "function output")
        plt.hist(test_array, bins=np.linspace(0,3*sigma,nbins), alpha = 0.5, label = "expected value")
        plt.title("velocity distribution")
        plt.xlabel("velocity value")
        plt.ylabel("distribution")
        plt.legend(loc='upper right')
        plt.show()
        plt.close()
    if save:
        plt.hist(velocity_array, bins=np.linspace(0,3*sigma,nbins), alpha = 0.5, label = "function output")
        plt.hist(test_array, bins=np.linspace(0,3*sigma,nbins), alpha = 0.5, label = "expected value")
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
        plt.hist(true_anomaly_array, bins=np.linspace(0,2*np.pi,nbins), alpha = 0.5, label = "function output")
        plt.hist(test_array, bins=np.linspace(0,2*np.pi,nbins), alpha = 0.5, label = "expected value")
        plt.title("true anomaly distribution")
        plt.xlabel("true anomaly value")
        plt.ylabel("distribution")
        plt.legend(loc='upper right')
        plt.show()
        plt.close()
    if save:
        plt.hist(true_anomaly_array, bins=np.linspace(0,2*np.pi,nbins), alpha = 0.5, label = "function output")
        plt.hist(test_array, bins=np.linspace(0,2*np.pi,nbins), alpha = 0.5, label = "expected value")
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
        semi_major, e, angle, boolean = kicks.post_explosion_params_circular(Ai, M1, M2, Mns,theta,phi,Vk)
        
        #calculating the momentum using the results from the function
        Momentum_function = angular_momentum(semi_major,Mns,M2,e)
        
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
        semi_major, e, angle, boolean = kicks.post_explosion_params_circular(test_Ai, test_M1, test_M2, test_Mns, kicks.rand_theta(),kicks.rand_phi(),kicks.rand_velocity(test_sigma))
        testing_function[i][0] = semi_major
        testing_function[i][1] = e
    
    theta = np.linspace(0,3.14,npoints)
    velocity = np.linspace(0,400,npoints)
    
    for j in range(len(constant_velocity)):
        semi_major, e, angle, boolean = kicks.post_explosion_params_circular(test_Ai, test_M1, test_M2, test_Mns,theta[j],0,sample_velocity)
        constant_velocity[j][0] = semi_major 
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


def testing_eccentric_function_graph(test_sigma = 100, test_M1 = 5.5, test_M2 = 55, test_Ai = 133, test_Mns = 1.4, seed="David Berne",sample_velocity = 100,npoints=10000, plot=True, save =False):
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
        semi_major, e, boolean = kicks.post_explosion_params_general(test_Ai,\
        test_M1, test_M2, test_Mns, 0, kicks.rand_theta(), kicks.rand_phi(),\
        kicks.rand_velocity(test_sigma),kicks.rand_true_anomaly(0))
        testing_function[i][0] = semi_major
        testing_function[i][1] = e
    
    theta = np.linspace(0,3.14,npoints)
    velocity = np.linspace(0,400,npoints)
    
    for j in range(len(constant_velocity)):
        semi_major, e, boolean = kicks.post_explosion_params_general(test_Ai, test_M1, test_M2, test_Mns,0,theta[j],0,sample_velocity,0)
        constant_velocity[j][0] = semi_major 
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

def testing_eccentric_function_momentum(Ai=133, M1=5.5, M2=55, Mns=1.4, test_sigma=100, num_sample=10000, seed = "Clara", tolerance=1e-3):
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
        true_anomaly = kicks.rand_true_anomaly(0)
        phi = kicks.rand_phi()
        
        #getting values from the post_explosion_params_general function
        semi_major, e, boolean = kicks.post_explosion_params_general(Ai, M1, M2, Mns,0,theta,phi,Vk,true_anomaly)
        
        #calculating the momentum using the results from the function
        Momentum_function = angular_momentum(semi_major,Mns,M2,e)
        
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
            print(Vk, theta, phi, true_anomaly, Momentum_calculated, Momentum_function, i)
            return False
            
    rd.seed()
        
    return True
        

def testing_eccentric_kick(Ai=133, M1=5.5, M2=55, Mns=1.4, num_sample=100, seed = "Guarnaschelli",tolerance=1e-4):
    """Test for the posta-explosion_params_general that calulates
    the necessary kick at perigee and appgee to cause a circular 
    orbital, then plug that result into the function to ensure it
    returns an eccentricity of zero.
    Arguments:
    - Ai: the initial semi-major axis of the system
    - M1: solar mass of the first mass pre-explosion
    - M2: solar mass of the second mass
    - Mns: solar mass of the first mass post-explosion
    - num_sample: number of points sampled
    - seed: the seed used for the random number generator
    - tolerance: tolerance for the test
    Returns: True or False as to whether the test was successful
    """  
    rd.seed(seed)
    M1 = M1*Msun
    M2 = M2*Msun
    Mns = Mns*Msun
    Ai = Ai*Rsun
        
    
    e_samples = np.linspace(0,.99,num_sample)
    
    for e in e_samples:

        V_apogee = np.sqrt(cgrav*(M1+M2)*(1-e)/(Ai*(1+e)))        
        V_perigee = np.sqrt(cgrav*(M1+M2)*(1+e)/(Ai*(1-e)))

        V_circular_apogee = np.sqrt(cgrav*(Mns+M2)/(Ai*(1+e)))
        V_circular_perigee = np.sqrt(cgrav*(Mns+M2)/(Ai*(1-e)))        
        
        V_kick_apogee = np.absolute((V_apogee - V_circular_apogee)*1e-5)
        V_kick_perigee = np.absolute((V_circular_perigee - V_perigee)*1e-5)
        
        theta_apogee = np.pi
        theta_perigee = np.pi
        
        if V_circular_apogee > V_apogee:
            theta_apogee = 0
        if V_circular_perigee > V_perigee:
            theta_perigee = 0
        
        semi_major_a, e_a, boulean_a = kicks.post_explosion_params_general(Ai/Rsun,M1/Msun,M2/Msun,Mns/Msun,e,theta_apogee,0,V_kick_apogee,np.pi)
        semi_major_p, e_p, boulean_p = kicks.post_explosion_params_general(Ai/Rsun,M1/Msun,M2/Msun,Mns/Msun,e,theta_perigee,0,V_kick_perigee,0)
        
        
        if e_a > tolerance or e_p > tolerance:
            return False
        
    rd.seed()
        
    return True


def testing_inverse_kick(Ai=133, M1=5.5, M2=55, Mns=1.4, test_sigma=1000, num_sample=100, seed="Tamlin",tolerance=1e-4):
    """Test for the post_explosions_params_general function that kicks
    a circular system with mass loss, then reverses that with mass 
    gain back into a circular orbit. There are four different possible
    ways to reverse the kick that are dependent on the true anomaly 
    and theta. 1) the inital kick sends the masses into an eccentric
    orbit in the same initial direction, with a true anomaly between
    0 and pi. 2) the inital kick sends the masses into an eccentric
    orbit in the same initial direction, with a true anomaly between
    pi and 2pi. 3)the inital kick sends the masses into an eccentric
    orbit in the opposit direction, with a true anomaly between
    0 and pi. 4) the inital kick sends the masses into an eccentric
    orbit in the opposit direction, with a true anomaly between
    pi and 2pi.
        
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
        theta = kicks.rand_theta()
        V_kick = kicks.rand_velocity(test_sigma)
        
        semi_major_i, e_i, boulean_i = kicks.post_explosion_params_general(Ai,M1,M2,Mns,0,theta,0,V_kick,0)
        k = semi_major_i*(1-e_i**2)/(e_i*Ai) - 1/e_i
        
        true_anomaly = np.arccos(k)
        
        semi_major_f, e_f, boulean_f = kicks.post_explosion_params_general(semi_major_i,Mns,M2,M1,e_i,np.pi-theta,np.pi,V_kick,true_anomaly)
        
        Worked = True
        
        if e_f > tolerance:
            true_anomaly = 2*np.pi - true_anomaly
            semi_major_f, e_f, boulean_f = kicks.post_explosion_params_general(semi_major_i,Mns,M2,M1,e_i,np.pi-theta,np.pi,V_kick,true_anomaly)
            if e_f > tolerance:
                semi_major_f, e_f, boulean_f = kicks.post_explosion_params_general(semi_major_i,Mns,M2,M1,e_i,theta,np.pi,V_kick,true_anomaly)
                if e_f > tolerance:
                    true_anomaly = 2*np.pi - true_anomaly
                    semi_major_f, e_f, boulean_f = kicks.post_explosion_params_general(semi_major_i,Mns,M2,M1,e_i,theta,np.pi,V_kick,true_anomaly)
                    if e_f > tolerance:
                        Worked = False       
                    
    rd.seed()
        
    return Worked
        
      

def testing_momentum_full_eccentric(Ai=133, M1=5.5, M2=55, Mns=1.4, test_sigma=15, num_sample=100, seed="Lucien",tolerance=1e-4):
    """Test that the post_explosion_params_general function produces
    a correct angular momentum against a calculated angular momentum. This 
    angular momentum is calculated by first finding the velocity of M1 pre-
    explosion, then adding to that the components of the kick to get a velocity
    in the theta, phi, and radial direction. This velocity is then changed to 
    x and y components. Next the center of mass possition and velocity are 
    calculated, and using those values the relative velocities and postion are
    calculated. From there, the angular momentum is calculated using the cross-
    product method, and adding the resulting values.
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
    e_samples = np.linspace(0,.99,num_sample)
    
    for e in e_samples:
        
        #establishing random parameters
        Vk = kicks.rand_velocity(test_sigma)*1e5
        theta = kicks.rand_theta()
        phi = kicks.rand_phi()
        true_anomaly = kicks.rand_true_anomaly(e) 
        separation_i = Rsun*separation(Ai,e,true_anomaly)
    
        #getting values from the post_explosion_params_general function
        semi_major_f, e_f, boolean = kicks.post_explosion_params_general(Ai, M1, M2, Mns,e,theta,phi,Vk*1e-5,true_anomaly)
        
        #calculating the momentum using the results from the function
        Momentum_function = angular_momentum(semi_major_f,Mns,M2,e_f)
    
        V_theta_i = np.sqrt(cgrav*(M1+M2)*Msun*Rsun*Ai*(1-e**2))/separation_i
        V_radius_i = np.sqrt(cgrav*(M1+M2)*Msun*(2/separation_i-1/(Rsun*Ai)-Ai*Rsun*(1-e**2)/(separation_i**2)))
        
        V_radius = V_radius_i + Vk*np.sin(theta)*np.cos(phi)
        V_theta = V_theta_i + Vk*np.cos(theta)
        V_phi = Vk*np.sin(theta)*np.sin(phi)
        
        V1x = V_radius
        V1y = np.sqrt(V_theta**2+V_phi**2)
        
        
        R_cm = Mns*separation_i/(Mns+M2) #x direction
        V_cm_x = Mns*V1x/(Mns+M2) #x dirrection
        V_cm_y = Mns*V1y/(Mns+M2) #y dirrection
        
        Vx1_prime = V1x - V_cm_x
        Vy1_prime = V1y - V_cm_y
        
        Rx1_prime = separation_i - R_cm #+x direction
        Ry1_prime = 0 
        
        Vx2_prime = -V_cm_x
        Vy2_prime = -V_cm_y
        
        Rx2_prime = 0 - R_cm #-x direction
        Ry2_prime = 0
        
        momentum1x = Vx1_prime*Mns*Msun
        momentum1y = Vy1_prime*Mns*Msun
        
        momentum2x = Vx2_prime*M2*Msun
        momentum2y = Vy2_prime*M2*Msun
        
        angular1 = Rx1_prime*momentum1y - Ry1_prime*momentum1x #z direction
        angular2 = Rx2_prime*momentum2y - Ry2_prime*momentum2x #z direction
        
        Momentum_calculated = (angular1 + angular2)
        
        if abs(Momentum_function - Momentum_calculated)/Momentum_function>tolerance:
            print(Vk, theta, phi, true_anomaly, Momentum_calculated, Momentum_function, e)
            return False


    rd.seed()
        
    return True



  
