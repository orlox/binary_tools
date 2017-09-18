#!/usr/bin/env python
import numpy as np
from binary_tools.utils.constants import *
from binary_tools.io.mesa_data import mesa_data
from binary_tools.utils.Keplers_laws import *
from binary_tools.utils.kicks import *

__author__ = "Pablo Marchant"
__credits__ = ["Pablo Marchant"]
__license__ = "GPL"
__version__ = "3.0"
__maintainer__ = "Pablo Marchant"
__email__ = "pamarca@gmail.com"

"""Simple population synthesis by post-processing MESA
TODO: properly explain
"""

def evolve(hist1="LOGS1/history.data", hist2="LOGS2/history.data",bhist="binary_history.data",is_hdf5=False,theta=0,phi=0,vkick=0,outfile="output.data"):
    """TODO: explain
    """

    hist1 = mesa_data(hist1, is_hdf5 = is_hdf5)
    hist2 = mesa_data(hist2, is_hdf5 = is_hdf5)
    bhist = mesa_data(bhist, is_hdf5 = is_hdf5)

    #find pre-explosion conditions
    for index, pmi in enumerate(bhist.get("point_mass_index")):
        if pmi>0:
            break

    initial_separation = bhist.get("binary_separation")[index-1]
    initial_m1 = bhist.get("star_1_mass")[index-1]
    initial_m2 = bhist.get("star_2_mass")[index-1]
    initial_period = keplers_third_law(initial_separation*Rsun/1e5, initial_m1, initial_m2)

    post_explosion_params = post_explosion_params_circular(
            initial_separation,
            initial_m1,
            initial_m2,
            1.4,
            theta,
            phi,
            vkick)

    print(initial_period, initial_separation, initial_m1, initial_m2)
    print("####Initial parameters####")
    print("M1(Msun):"+'{0:>5e}'.format(bhist.get("star_1_mass")[0]))
    print("M2(Msun):"+'{0:>5e}'.format(bhist.get("star_2_mass")[0]))
    print("a(Rsun):"+'{0:>5e}'.format(bhist.get("binary_separation")[0]))
    print("P(days):"+'{0:>5e}'.format(bhist.get("period_days")[0]))
    print("####Pre-explosion parameters####")
    print("M1(Msun):"+'{0:>5e}'.format(bhist.get("star_1_mass")[index-1]))
    print("M2(Msun):"+'{0:>5e}'.format(bhist.get("star_2_mass")[index-1]))
    print("a(Rsun):"+'{0:>5e}'.format(bhist.get("binary_separation")[index-1]))
    print("P(days):"+'{0:>5e}'.format(bhist.get("period_days")[index-1]))

    post_exp_m1 = 1.4

    angular_momentum = initial_m2*Msun*post_exp_m1*Msun*np.sqrt(
            cgrav*post_explosion_params[0]*Rsun*(1-post_explosion_params[1]**2)/((initial_m2+post_exp_m1)*Msun))

    m2 = hist2.get("star_mass")[index:]
    m2_dot = hist2.get("star_mdot")[index:]
    age = hist2.get("star_age")[index:]
    radius = np.power(10,hist2.get("log_R")[index:])

    am = np.zeros(len(m2))
    am[0] = angular_momentum
    separation = np.zeros(len(m2))
    separation[0] = post_explosion_params[0]
    eccentricity = np.zeros(len(m2))
    eccentricity[0] = post_explosion_params[1]
    m1 = np.zeros(len(m2))
    m1[0] = post_exp_m1
    m1_dot = np.zeros(len(m2))
    m1_dot[0] = 0
    RL2 = np.zeros(len(m2))
    RL2[0] = eval_roche_lobe(m2[0],m1[0],separation[0])
    period = np.zeros(len(m2))
    period[0] = keplers_third_law(separation[0], m1[0], m2[0])

    print("####Post-explosion parameters####")
    print("M1(Msun):"+'{0:>5e}'.format(post_exp_m1))
    print("M2(Msun):"+'{0:>5e}'.format(bhist.get("star_2_mass")[index-1]))
    print("a(Rsun):"+'{0:>5e}'.format(post_explosion_params[0]))
    print("P(days):"+'{0:>5e}'.format(period[0]))
    print("eccentricity:"+'{0:>5e}'.format(eccentricity[0]))



    f = open(outfile,'w')
    f.write('{0:>20}'.format("step")+
            '{0:>20}'.format("age")+
            '{0:>20}'.format("m1")+
            '{0:>20}'.format("m2")+
            '{0:>20}'.format("m1_dot")+
            '{0:>20}'.format("m2_dot")+
            '{0:>20}'.format("m2")+
            '{0:>20}'.format("separation")+
            '{0:>20}'.format("period")+
            '{0:>20}'.format("eccentricity")+
            '{0:>20}'.format("r_div_rl")+
            '{0:>20}'.format("angular_momentum"))
    f.write('\n')
    f.write('{0:>20}'.format(0)+
            '{0:>20e}'.format(age[0])+
            '{0:>20e}'.format(m1[0])+
            '{0:>20e}'.format(m2[0])+
            '{0:>20e}'.format(m1_dot[0])+
            '{0:>20e}'.format(m2_dot[0])+
            '{0:>20e}'.format(m2[0])+
            '{0:>20e}'.format(separation[0])+
            '{0:>20e}'.format(period[0])+
            '{0:>20e}'.format(eccentricity[0])+
            '{0:>20e}'.format(radius[0]/RL2[0])+
            '{0:>20e}'.format(am[0]))
    f.write('\n')

    if eccentricity[0]>1:
        print("Unbound system")
        f.close()
        return

    for i in range(1,len(m1)):

        dt = (age[i] - age[i-1])*secyear
        
        #compute wind accretion
        alpha_wind = 1.5
        beta_wind = 1.
        vwind = np.sqrt(beta_wind*2*cgrav*m2[i-1]*Msun/(radius[i-1]*Rsun))
        vorb = np.sqrt(cgrav*(m1[i-1]+m2[i-1])*Msun/(separation[i-1]*Rsun))
        vfrac = vorb/vwind
        frac_wind = 1/np.sqrt(1-eccentricity[i-1]**2)\
                *(cgrav*m2[i]*Msun/vwind**2)**2\
                *alpha_wind/(2*separation[i-1]**2*Rsun**2)\
                *1/np.power(1+vfrac**2,1.5)

        m1_dot[i] = -m2_dot[i]*frac_wind
        m1[i] = m1[i-1] + m1_dot[i]*dt/secyear
        am[i] = am[i-1] + (m1[i-1]/(m1[i-2]+m2[i-1])*separation[i-1]*Rsun)**2*2*np.pi/(period[i-1]*3600*24) \
            *np.sqrt(1 - eccentricity[i-1]**2)*dt*(m2_dot[i-1]*Msun/secyear)
        eccentricity[i] = eccentricity[i-1]

        separation[i] =  ((am[i]/(m1[i]*m2[i]*Msun**2))**2) \
                *(m1[i]+m2[i])*Msun / cgrav / (1 - eccentricity[i]**2)
        separation[i] = separation[i]/Rsun

        period[i] = keplers_third_law(separation[i]*Rsun, m1[i], m2[i])
        RL2[i] = eval_roche_lobe(m2[i],m1[i],separation[i])

        #print(i-1, age[i-1], m1[i-1], m2[i-1], separation[i-1], period[i-1], eccentricity[i-1],am[i-1])
        f.write('{0:>20}'.format(i)+
                '{0:>20e}'.format(age[i])+
                '{0:>20e}'.format(m1[i])+
                '{0:>20e}'.format(m2[i])+
                '{0:>20e}'.format(m1_dot[i])+
                '{0:>20e}'.format(m2_dot[i])+
                '{0:>20e}'.format(m2[i])+
                '{0:>20e}'.format(separation[i])+
                '{0:>20e}'.format(period[i])+
                '{0:>20e}'.format(eccentricity[i])+
                '{0:>20e}'.format(radius[i]/RL2[i])+
                '{0:>20e}'.format(am[i]))
        f.write('\n')

    print("####final parameters####")
    print("M1(Msun):"+'{0:>5e}'.format(m1[-1]))
    print("M2(Msun):"+'{0:>5e}'.format(m2[-1]))
    print("a(Rsun):"+'{0:>5e}'.format(separation[-1]))
    print("P(days):"+'{0:>5e}'.format(period[-1]))
    print("eccentricity:"+'{0:>5e}'.format(eccentricity[-1]))
    f.close()


def eval_roche_lobe(m1,m2,a):
    q = np.power(m1/m2,1/3.)
    return a*0.49*q*q/(0.6*q*q + np.log(q))
