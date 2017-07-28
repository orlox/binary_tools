#!/usr/bin/env python
from binary_tools.utils.tests import test_kicks

__author__ = "Kaliroe Pappas"
__credits__ = ["Kaliroe Pappas"]
__license__ = "GPL"
__version__ = "3.0"
__maintainer__ = "Kaliroe Pappas"
__email__ = "kmpappa2@illinois.edu"

"""Runs all tests in this package
"""

def run_tests():
    """Runs all the tests and saves the graphs
    -using a sigma value of 100
    Returns True
    """

    print("Run test_rand_phi")
    result1 = test_kicks.test_rand_phi()
    if result1:
        print("test_rand_phi worked")
    else:
        print("test_rand_phi FAILED!!!!!")
            
    print("Run test_rand_theta")
    result2 = test_kicks.test_rand_theta()
    if result2:
        print("test_rand_theta worked")
    else:
        print("test_rand_theta FAILED!!!!!")
        
    print("Run test_rand_velocity")
    result3 = test_kicks.test_rand_velocity(100)
    if result3:
        print("test_rand_velocity worked")
    else:
        print("test_rand_velocity FAILED!!!!!")

    print("Run testing_circular_function_graph")
    result4 = test_kicks.testing_circular_function_graph()
    if result4:
        print("testing_circular_function_graph worked")
    else:
        print("testing_circular_function_graph FAILED!!!!!")
    return True

    print("Run testing_circular_function_momentum")
    result4 = test_kicks.testing_circular_function_momentum()
    if result4:
        print("testing_circular_function_momentum worked")
    else:
        print("testing_circular_function_graph FAILED!!!!!")
    return True

    print("Run testing_eccentric_function_graph")
    result4 = test_kicks.testing_eccentric_function_graph()
    if result4:
        print("testing_eccentric_function_graph worked")
    else:
        print("testing_eccentric_function_graph FAILED!!!!!")
    return True

    print("Run testing_eccentric_function_momentum")
    result4 = test_kicks.testing_eccentric_function_momentum()
    if result4:
        print("testing_eccentric_function_momentum worked")
    else:
        print("testing_eccentric_function_momentum FAILED!!!!!")
    return True

    print("Run testing_inverse_kick")
    result5 = test_kicks.testing_inverse_kick()
    if result5:
        print("testing_inverse_kick worked")
    else:
        print("testing_inverse_kick FAILED!!!!!")
    return True

    print("Run testing_momentum_full_eccentric")
    result6 = test_kicks.testing_momentum_full_eccentric()
    if result6:
        print("testing_momentum_full_eccentric worked")
    else:
        print("testing_momentum_full_eccentric FAILED!!!!!")
    return True





