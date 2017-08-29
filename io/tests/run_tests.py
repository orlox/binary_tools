#!/usr/bin/env python
from binary_tools.io.tests import test_mesa_data

__author__ = "Pablo Marchant"
__credits__ = ["Pablo Marchant"]
__license__ = "GPL"
__version__ = "3.0"
__maintainer__ = "Pablo Marchant"
__email__ = "pamarca@gmail.com"

"""Runs all tests in this package
"""

def run_tests():
    """Runs all the tests and saves the graphs
    """

    print("Run test_history")
    result1 = test_mesa_data.test_history("history.data")
    if result1:
        print("test_history worked")
    else:
        print("test_history FAILED!!!!!")
            
    print("Run test_profile")
    result2 = test_mesa_data.test_profile("profile1.data")
    if result2:
        print("test_profile worked")
    else:
        print("test_profile FAILED!!!!!")

    return result1 and result2

run_tests()
