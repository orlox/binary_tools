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

print("Run test_rand_phi")
result = test_kicks.test_rand_phi()
if result:
    print("test_rand_phi worked")
else:
    print("test_rand_phi FAILED!!!!!")
