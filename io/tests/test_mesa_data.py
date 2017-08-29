#!/usr/bin/env python
from binary_tools.io import mesa_data as md
import matplotlib.pyplot as plt
import filecmp
import numpy as np

__author__ = "Pablo Marchant"
__credits__ = ["Pablo Marchant"]
__license__ = "GPL"
__version__ = "3.0"
__maintainer__ = "Pablo Marchant"
__email__ = "pamarca@gmail.com"

"""Test the functions in the mesa_data.py module
"""

def test_history(filename, plot=False, save=True):
    """Test that a MESA history.data file can be loaded properly
    Arguments:
        - filename: path to the history file
        - plot: if true, plot results and show them
        - save: if true, saves the plot
    Returns: True if the test is succesful, False otherwise
    """

    #load a full original history.data
    data = md.mesa_data(filename)
    data.save_as_ascii("test.data")
    data.save_as_hdf5("test.hdf5")
    mn1 = data.get("model_number")
    lcT1 = data.get("log_center_T")
    #load the saved ascii data
    data = md.mesa_data("test.data")
    data.save_as_ascii("test.data2")
    data.save_as_hdf5("test.hdf52")
    mn2 = data.get("model_number")
    lcT2 = data.get("log_center_T")
    #load the saved hdf5 data
    data = md.mesa_data("test.hdf5", is_hdf5=True)
    data.save_as_ascii("test.data3")
    data.save_as_hdf5("test.hdf53")
    mn3 = data.get("model_number")
    lcT3 = data.get("log_center_T")
    #load only some columns of the original history.data
    data = md.mesa_data(filename, read_data_cols = ["model_number","log_Teff","log_L","log_center_T","log_center_Rho"])
    data.save_as_ascii("test.data4")
    data.save_as_hdf5("test.hdf54")
    mn4 = data.get("model_number")
    lcT4 = data.get("log_center_T")
    #load the truncated hdf5
    data = md.mesa_data("test.hdf54", is_hdf5=True)
    data.save_as_ascii("test.data5")
    data.save_as_hdf5("test.hdf55")
    mn5 = data.get("model_number")
    lcT5 = data.get("log_center_T")
    #load the saved truncated hdf5
    data = md.mesa_data("test.hdf55", is_hdf5=True)
    data.save_as_ascii("test.data6")
    data.save_as_hdf5("test.hdf56")
    mn6 = data.get("model_number")
    lcT6 = data.get("log_center_T")
    #and reload it
    data = md.mesa_data("test.hdf56", is_hdf5=True)
    data.save_as_ascii("test.data7")
    data.save_as_hdf5("test.hdf57")
    mn7 = data.get("model_number")
    lcT7 = data.get("log_center_T")
    #and finally, reload it as ascii
    data = md.mesa_data("test.data7")
    data.save_as_ascii("test.data8")
    data.save_as_hdf5("test.hdf58")
    mn8 = data.get("model_number")
    lcT8 = data.get("log_center_T")

    #check if ascii files come out right
    success = filecmp.cmp('test.data', 'test.data2')
    print(success)
    success = success and filecmp.cmp('test.data', 'test.data3')
    print(success)
    success = success and filecmp.cmp('test.data4', 'test.data5')
    print(success)
    success = success and filecmp.cmp('test.data4', 'test.data6')
    print(success)
    success = success and filecmp.cmp('test.data4', 'test.data7')
    print(success)
    success = success and filecmp.cmp('test.data4', 'test.data8')
    print(success)

    success = success and np.array_equal(mn1,mn2)
    print(success)
    success = success and np.array_equal(mn1,mn3)
    print(success)
    success = success and np.array_equal(mn1,mn4)
    print(success)
    success = success and np.array_equal(mn1,mn5)
    print(success)
    success = success and np.array_equal(mn1,mn6)
    print(success)
    success = success and np.array_equal(mn1,mn7)
    print(success)
    success = success and np.array_equal(mn1,mn8)
    print(success)

    success = success and np.array_equal(lcT1,lcT2)
    print(success)
    success = success and np.array_equal(lcT1,lcT3)
    print(success)
    success = success and np.array_equal(lcT1,lcT4)
    print(success)
    success = success and np.array_equal(lcT1,lcT5)
    print(success)
    success = success and np.array_equal(lcT1,lcT6)
    print(success)
    success = success and np.array_equal(lcT1,lcT7)
    print(success)
    success = success and np.array_equal(lcT1,lcT8)
    print(success)

    #plot a logRho_c-logT_c diagram
    plt.plot(data.get("log_center_Rho"),data.get("log_center_T"))
    plt.gca().set_ylabel("log T center")
    plt.gca().set_xlabel("log Rho center")
    plt.show()

    return success


def test_profile(filename, plot=False, save=True):
    """Test that a MESA profile file can be loaded properly
    Arguments:
        - filename: path to the history file
        - plot: if true, plot results
        - save: if true, saves the plot
    Returns: True if the test is succesful, False otherwise
    """

    #repeat for a profile
    data = md.mesa_data(filename)
    data.save_as_ascii("test.data")
    data.save_as_hdf5("test.hdf5")
    zone1 = data.get("zone")
    lT1 = data.get("logT")
    #load the saved ascii data
    data = md.mesa_data("test.data")
    data.save_as_ascii("test.data2")
    data.save_as_hdf5("test.hdf52")
    zone2 = data.get("zone")
    lT2 = data.get("logT")
    #load the saved hdf5 data
    data = md.mesa_data("test.hdf5", is_hdf5=True)
    data.save_as_ascii("test.data3")
    data.save_as_hdf5("test.hdf53")
    zone3 = data.get("zone")
    lT3 = data.get("logT")
    #load only some columns of the original history.data
    data = md.mesa_data(filename, read_data_cols = ["zone","logT","logRho"])
    data.save_as_ascii("test.data4")
    data.save_as_hdf5("test.hdf54")
    zone4 = data.get("zone")
    lT4 = data.get("logT")
    #load the truncated hdf5
    data = md.mesa_data("test.hdf54", is_hdf5=True)
    data.save_as_ascii("test.data5")
    data.save_as_hdf5("test.hdf55")
    zone5 = data.get("zone")
    lT5 = data.get("logT")
    #load the saved truncated hdf5
    data = md.mesa_data("test.hdf55", is_hdf5=True)
    data.save_as_ascii("test.data6")
    data.save_as_hdf5("test.hdf56")
    zone6 = data.get("zone")
    lT6 = data.get("logT")
    #and reload it
    data = md.mesa_data("test.hdf56", is_hdf5=True)
    data.save_as_ascii("test.data7")
    data.save_as_hdf5("test.hdf57")
    zone7 = data.get("zone")
    lT7 = data.get("logT")
    #and finally, reload it as ascii
    data = md.mesa_data("test.data7")
    data.save_as_ascii("test.data8")
    data.save_as_hdf5("test.hdf58")
    zone8 = data.get("zone")
    lT8 = data.get("logT")

    #check if ascii files come out right
    success = filecmp.cmp('test.data', 'test.data2')
    print(success)
    success = success and filecmp.cmp('test.data', 'test.data3')
    print(success)
    success = success and filecmp.cmp('test.data4', 'test.data5')
    print(success)
    success = success and filecmp.cmp('test.data4', 'test.data6')
    print(success)
    success = success and filecmp.cmp('test.data4', 'test.data7')
    print(success)
    success = success and filecmp.cmp('test.data4', 'test.data8')
    print(success)

    success = success and np.array_equal(zone1,zone2)
    print(success)
    success = success and np.array_equal(zone1,zone3)
    print(success)
    success = success and np.array_equal(zone1,zone4)
    print(success)
    success = success and np.array_equal(zone1,zone5)
    print(success)
    success = success and np.array_equal(zone1,zone6)
    print(success)
    success = success and np.array_equal(zone1,zone7)
    print(success)
    success = success and np.array_equal(zone1,zone8)
    print(success)

    success = success and np.array_equal(lT1,lT2)
    print(success)
    success = success and np.array_equal(lT1,lT3)
    print(success)
    success = success and np.array_equal(lT1,lT4)
    print(success)
    success = success and np.array_equal(lT1,lT5)
    print(success)
    success = success and np.array_equal(lT1,lT6)
    print(success)
    success = success and np.array_equal(lT1,lT7)
    print(success)
    success = success and np.array_equal(lT1,lT8)
    print(success)

    #plot a logRho-logT diagram
    plt.plot(data.get("logRho"),data.get("logT"))
    plt.gca().set_ylabel("log T")
    plt.gca().set_xlabel("log Rho")
    plt.show()

    return success
