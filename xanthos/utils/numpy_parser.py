"""
Helper functions for reading in numpy arrays.

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute
"""

import numpy as np


def GetArrayCSV(filename, header_num):

    try:
        with open(filename, "rU") as f:
            return np.genfromtxt(f, delimiter=",", skip_header=header_num, filling_values="0")

    except IOError:
        pass


def GetArrayTXT(filename, header_num):

    try:
        with open(filename, "rU") as f:

            a = np.__version__

            if int(a.split('.')[1]) >= 10:

                return np.genfromtxt(f, delimiter=" ", skip_header=header_num, filling_values="0")

            else:
                return np.genfromtxt(f, delimiter=" ", skiprows=header_num, filling_values="0")

    except IOError:
        pass
