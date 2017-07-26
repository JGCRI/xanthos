"""
License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute
"""

import csv
import numpy


def GetArrayCSV(filename, headerNum):
    try:
        f = open(filename, "rU")
        data = numpy.genfromtxt(f, delimiter=",", skip_header=headerNum, filling_values="0")
        return data
        f.close()
    except IOError:
        pass


def GetArrayTXT(filename, headerNum):
    try:
        f = open(filename, "rU")
        a = numpy.__version__
        if int(a.split('.')[1])>= 10:
            data = numpy.genfromtxt(f, delimiter=" ", skip_header=headerNum, filling_values="0")
        else:
            data = numpy.genfromtxt(f, delimiter=" ", skiprows=headerNum, filling_values="0")
        return data
        f.close()
    except IOError:
        pass