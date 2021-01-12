import os
import logging

import numpy as np

from scipy import io as sio


class ValidationException(Exception):
    """Exception for invalid configuration options."""

    def __init__(self, *args, **kwargs):
        """Exception for invalid configuration options."""
        Exception.__init__(self, *args, **kwargs)


class DataUtils:
    """Methods to process input data.  Loads the configuration object."""

    def __init__(self, nmonths, ncells=67420):

        self.nmonths = nmonths
        self.ncells = ncells

    def load_to_array(self, f, var_name=None, neg_to_zero=False, nan_to_num=False, warn_nan=False):
        """Load and validate monthly input data.

        Dimension: 67420 x number of years*12, for example:
        Historical: 1950-2005  672 months
        Future: 2006-2100  1140 months

        :param f:              file path with extension
        :param var_name:       NetCDF variable name
        :param neg_to_zero:    convert negative values to zero
        :param nan_to_num:     convert nan to zero and inf to finite numbers
        :param warn_nan:       warn if input data contains nan values

        :return:               array

        """

        # load data to array from file, unless data is already an array
        if isinstance(f, np.ndarray):
            arr = f
            f = 'in memory'
        else:
            arr = self.load_data(f, 0, var_name)

        if var_name is None:
            var_name = os.path.splitext(os.path.basename(f))[0]

        if neg_to_zero:
            arr[np.where(arr < 0)] = 0

        if warn_nan and np.any(np.isnan(arr)):
            logging.warning("NaNs found in input file {}".format(var_name))

        if nan_to_num:
            arr = np.nan_to_num(arr)

        return self.validate(arr, text=var_name)

    def validate(self, arr, text):
        """Check array size of input and check to make sure the total number of months can be split into years.

        :param arr:             input array
        :param text:            name of target variable

        """

        err = "Error: Inconsistent {0} data grid size. Expecting size: {1}. Received size: {2}"

        if not arr.shape[0] == self.ncells:
            raise ValidationException(err.format(text, self.ncells, arr.shape[0]))

        if not arr.shape[1] == self.nmonths:
            raise ValidationException(err.format(text, self.nmonths, arr.shape[1]))

        return arr

    @staticmethod
    def load_data(fn, header_num=0, key=None):
        """Load data from a file.

        :param fn:              name of file to load
        :param header_num:      number of lines in file to skip, if text or csv file
        :param key:             name of variable to extract, if matlab or NetCDF file

        """

        if not os.path.isfile(fn):
            raise IOError("Error: File does not exist:", fn)

        # for MATLAB files
        if fn.endswith('.mat'):
            data = sio.loadmat(fn)[key]

        # for Numpy pickled files
        elif fn.endswith('.npy'):
            data = np.load(fn)

        # for text files
        elif fn.endswith('.txt'):
            try:
                data = np.genfromtxt(fn, delimiter=" ", skip_header=header_num, filling_values="0")

            except:
                with open(fn, 'r') as f:
                    data = np.array(f.read().splitlines())

        # for CSV files
        elif fn.endswith('.csv'):
            data = np.genfromtxt(fn, delimiter=",", skip_header=header_num, filling_values="0")

        # for NetCDF classic files
        elif fn.endswith('.nc'):
            datagrp = sio.netcdf.netcdf_file(fn, 'r', mmap=False)

            # copy() added to handle numpy 'ValueError:assignment destination is read-only' for non-contiguous memory
            data = datagrp.variables[key][:].copy()

            datagrp.close()

            # we only support data in little-endian format, so convert it if the NetCDF is big-endian
            if data.dtype.byteorder == ">":
                data = data.byteswap().newbyteorder()

        else:
            raise ValidationException("File {} has unrecognized extension".format(fn))

        return data
