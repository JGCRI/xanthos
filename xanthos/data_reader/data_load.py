"""
Module to load input data files
Created 8/8/2016

Modified:
@Date 10/05/2016
@author: Xinya Li (xinya.li@pnl.gov), Chris R. Vernon (chris.vernon@pnnl.gov)
@Project: Xanthos V2.0


License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute
"""

import os
import sys
from scipy import io as spio
import numpy as np

from xanthos.utils.numpy_parser import GetArrayCSV, GetArrayTXT


class LoadReferenceData:
    """
    Load reference data.

    :param settings:        settings object from configuration
    """
    def __init__(self, settings):

        # Area value for each land grid cell: 67420 x 1, convert from ha to km2
        self.area = load_const_griddata(settings.Area) * 0.01

        # Coordinates for flattened grid:  67420 x 5, the columns are ID#, lon, lat, ilon, ilat
        self.coords = load_const_griddata(settings.Coord)

        # Basin ID Map: 67420 x 1, 235 Basins
        self.basin_ids = load_const_griddata(settings.BasinIDs, 1).astype(int)

        # Corresponding to Basin ID Map, 235 Basin Names: 1D String Array
        self.basin_names = load_const_griddata(settings.BasinNames)

        # GCAM Region ID Map :  67420 x 1 (The nonag region table will be the 'primary' region assignment)
        self.region_ids = load_const_griddata(settings.GCAMRegionIDs, 1).astype(int)

        # Corresponding to GCAM Region ID Map
        with open(settings.GCAMRegionNames, 'r') as f:
            f.readline()
            temp = f.read().split('\n')
            self.region_names = np.array([i.split(',') for i in temp])[:, 0]

        # Country ID Map : 67420 x 1 (249 countries: 1-249)
        self.country_ids = load_const_griddata(settings.CountryIDs, 1).astype(int)

        # Corresponding to Country ID Map, 0-248 index number and 249 Country Names: 2D String Array
        with open(settings.CountryNames, 'r') as f:
            temp = f.read().splitlines()
            self.country_names = np.array([i.split(',') for i in temp])[:, 1]

        if settings.runoff_module == 'gwam':
            # Max Soil Moisture Map (mm/month): 67420 x 1
            self.max_soil_moist = load_const_griddata(settings.MaxSoilMois, 1)

            # Water Bodies: assign MSM = 999, 306 x 2, Col 1 is the cell number in 67420
            self.lakes_msm = load_const_griddata(settings.LakesMSM).astype(int)
            self.lakes_msm[:, 0] -= 1

            #    ''' Rivers:  assign MSM = 999, 4198 x 2, Col 1 is the cell number in 67420
            #    constants['RiversMSM']      = load_const_griddata(settings.RiversMSM).astype(int)
            #    constants['RiversMSM'][:,0]  -= 1

            # Additional water bodies: assign MSM = 999, 421 x 2,  Col 1 is the cell number in 67420
            self.addit_water_msm = load_const_griddata(settings.AdditWaterMSM).astype(int)
            self.addit_water_msm[:, 0] -= 1


def load_climate_data(fle, var_name, n_cells, n_months, neg_to_zero=False):
    """
    Loads and checks input climate data.

    Dimension: 67420 x number of years*12, for example:
    Historical: 1950-2005  672 months
    Future: 2006-2100  1140 months

    @:param fle:            file path with extension
    @:param var_name:       NetCDF variable name
    @:param neg_to_zero:    convert negative values to zero
    @:param n_cells:        number of cells
    @:param n_months:       number of months

    @:return:               array
    """
    a = load_const_griddata(fle, 0, var_name)

    if neg_to_zero:
        a[np.where(a < 0)] = 0

    return check_climate_data(a, n_cells=n_cells, n_months=n_months, text=var_name)


def load_routing_data(fle, ngridrow, ngridcol, map_index, skip=68, rep_val=None):
    """
    Load routing data.

    DRT data, 280 x 720, -9999 for missing values, convert to 67420 X 1

    @:param fle             file to load
    @:param ngridrow        number of grids per row
    @:param ngridcol        number of grids per column
    @:param map_index
    @:param skip
    @:param rep_val         value to replace with when less than value
    """
    fd = load_const_griddata(fle)
    v = vectorize(fd, ngridrow, ngridcol, map_index, skip=skip)

    if rep_val is None:
        return v

    else:
        v[np.where(v < rep_val)[0]] = rep_val
        return v


def load_soil_data(settings):
    """
    Load soil moisture file into array if in future mode, else stage zeros array.
    """
    try:
        # Initialize channel storage/soil moisture.
        if settings.HistFlag == "True":
            return np.zeros((settings.ncell,), dtype=float)

        # For future runs, initialize with the last value of the historical channel storage/soil moisture
        else:
            return load_const_griddata(settings.SavFile, 0, settings.SavVarName)[:, -1]

    # if not in use
    except AttributeError:
        return np.zeros((settings.ncell,), dtype=float)


def load_chs_data(settings):
    """
    Load channel velocity file into array if in future mode, else stage zeros array.
    """
    try:

        # Initialize channel storage/soil moisture.
        if settings.HistFlag == "True":
            return np.zeros((settings.ncell,), dtype=float)

        # For future runs, initialize with the last value of the historical channel storage/soil moisture
        else:
            return load_const_griddata(settings.ChStorageFile, 0, settings.ChStorageVarName)[:, -1]
    except AttributeError:
        return np.zeros((settings.ncell,), dtype=float)


def load_gcm_var(fn, varname):
    """
    Loads climate data from the specified GCM
    """

    if not os.path.isfile(fn):
        print("Error: File does not exist:  {}".format(fn))
        sys.exit()

    temp = spio.loadmat(fn)
    data = temp[varname]

    return data


def check_climate_data(data, n_cells, n_months, text):
    """
    Check array size of input and check to make sure the total number of months can be split into years.

    :param data:            input array
    :param n_cells:         number of cells
    :param n_months:        number of months
    :param text:            name of target variable
    """
    err_cell = "Error: Inconsistent {0} data grid size. Expecting size: {1}. Received size: {2}".format(text, n_cells, data.shape[0])
    err_mth = "Error: Inconsistent {0} data grid size. Expecting size: {1}. Received size: {2}".format(text, n_months, data.shape[1])

    if not data.shape[0] == n_cells:
        print(err_cell)
        sys.exit(1)

    if not data.shape[1] == n_months:
        print(err_mth)
        sys.exit(1)

    if not data.shape[1] % 12 == 0:
        print("Error: Number of months in climate data can not be converted into integral years.")
        sys.exit(1)

    return data


def load_const_griddata(fn, headerNum=0, key=" "):
    """
    Load constant grid data stored in files defined in GRID_CONSTANTS.
    """

    # for MATLAB files
    if fn.endswith('.mat'):
        data = load_gcm_var(fn, key)

    # for Numpy pickled files
    elif fn.endswith('.npy'):
        data = np.load(fn)

    # for text files
    elif fn.endswith('.txt'):

        if not os.path.isfile(fn):
            print("Error: File does not exist:", fn)
            sys.exit()

        try:
            data = GetArrayTXT(fn, headerNum)

        except:
            with open(fn, 'r') as f:
                data = np.array(f.read().splitlines())

    # for CSV files
    elif fn.endswith('.csv'):

        if not os.path.isfile(fn):
            print("Error: File does not exist:", fn)
            sys.exit()

        data = GetArrayCSV(fn, headerNum)

    # for NetCDF classic files
    elif fn.endswith('.nc'):

        if not os.path.isfile(fn):
            print("Error: File does not exist:", fn)
            sys.exit()

        datagrp = spio.netcdf.netcdf_file(fn, 'r')

        # copy() added to handle numpy 'ValueError:assignment destination is read-only' related to non-contiguous memory
        try:
            data = datagrp.variables[key][:, :].copy()

        except:
            data = datagrp.variables[key][:].copy()

        datagrp.close()

    return data


def vectorize(data, ngridrow, ngridcol, map_index, skip):
    """
    Convert 2D Map (360 x 720) Matrix to 1D Map(67420)
    """
    new = np.zeros((ngridrow, ngridcol), dtype=float) - 9999

    for i in range(0, data.shape[0]):
        new[i + skip, :] = data[data.shape[0] - 1 - i, :]

    new = new.reshape((ngridrow * ngridcol,), order='F')

    return new[map_index]


def load_soil_moisture(d, ngrids, missing=-9999):
    data = np.zeros((ngrids, 5), order='F')

    data[:, 0] = d.area
    data[:, 1] = d.region_ids
    data[:, 2] = d.max_soil_moist

    # add max value (999) where water is
    data[d.lakes_msm[:, 0], 2] = d.lakes_msm[:, 1]
    data[d.addit_water_msm[:, 0], 2] = d.addit_water_msm[:, 1]

    country = d.country_ids[:]
    basin = d.basin_ids[:]

    # Ignore all the cells in which we are missing an ID value for soil moisture, country, or basin.
    # Thus, country and basin coverage must be consistent.
    # Basins coverage is smaller, and GCAM region ignores Greenland.
    invalid = np.where((data[:, 2] == 0) | (country == 0) | (basin == 0))[0]

    # should this be 0:2
    data[invalid, 1:2] = 0

    # should these be returned?
    country[invalid] = missing
    basin[invalid] = missing

    return data
