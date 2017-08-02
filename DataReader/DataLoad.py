'''
Module to load input data files
Created 8/8/2016

Modified:
@Date 10/05/2016
@author: Xinya Li (xinya.li@pnl.gov)
@Project: Xanthos V1.0


License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute

'''

import os, sys
from scipy import io as spio
import numpy as np
#from netCDF4 import Dataset

from Utils.NumpyParser import GetArrayCSV, GetArrayTXT
from Utils.Math import sub2ind


def load_gcm_data(settings):
    """Loads and checks input climate data (monthly average precipitation,
    monthly average temperature, and average daily temperature range

    Dimension: 67420 x number of years*12, for example:
    Historical: 1950-2005  672 months
    Future: 2006-2100  1140 months"""

    pre = load_const_griddata(settings.PrecipitationFile, 0, settings.PrecipVarName)  # mm/month
    temp = load_const_griddata(settings.TemperatureFile, 0, settings.TempVarName)  # deg C
    dtr = load_const_griddata(settings.DailyTemperatureRangeFile, 0, settings.DTRVarName)  # deg C
    # DTR is defined as the range between maximum and minimum daily temperatures in Celsius.
    # So should be strictly a positive value. Defaults the value to zero.
    dtr[np.where(dtr < 0)] = 0  # NaNs in dtr -> warning RuntimeWarning: invalid value encountered

    # Initialize channel storage/soil moisture.
    if settings.HistFlag == "True":  # For historical runs just use zeros
        inic = np.zeros((settings.ncell,), dtype=float)
        inis = np.zeros((settings.ncell,), dtype=float)
    else:
        # For future runs, we want to initialize with the last value of the historical channel storage/soil moisture.
        inic = load_const_griddata(settings.ChStorageFile, 0, settings.ChStorageVarName)[:, -1]
        inis = load_const_griddata(settings.SavFile, 0, settings.SavVarName)[:, -1]

    check_size((pre, temp, dtr), settings.ncell, settings.nmonths)
    check_n_years(pre)

    return pre, temp, dtr, inis, inic


def load_gcm_var(fn, varname):
    """Loads climate data from the specified GCM"""

    if not os.path.isfile(fn):
        print "Error! File does not exist:", fn
        sys.exit()

    temp = spio.loadmat(fn)
    data = temp[varname]

    return data


def check_size(data, ng, nm):
    n = 0
    m = 0
    var = ['precip', 'temp', 'dtr']
    for d in data:
        if not d.shape[0] == ng:
            print "Error! Inconsistent " + var[m] + " data grid size, which is not", ng
            n += 1
        if not d.shape[1] == nm:
            print "Error! Inconsistent " + var[m] + " data month size, which is not", nm
            n += 1
        m += 1

    if n > 0:
        sys.exit()


def check_n_years(data):
    if not data.shape[1] % 12 == 0:
        print "Error! Number of months in climate data can not be converted into integral years."
        sys.exit()


def load_map_data(settings):
    constants = {}

    '''' Area value for each land grid cell: 67420 x 1, convert from ha to km2 '''
    constants['Area'] = load_const_griddata(settings.Area) * 0.01

    '''' Coordinates for flattened grid:  67420 x 5, the columns are ID#, lon, lat, ilon, ilat '''
    constants['Coord'] = load_const_griddata(settings.Coord)
    settings.mapindex = sub2ind([settings.ngridrow, settings.ngridcol], constants['Coord'][:, 4].astype(int) - 1,
                                constants['Coord'][:, 3].astype(int) - 1)

    ''' DRT data, 280 x 720, -9999 for missing values, convert to 67420 X 1'''
    orig = load_const_griddata(settings.FlowDis)
    constants['FlowDis'] = vectorizeGridData(orig, settings, 68)

    orig = load_const_griddata(settings.FlowDir)
    constants['FlowDir'] = vectorizeGridData(orig, settings, 68)

    ''' Basin ID Map: 67420 x 1, 235 Basins '''
    constants['BasinIDs'] = load_const_griddata(settings.BasinIDs, 1).astype(int)

    ''' Corresponding to Basin ID Map, 235 Basin Names: 1D String Array'''
    constants['BasinNames'] = load_const_griddata(settings.BasinNames)

    ''' GCAM Region ID Map :  67420 x 1 (The nonag region table will be the 'primary' region assignment)'''
    constants['GCAMRegionIDs'] = load_const_griddata(settings.GCAMRegionIDs, 1).astype(int)

    ''' Corresponding to GCAM Region ID Map'''
    with open(settings.GCAMRegionNames, 'r') as f:
        f.readline()
        temp = f.read().split('\n')
        constants['GCAMRegionNames'] = np.array([i.split(',') for i in temp])[:, 0]

    '''Country ID Map : 67420 x 1 (249 countries: 1-249)'''
    constants['CountryIDs'] = load_const_griddata(settings.CountryIDs, 1).astype(int)

    ''' Corresponding to Country ID Map, 0-248 index number and 249 Country Names: 2D String Array'''
    with open(settings.CountryNames, 'r') as f:
        temp = f.read().splitlines()
        constants['CountryNames'] = np.array([i.split(',') for i in temp])[:, 1]

    ''' Max Soil Moisture Map (mm/month): 67420 x 1'''
    constants['MaxSoilMois'] = load_const_griddata(settings.MaxSoilMois, 1)

    ''' Water Bodies: assign MSM = 999, 306 x 2, Col 1 is the cell number in 67420'''
    constants['LakesMSM'] = load_const_griddata(settings.LakesMSM).astype(int)
    constants['LakesMSM'][:, 0] -= 1

    #    ''' Rivers:  assign MSM = 999, 4198 x 2, Col 1 is the cell number in 67420'''
    #    constants['RiversMSM']      = load_const_griddata(settings.RiversMSM).astype(int)
    #    constants['RiversMSM'][:,0]  -= 1

    ''' Additional water bodies: assign MSM = 999, 421 x 2,  Col 1 is the cell number in 67420'''
    constants['AdditWaterMSM'] = load_const_griddata(settings.AdditWaterMSM).astype(int)
    constants['AdditWaterMSM'][:, 0] -= 1

    #     for key in constants:
    #         print key, constants[key].shape, constants[key].dtype

    return constants


def load_const_griddata(fn, headerNum=0, key=" "):
    """ Load constant grid data stored in files defined in GRID_CONSTANTS."""

    if fn.endswith('.mat'):
        data = load_gcm_var(fn, key)

    elif fn.endswith('.txt'):
        if not os.path.isfile(fn):
            print "Error! File does not exist:", fn
            sys.exit()

        try:
            data = GetArrayTXT(fn, headerNum)
        except:
            with open(fn, 'r') as f:
                data = np.array(f.read().splitlines())

    elif fn.endswith('.csv'):
        if not os.path.isfile(fn):
            print "Error! File does not exist:", fn
            sys.exit()

        data = GetArrayCSV(fn, headerNum)

    elif fn.endswith('.nc'):
        if not os.path.isfile(fn):
            print "Error! File does not exist:", fn
            sys.exit()

#        datagrp = Dataset(fn, 'r', format='NETCDF4')
        datagrp = spio.netcdf.netcdf_file(fn, 'r')

        # copy() added to handle numpy 'ValueError:assignment destination is read-only' related to non-contiguous memory
        try:
#            data = datagrp[key][:, :]
            data = datagrp.variables[key][:, :].copy()
        except:
#            data = datagrp[key][:]            
            data = datagrp.variables[key][:].copy()

        datagrp.close()

    return data


def vectorizeGridData(data, settings, skip):
    """Convert 2D Map (360 x 720) Matrix to 1D Map(67420)"""

    new = np.zeros((settings.ngridrow, settings.ngridcol), dtype=float) - 9999
    for i in range(0, data.shape[0]):
        new[i + skip, :] = data[data.shape[0] - 1 - i, :]

    new = new.reshape((settings.ngridrow * settings.ngridcol,), order='F')

    return new[settings.mapindex]


def get_MaxSoilMoisture_matrix(d, ngrids, missing=-9999):
    data = np.zeros((ngrids, 5), order='F')

    data[:, 0] = d['Area']
    data[:, 1] = d['GCAMRegionIDs']
    data[:, 2] = d['MaxSoilMois']
    data[d['LakesMSM'][:, 0], 2] = d['LakesMSM'][:, 1]
    # data[d['RiversMSM '][:,0], 2]    = d['RiversMSM '][:,1]
    data[d['AdditWaterMSM'][:, 0], 2] = d['AdditWaterMSM'][:, 1]

    country = d['CountryIDs'][:]
    basin = d['BasinIDs'][:]

    # Ignore all the cells in which we are missing an ID value for Sm, country, or basin.
    # Thus, country and basin coverage must be consistent.
    # Basins coverage is smaller, and GCAM region ignores Greenland.
    invalid = np.where((data[:, 2] == 0) | (country == 0) | (basin == 0))[0]

    data[invalid, 1:2] = 0
    country[invalid] = missing
    basin[invalid] = missing

    return data
