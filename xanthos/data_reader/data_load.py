"""
Module to load input data files.

Created 8/8/2016

Modified:
@Date 12/21/2018
@author: Xinya Li (xinya.li@pnl.gov), Chris R. Vernon (chris.vernon@pnnl.gov)
@Project: Xanthos V2.0


License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2018, Battelle Memorial Institute
"""

import os
import logging
import numpy as np
from scipy import io as sio
from xanthos.utils.math import sub2ind


class ValidationException(Exception):
    """Exception for invalid configuration options."""

    def __init__(self, *args, **kwargs):
        """Exception for invalid configuration options."""
        Exception.__init__(self, *args, **kwargs)


class DataLoader:
    """Load system-wide input data."""

    def __init__(self, config_obj):
        """
        Validate the configuration set up.

        :param config_obj:      A ConfigObj object that has data from a config file.
        """
        self.s = config_obj

        # ====================================================
        # |        Load data for all configurations          |
        # ====================================================

        # Area value for each land grid cell: 67420 x 1, convert from ha to km2
        self.area = self.load_data(self.s.Area) * 0.01

        # Coordinates for flattened grid:  67420 x 5, the columns are ID#, lon, lat, ilon, ilat
        self.coords = self.load_data(self.s.Coord)

        # Basin ID Map: 67420 x 1, 235 Basins
        self.basin_ids = self.load_data(self.s.BasinIDs, 1).astype(int)

        # Corresponding to Basin ID Map, 235 Basin Names: 1D String Array
        self.basin_names = self.load_data(self.s.BasinNames)

        # GCAM Region ID Map :  67420 x 1 (The nonag region table will be the 'primary' region assignment)
        self.region_ids = self.load_data(self.s.GCAMRegionIDs, 1).astype(int)

        # Corresponding to GCAM Region ID Map
        self.region_names = self.get_region_names()

        # Country ID Map : 67420 x 1 (249 countries: 1-249)
        self.country_ids = self.load_data(self.s.CountryIDs, 1).astype(int)

        # Corresponding to Country ID Map, 0-248 index number and 249 Country Names: 2D String Array
        self.country_names = self.get_country_names()

        self.latitude = np.copy(self.coords[:, 2])
        self.lat_radians = np.radians(self.latitude)

        # ====================================================
        # |        Load data for selected PET module         |
        # ====================================================
        if self.s.pet_module == 'hargreaves':

            # average monthly temperature in degree C
            self.temp = self.load_to_array(self.s.TemperatureFile, var_name=self.s.TempVarName)

            # monthly average of daily temperature range in degree c
            self.dtr = self.load_to_array(self.s.DailyTemperatureRangeFile, var_name=self.s.DTRVarName,
                                          neg_to_zero=True)

        elif self.s.pet_module == 'hs':

            self.hs_tas = self.load_to_array(self.s.hs_tas)
            self.hs_tmin = self.load_to_array(self.s.hs_tmin)
            self.hs_tmax = self.load_to_array(self.s.hs_tmax)

        elif self.s.pet_module == 'pm':
            # values from literature
            ETpara = np.genfromtxt(self.s.pm_params, delimiter=',')
            self.cL = ETpara[:, 0]
            self.beta = ETpara[:, 1]
            self.rslimit = ETpara[:, 2]

            # correlation coefficient for calculating emissivity (range from 0.34 to 0.44), not the ABCD parameter a
            self.ae = ETpara[:, 3]

            # correlation coefficient for calculating emissivity (range from -0.14 to -0.25), not the ABCD parameter b
            self.be = ETpara[:, 4]
            self.Tminopen = ETpara[:, 5]
            self.Tminclose = ETpara[:, 6]
            self.VPDclose = ETpara[:, 7]
            self.VPDopen = ETpara[:, 8]
            self.RBLmin = ETpara[:, 9]
            self.RBLmax = ETpara[:, 10]
            self.rc = ETpara[:, 11]
            self.emiss = ETpara[:, 12]

            # 2-d, rows:11 land cover types, cols:12 months
            self.alpha = np.genfromtxt(self.s.pm_alpha, delimiter=',')
            self.lai = np.genfromtxt(self.s.pm_lai, delimiter=',')
            self.laimax = np.genfromtxt(self.s.pm_laimax, delimiter=',')
            self.laimin = np.genfromtxt(self.s.pm_laimin, delimiter=',')

            # convert missing values to 0
            self.tair_load = self.load_to_array(self.s.pm_tas, 'pm_tas', nan_to_num=True)
            self.TMIN_load = self.load_to_array(self.s.pm_tmin, 'pm_tmin', nan_to_num=True)
            self.rhs_load = self.load_to_array(self.s.pm_rhs, 'pm_rhs', nan_to_num=True)
            self.wind_load = self.load_to_array(self.s.pm_wind, 'pm_wind', nan_to_num=True)
            self.rsds_load = self.load_to_array(self.s.pm_rsds, 'pm_rsds', nan_to_num=True)
            self.rlds_load = self.load_to_array(self.s.pm_rlds, 'pm_rlds', nan_to_num=True)

            # set previous air temp value leaving the first value at 0
            self.tairprev_load = np.zeros_like(self.tair_load)
            self.tairprev_load[1:, :] = self.tair_load[:-1, :]

            # use land cover for each target year (nan to 0);  1-d:67420 cells, 2-d: land cover type, 3-d: years
            self.lct_load = np.nan_to_num(np.load(self.s.pm_lct))

            # static data
            self.elev = np.nan_to_num(np.load(self.s.pm_elev))

        elif self.s.pet_module == 'thornthwaite':
            self.tair = self.load_to_array(self.s.trn_tas, 'trn_tas', nan_to_num=True, warn_nan=True)

        elif self.s.pet_module == 'none':
            # load user supplied PET data
            self.pet_out = self.load_to_array(self.s.pet_file)

        # ====================================================
        # |      Load data for selected runoff module        |
        # ====================================================
        if self.s.runoff_module == 'gwam':
            # monthly precipitation mm/mth
            self.precip = self.load_to_array(self.s.PrecipitationFile, var_name=self.s.PrecipVarName)

            # Max Soil Moisture Map (mm/month): 67420 x 1
            self.max_soil_moist = self.load_data(self.s.max_soil_moisture, 1)

            # Water Bodies: assign MSM = 999, 306 x 2, Col 1 is the cell number in 67420
            self.lakes_msm = self.load_data(self.s.lakes_msm).astype(int)
            self.lakes_msm[:, 0] -= 1

            # Additional water bodies: assign MSM = 999, 421 x 2,  Col 1 is the cell number in 67420
            self.addit_water_msm = self.load_data(self.s.addit_water_msm).astype(int)
            self.addit_water_msm[:, 0] -= 1

            # create a matrix (MSMC: Maximum Soil Moisture Capacity) with all data
            # 1: area; 2: region; 3: Max Soil Moisture (mm/month)
            msmc = self.load_soil_moisture()

            # harmonized grid area
            self.grid_area = np.copy(msmc[:, 0])

            # maximum soil moisture
            self.soil_moisture = np.copy(msmc[:, 2])

            # harmonized grid area
            self.grid_area = np.copy(msmc[:, 0])

            # maximum soil moisture
            self.soil_moisture = np.copy(msmc[:, 2])

            # load soil moisture file if running future
            if self.s.HistFlag.lower() == "true":
                self.sm_prev = 0.5 * self.soil_moisture
            else:
                self.sm_prev = self.load_soil_data()

        elif self.s.runoff_module == 'abcd':
            # monthly precipitation mm/mth
            self.precip = self.load_to_array(self.s.PrecipitationFile, var_name=self.s.PrecipVarName, warn_nan=True)

            # monthly average minimum daily temperature degree C (optional)
            if self.s.TempMinFile is None:
                logging.info(
                    'TempMinFile variable not found for the ABCD runoff module; Snowmelt will not be accounted for.')
                self.tmin = None
            else:
                self.tmin = self.load_to_array(self.s.TempMinFile, var_name=self.s.TempMinVarName,
                                               nan_to_num=True, warn_nan=True)

        # ====================================================
        # |     Load data for selected routing module        |
        # ====================================================
        if self.s.routing_module == 'mrtm':
            map_index = sub2ind([self.s.ngridrow, self.s.ngridcol],
                                self.coords[:, 4].astype(int) - 1,
                                self.coords[:, 3].astype(int) - 1)
            self.flow_dist = self.load_routing_data(
                self.s.flow_distance, self.s.ngridrow, self.s.ngridcol, map_index, rep_val=1000)
            self.flow_dir = self.load_routing_data(self.s.flow_direction, self.s.ngridrow, self.s.ngridcol, map_index)
            self.str_velocity = self.load_routing_data(
                self.s.strm_veloc, self.s.ngridrow, self.s.ngridcol, map_index, rep_val=0)

            self.instream_flow = np.zeros((self.s.ncell,), dtype=float)
            self.chs_prev = self.load_chs_data()

        # ====================================================
        # |     Load data for diagnostics                    |
        # ====================================================
        if self.s.PerformDiagnostics:
            self.vic = self.load_data(self.s.VICDataFile, 0, "q")
            self.unh = self.load_data(self.s.UNHDataFile, 0, "q")
            self.wbmd = self.load_data(self.s.WBMDataFile, 0, 'q')
            self.wbmc = self.load_data(self.s.WBMCDataFile, 0, "q")

        if self.s.calibrate:
            # use basin-level flow as target for calibration; select only columns for basin number and runoff
            self.cal_obs = self.load_data(self.s.cal_observed, 0)[:, [0, 3]]

    def load_soil_data(self):
        """Load soil moisture file into array if in future mode, else stage zeros array."""
        try:
            # Initialize channel storage/soil moisture.
            if self.s.HistFlag.lower() == "true":
                return np.zeros((self.s.ncell,), dtype=float)

            # For future runs, initialize with the last value of the historical channel storage/soil moisture
            else:
                return self.load_data(self.s.SavFile, 0, self.s.SavVarName)[:, -1]

        # if not in use
        except AttributeError:
            return np.zeros((self.s.ncell,), dtype=float)

    def load_soil_moisture(self, missing=-9999):
        """
        Load soil moisture data.

        Assign max soil moisture (mm/month) [2] to Sm.  For historic data use 0.5 * sm to an initial value to pass to
        runoff model. If future mode, read values from historical file.
        """
        data = np.zeros((self.s.ncell, 5), order='F')

        data[:, 0] = self.area
        data[:, 1] = self.region_ids
        data[:, 2] = self.max_soil_moist

        # add max value (999) where water is
        data[self.lakes_msm[:, 0], 2] = self.lakes_msm[:, 1]
        data[self.addit_water_msm[:, 0], 2] = self.addit_water_msm[:, 1]

        country = self.country_ids[:]
        basin = self.basin_ids[:]

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

    def get_country_names(self):
        """Get an array of country names corresponding to GCAM countries."""
        with open(self.s.CountryNames, 'r') as f:
            country = f.read().splitlines()
            return np.array([i.split(',') for i in country])[:, 1]

    def get_region_names(self):
        """Get an array of region names corresponding to the GCAM region id map."""
        with open(self.s.GCAMRegionNames, 'r') as f:
            f.readline()
            region = f.read().split('\n')
            return np.array([i.split(',') for i in region])[:, 0]

    def load_to_array(self, f, var_name=None, neg_to_zero=False, nan_to_num=False, warn_nan=False):
        """
        Load and validate monthly input data.

        Dimension: 67420 x number of years*12, for example:
        Historical: 1950-2005  672 months
        Future: 2006-2100  1140 months

        @:param f:              file path with extension
        @:param var_name:       NetCDF variable name
        @:param neg_to_zero:    convert negative values to zero
        @:param nan_to_num:     convert nan to zero and inf to finite numbers
        @:param warn_nan:       warn if input data contains nan values

        @:return:               array
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
        """
        Check array size of input and check to make sure the total number of months can be split into years.

        :param arr:             input array
        :param text:            name of target variable
        """
        err = "Error: Inconsistent {0} data grid size. Expecting size: {1}. Received size: {2}"

        if not arr.shape[0] == self.s.ncell:
            raise ValidationException(err.format(text, self.s.ncell, arr.shape[0]))

        if not arr.shape[1] == self.s.nmonths:
            raise ValidationException(err.format(text, self.s.nmonths, arr.shape[1]))

        return arr

    @staticmethod
    def load_data(fn, header_num=0, key=None):
        """
        Load data from a file.

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
            raise RuntimeError("File {} has unrecognized extension".format(fn))

        return data

    def load_routing_data(self, fn, ngridrow, ngridcol, map_index, skip=68, rep_val=None):
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
        fd = self.load_data(fn)
        v = self.vectorize(fd, ngridrow, ngridcol, map_index, skip=skip)

        if rep_val is None:
            return v

        else:
            v[np.where(v < rep_val)[0]] = rep_val
            return v

    @staticmethod
    def vectorize(data, ngridrow, ngridcol, map_index, skip):
        """Convert 2D Map (360 x 720) Matrix to 1D Map(67420)."""
        new = np.zeros((ngridrow, ngridcol), dtype=float) - 9999

        for i in range(0, data.shape[0]):
            new[i + skip, :] = data[data.shape[0] - 1 - i, :]

        new = new.reshape((ngridrow * ngridcol,), order='F')

        return new[map_index]

    def load_chs_data(self):
        """Load channel velocity file into array if in future mode, else stage zeros array."""
        try:
            # Initialize channel storage/soil moisture.
            if self.s.HistFlag == "True":
                return np.zeros((self.s.ncell,), dtype=float)

            # For future runs, initialize with the last value of the historical channel storage/soil moisture
            else:
                return self.load_data(self.s.ChStorageFile, 0, self.s.ChStorageVarName)[:, -1]
        except AttributeError:
            return np.zeros((self.s.ncell,), dtype=float)
