"""
Module to write output data files.

Output settings:
OutputFormat:  = 0(default, netcdf file); = 1(csv file)
OutputUnit:    = 0(default, mm); = 1(km3)
OutputInYear:  = 0(default, per month); = 1(per year, the output will combine 12-month results into annual result)

Created on Oct 11, 2016

@author: lixi729
@Project: Xanthos V1.0

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute
"""

import os
import logging
import numpy as np
import pandas as pd
from scipy import io as spio

FORMAT_NETCDF = 0
FORMAT_CSV = 1
FORMAT_MAT = 2
FORMAT_PARQUET = 3

UNIT_MM_MTH = 0
UNIT_KM3_MTH = 1

NMONTHS = 12


class OutWriter:
    """Write out main Xanthos output variables."""

    def __init__(self, settings, grid_areas, all_outputs):
        """
        Format and call appropriate writer for output variables.

        :param settings:        parsed settings from input configuration file
        :param grid_areas:      map of basin indices to grid cell area, in km2 (numpy array)
        :param all_outputs:     dictionary mapping all output names (strings) to their values (numpy arrays)
        """
        self.output_names = settings.output_vars
        self.outputs = [pd.DataFrame(all_outputs[out_name]) for out_name in self.output_names]

        # output options
        self.proj_name = settings.ProjectName
        self.out_folder = settings.OutputFolder
        self.out_format = settings.OutputFormat
        self.out_unit_str = settings.OutputUnitStr
        self.output_in_year = settings.OutputInYear
        self.ChStorageNameStr = settings.OutputNameStr

        self.start_year = settings.StartYear
        self.end_year = settings.EndYear

        # parameter checks
        if self.out_format not in [FORMAT_NETCDF, FORMAT_CSV, FORMAT_MAT, FORMAT_PARQUET]:
            logging.warning("Output format {} is invalid; writing output as .csv".format(self.out_format))
            self.out_format = FORMAT_CSV

        try:
            SO = self.outputs[self.output_names.index('soilmoisture')]
            self.SO = np.copy(SO)
        except ValueError:
            self.SO = None

        if self.output_in_year:
            logging.debug("Outputting data annually")
            self.outputs = [self.agg_to_year(df) for df in self.outputs]

        if settings.OutputUnit == UNIT_KM3_MTH:
            conversion = grid_areas / 1e6  # mm -> km3
            self.outputs = [df.multiply(conversion, axis=0) for df in self.outputs]

        logging.debug("Unit is {}".format(settings.OutputUnitStr))
        logging.debug("Output dimension is {}".format(self.outputs[0].shape))

        for var, data in zip(self.output_names, self.outputs):
            self.write_data(var, data)

    def write_data(self, var, data):
        """Save output data as a NetCDF or .csv, .mat, or parquet file."""
        if var == 'avgchflow':
            unit = 'm3persec'
        else:
            unit = self.out_unit_str

        filename = '{}_{}_{}'.format(var, unit, '_'.join(self.proj_name.split(' ')))
        filename = os.path.join(self.out_folder, filename)

        if self.out_format == FORMAT_NETCDF:
            self.save_netcdf(filename, data, var)

        elif self.out_format == FORMAT_CSV:
            self.save_csv(filename, data)

        elif self.out_format == FORMAT_MAT:
            self.save_mat(filename, data, var)

        elif self.out_format == FORMAT_PARQUET:
            self.save_parquet(filename, data, var)

    def save_mat(self, filename, data, varstr):
        """Write output data in the .mat format."""
        filename += ".mat"

        spio.savemat(filename, {varstr: data})

    def save_csv(self, filename, df):
        """Write pandas DataFrame as a .csv file."""
        filename += ".csv"

        # add in index as basin or grid cell number
        df.insert(loc=0, column='id', value=df.index.copy() + 1)

        year_range = range(self.start_year, self.end_year + 1, 1)

        if self.output_in_year:
            cols = ','.join(['{}'.format(i) for i in year_range])
        else:
            col_list = []
            for i in year_range:
                for m in range(1, NMONTHS + 1):
                    if m < 10:
                        mth = '0{}'.format(m)
                    else:
                        mth = m
                    col_list.append('{}{}'.format(i, mth))
            cols = ','.join(col_list)

        # set header
        hdr = 'id,{}'.format(cols)

        try:
            df.columns = hdr.split(',')
        except ValueError:
            raise

        df.to_csv(filename, index=False)

    def save_netcdf(self, filename, data, varstr):
        """Write numpy array as a NetCDF."""
        filename += ".nc"

        # open
        datagrp = spio.netcdf.netcdf_file(filename, 'w')
        (nrows, ncols) = data.shape

        # dimensions
        datagrp.createDimension('index', nrows)

        if self.output_in_year:
            datagrp.createDimension('year', ncols)
            griddata = datagrp.createVariable('data', 'f4', ('index', 'year'))
        else:
            datagrp.createDimension('month', ncols)
            griddata = datagrp.createVariable('data', 'f4', ('index', 'month'))

        # variables
        unit = self.out_unit_str
        griddata.units = unit
        griddata.description = varstr + "_" + unit

        # data
        griddata[:, :] = data[:, :].copy()

        # close
        datagrp.close()

    def save_parquet(self, filename, df, varstr):
        """Write pandas DataFrame to parquet file."""
        filename += ".parquet"

    def agg_to_year(self, df):
        """Aggregate an array(cells x months) to(cells x years)."""
        return df.groupby(np.arange(len(df.columns)) // NMONTHS, axis=1).sum()
