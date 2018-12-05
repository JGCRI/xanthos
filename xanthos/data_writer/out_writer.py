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
        self.output_names = settings.output_vars
        self.outputs = [pd.DataFrame(all_outputs[out_name]) for out_name in self.output_names]

        self.proj_name = settings.ProjectName
        self.out_folder = settings.OutputFolder
        self.out_format = settings.OutputFormat
        self.out_unit_str = settings.OutputUnitStr
        self.output_in_year = settings.OutputInYear
        self.ChStorageNameStr = settings.OutputNameStr

        self.start_year = settings.StartYear
        self.end_year = settings.EndYear

        try:
            SO = self.outputs[self.output_names.index('soilmoisture')]
            self.SO = np.copy(SO)
        except ValueError:
            self.SO = None

        if self.out_format == FORMAT_NETCDF:
            logging.debug("Save in NetCDF files")
        else:
            logging.debug("Save in CSV files")

        if self.output_in_year:
            logging.debug("Outputting data annually")
            self.outputs = [self.agg_to_year(df) for df in self.outputs]

        if settings.OutputUnit == UNIT_KM3_MTH:
            conversion = grid_areas / 1e6  # mm -> km3
            self.outputs = [df.multiply(conversion, axis=0) for df in self.outputs]

        map(self.write_data, self.output_names, self.outputs)

        # if settings.OutputInYear:
        #     ny = int(settings.EndYear - settings.StartYear + 1)
        #     pet = np.zeros((settings.ncell, ny), dtype=float)
        #     aet = np.zeros((settings.ncell, ny), dtype=float)
        #     q = np.zeros((settings.ncell, ny), dtype=float)
        #     sav = np.zeros((settings.ncell, ny), dtype=float)
        #     ac = np.zeros((settings.ncell, ny), dtype=float)
        #
        #     for i in range(ny):
        #         pet[:, i] = np.sum(self.PET[:, i * 12:(i + 1) * 12], axis=1)
        #         aet[:, i] = np.sum(self.AET[:, i * 12:(i + 1) * 12], axis=1)
        #         q[:, i] = np.sum(self.Q[:, i * 12:(i + 1) * 12], axis=1)
        #         sav[:, i] = np.sum(self.SAV[:, i * 12:(i + 1) * 12], axis=1)
        #         ac[:, i] = np.sum(self.Avg_ChFlow[:, i * 12:(i + 1) * 12], axis=1)
        #
        #     PET = np.copy(pet)
        #     AET = np.copy(aet)
        #     Q = np.copy(q)
        #     SAV = np.copy(sav)
        #     Avg_ChFlow = np.copy(ac)

        # if settings.OutputUnit == UNIT_KM3_MTH:  # convert the original unit mm/month to new unit km3/month
        #     conversion = self.grid_areas / 1e6  # mm -> km3
        #
        #     for j in range(Q.shape[1]):
        #         PET[:, j] *= conversion
        #         AET[:, j] *= conversion
        #         Q[:, j] *= conversion
        #         SAV[:, j] *= conversion
        #         Avg_ChFlow[:, j] = Avg_ChFlow[:, j] * conversion
        #
        logging.debug("Unit is {}".format(settings.OutputUnitStr))

        logging.debug("Output dimension is {}".format(self.outputs[0].shape))

        # self.write_data(settings, 'pet', PET, self.out_format)
        # self.write_data(settings, 'aet', AET, self.out_format)
        # self.write_data(settings, 'q', Q, self.out_format)
        # self.write_data(settings, 'soilmoisture', SAV, self.out_format)
        # self.write_data(settings, 'avgchflow', Avg_ChFlow, self.out_format)

        # if settings.HistFlag == 'True':
        #     logging.info("The following two files are saved as initialization data sets (latest month) for future mode:")
        #     logging.info("ChStorage: monthly output, unit is m^3, dimension is {}".format(self.ChStorage.shape))
        #     settings.OutputNameStr = self.ChStorageNameStr
        #
        #     logging.info("Soil column moisture: monthly output, unit is mm/month, dimension is {}".format(SO.shape))
        #     settings.OutputNameStr = self.ChStorageNameStr

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
            self.save_mat(filename, data)
        elif self.out_format == FORMAT_PARQUET:
            # not yet implemented
            pass

    def save_mat(self, filename, data, varstr):
        """Save output data in the .mat format."""
        filename = filename + ".mat"
        spio.savemat(filename, {varstr: data})

    def save_csv(self, filename, df):
        """Write numpy array as a csv."""
        # add in index as basin or grid cell number
        df.insert(loc=0, column='id', value=df.index.copy()+1)

        filename += ".csv"
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
        filename = filename + ".nc"
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

    def agg_to_year(self, df):
        """Aggregate an array (cells x months) to (cells x years)."""
        return df.groupby(np.arange(len(df.columns)) // NMONTHS, axis=1).sum()
