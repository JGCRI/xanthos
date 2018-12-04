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
FORMAT_PARQUET = 2

UNIT_MM_MTH = 0
UNIT_KM3_MTH = 1

NMONTHS = 12


class OutWriter:
    """Write out main Xanthos output variables."""

    def __init__(self, settings, grid_areas, outputs):
        self.PET, self.AET, self.Q, self.SAV, self.ChStorage, self.Avg_ChFlow = outputs
        self.settings = settings
        self.grid_areas = grid_areas

        ChStorageNameStr = settings.OutputNameStr
        SO = np.copy(self.SAV)

        flag = settings.OutputFormat
        if flag == FORMAT_NETCDF:
            logging.debug("Save in NetCDF files")
        else:
            logging.debug("Save in CSV files")

        if settings.OutputInYear:
            ny = int(settings.EndYear - settings.StartYear + 1)
            pet = np.zeros((settings.ncell, ny), dtype=float)
            aet = np.zeros((settings.ncell, ny), dtype=float)
            q = np.zeros((settings.ncell, ny), dtype=float)
            sav = np.zeros((settings.ncell, ny), dtype=float)
            ac = np.zeros((settings.ncell, ny), dtype=float)

            for i in range(ny):
                pet[:, i] = np.sum(self.PET[:, i * 12:(i + 1) * 12], axis=1)
                aet[:, i] = np.sum(self.AET[:, i * 12:(i + 1) * 12], axis=1)
                q[:, i] = np.sum(self.Q[:, i * 12:(i + 1) * 12], axis=1)
                sav[:, i] = np.sum(self.SAV[:, i * 12:(i + 1) * 12], axis=1)
                ac[:, i] = np.sum(self.Avg_ChFlow[:, i * 12:(i + 1) * 12], axis=1)

            PET = np.copy(pet)
            AET = np.copy(aet)
            Q = np.copy(q)
            SAV = np.copy(sav)
            Avg_ChFlow = np.copy(ac)

            logging.debug("Output data annually")

        if settings.OutputUnit == UNIT_KM3_MTH:  # convert the original unit mm/month to new unit km3/month
            conversion = self.grid_areas / 1e6  # mm -> km3

            for j in range(Q.shape[1]):
                PET[:, j] *= conversion
                AET[:, j] *= conversion
                Q[:, j] *= conversion
                SAV[:, j] *= conversion
                Avg_ChFlow[:, j] = Avg_ChFlow[:, j] * conversion

            if settings.OutputInYear == 1:
                settings.OutputUnitStr = "km3peryear"
            else:
                settings.OutputUnitStr = "km3permonth"
        else:
            if settings.OutputInYear == 1:
                settings.OutputUnitStr = "mmperyear"
            else:
                settings.OutputUnitStr = "mmpermonth"

        logging.debug("Unit is {}".format(settings.OutputUnitStr))

        logging.debug("Output dimension is {}".format(Q.shape))

        self.SaveData(settings, 'pet', PET, flag)
        self.SaveData(settings, 'aet', AET, flag)
        self.SaveData(settings, 'q', Q, flag)
        self.SaveData(settings, 'soilmoisture', SAV, flag)
        self.SaveData(settings, 'avgchflow', Avg_ChFlow, flag)

        if settings.HistFlag == 'True':
            logging.info("The following two files are saved as initialization data sets (latest month) for future mode:")
            logging.info("ChStorage: monthly output, unit is m^3, dimension is {}".format(self.ChStorage.shape))
            settings.OutputNameStr = ChStorageNameStr

            logging.info("Soil column moisture: monthly output, unit is mm/month, dimension is {}".format(SO.shape))
            settings.OutputNameStr = ChStorageNameStr

    def SaveData(self, settings, var, data, flag):
        """Save output data as a NetCDF or .csv."""
        if var == 'avgchflow':
            unit = 'm3persec'
        else:
            unit = settings.OutputUnitStr

        filename = '{}_{}_{}'.format(var, unit, '_'.join(settings.ProjectName.split(' ')))
        filename = os.path.join(settings.OutputFolder, filename)

        if flag == 0:
            self.SaveNetCDF(filename, data, settings, var)
        else:
            self.SaveCSV(filename, data, settings)

    def SaveMAT(self, filename, data, varstr):
        """Save output data in the .mat format."""
        filename = filename + ".mat"
        spio.savemat(filename, {varstr: data})

    def SaveCSV(self, filename, data, settings):
        """Write numpy array as a csv."""
        # convert to data frame to set header and basin number in file
        df = pd.DataFrame(data)

        # add in index as basin or grid cell number
        df.insert(loc=0, column='id', value=df.index.copy()+1)

        filename += ".csv"

        if settings.OutputInYear == 1:
            cols = ','.join(['{}'.format(i) for i in range(settings.StartYear, settings.EndYear + 1, 1)])
        else:
            col_list = []
            for i in range(settings.StartYear, settings.EndYear + 1, 1):
                for m in range(1, 13):
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

    def SaveNetCDF(self, filename, data, settings, varstr):
        """Write numpy array as a NetCDF."""
        filename = filename + ".nc"
        # open
        datagrp = spio.netcdf.netcdf_file(filename, 'w')
        (nrows, ncols) = data.shape

        # dimensions
        datagrp.createDimension('index', nrows)

        if settings.OutputInYear:
            datagrp.createDimension('year', ncols)
            griddata = datagrp.createVariable('data', 'f4', ('index', 'year'))
        else:
            datagrp.createDimension('month', ncols)
            griddata = datagrp.createVariable('data', 'f4', ('index', 'month'))

        # variables
        unit = settings.OutputUnitStr
        griddata.units = unit
        griddata.description = varstr + "_" + unit

        # data
        griddata[:, :] = data[:, :].copy()

        # close
        datagrp.close()

    def writecsvMap(self, filename, data, settings):
        """Write .csv map."""
        years = list(map(str, list(range(settings.StartYear, settings.EndYear + 1))))
        headerline = "id," + ",".join([year for year in years])

        with open(filename + '.csv', 'w') as outfile:
            np.savetxt(outfile, data, delimiter=',', header=headerline, fmt='%s', comments='')

    def agg_to_year(self, df):
        """Aggregate an array (cells x months) to (cells x years)."""
        return df.groupby(np.arange(len(df.columns)) // NMONTHS, axis=1).sum()
