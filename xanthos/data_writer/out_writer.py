"""
Created on Oct 11, 2016
@author: lixi729
@Project: Xanthos V1.0


License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute


Output Settings:
OutputFormat:  = 0(default, netcdf file); = 1(csv file)
OutputUnit:    = 0(default, mm); = 1(km3)
OutputInYear:  = 0(default, per month); =1(per year, the output will combine 12-month results into annual result)
"""

import os
import numpy as np
import pandas as pd
from scipy import io as spio


def OUTWriter(Settings, area, PET, AET, Q, SAV, ChStorage, Avg_ChFlow):

    ChStorageNameStr = Settings.OutputNameStr
    SO = np.copy(SAV)

    flag = Settings.OutputFormat
    if flag == 0:
        print("Save in NetCDF files")
    else:
        print("Save in CSV files")

    if Settings.OutputInYear == 1:
        ny = int(Settings.EndYear - Settings.StartYear + 1)
        pet = np.zeros((Settings.ncell, ny), dtype=float)
        aet = np.zeros((Settings.ncell, ny), dtype=float)
        q = np.zeros((Settings.ncell, ny), dtype=float)
        sav = np.zeros((Settings.ncell, ny), dtype=float)
        ac = np.zeros((Settings.ncell, ny), dtype=float)

        for i in range(ny):
            pet[:, i] = np.sum(PET[:, i * 12:(i + 1) * 12], axis=1)
            aet[:, i] = np.sum(AET[:, i * 12:(i + 1) * 12], axis=1)
            q[:, i] = np.sum(Q[:, i * 12:(i + 1) * 12], axis=1)
            sav[:, i] = np.sum(SAV[:, i * 12:(i + 1) * 12], axis=1)
            ac[:, i] = np.sum(Avg_ChFlow[:, i * 12:(i + 1) * 12], axis=1)

        del PET, AET, Q, SAV, Avg_ChFlow

        PET = np.copy(pet)
        AET = np.copy(aet)
        Q = np.copy(q)
        SAV = np.copy(sav)
        Avg_ChFlow = np.copy(ac)

        del pet, aet, q, sav, ac

        print("Output data annually")

    if Settings.OutputUnit == 1:  # convert the original unit mm/month to new unit km3/month
        conversion = area / 1e6  # mm -> km3

        for j in range(PET.shape[1]):
            PET[:, j] *= conversion
            AET[:, j] *= conversion
            Q[:, j] *= conversion
            SAV[:, j] *= conversion
            #Avg_ChFlow[:, j] = Avg_ChFlow[:, j] * conversion

        if Settings.OutputInYear == 1:
            Settings.OutputUnitStr = "km3peryear"
        else:
            Settings.OutputUnitStr = "km3permonth"
    else:
        if Settings.OutputInYear == 1:
            Settings.OutputUnitStr = "mmperyear"
        else:
            Settings.OutputUnitStr = "mmpermonth"

    print("Unit is {}".format(Settings.OutputUnitStr))

    print("Output dimension is {}".format(PET.shape))

    SaveData(Settings, 'pet', PET, flag)
    SaveData(Settings, 'aet', AET, flag)
    SaveData(Settings, 'q', Q, flag)
    SaveData(Settings, 'soilmoisture', SAV, flag)
    SaveData(Settings, 'avgchflow', Avg_ChFlow, flag)

    if Settings.HistFlag == 'True':
        print("The following two files are saved as initialization data sets (latest month) for future mode:")
        print("ChStorage: monthly output, unit is m^3, dimension is {}".format(ChStorage.shape))
        Settings.OutputNameStr = ChStorageNameStr

        print("Soil column moisture: monthly output, unit is mm/month, dimension is {}".format(SO.shape))
        Settings.OutputNameStr = ChStorageNameStr

    return Q, Avg_ChFlow


def SaveData(settings, var, data, flag):

    if var == 'avgchflow':
        unit = 'm3persec'
    else:
        unit = settings.OutputUnitStr

    filename = os.path.join(settings.OutputFolder, '{}_{}_{}'.format(var, unit, '_'.join(settings.ProjectName.split(' '))))

    if flag == 0:
        SaveNetCDF(filename, data, settings, var)
    else:
        SaveCSV(filename, data, settings)


def SaveMAT(filename, data, varstr):
    filename = filename + ".mat"
    spio.savemat(filename, {varstr: data})


def SaveCSV(filename, data, settings):

    # convert to data frame to set header and basin number in file
    df = pd.DataFrame(data)

    # add in index as basin or grid cell number
    df.insert(loc=0, column='id', value=df.index.copy()+1)

    filename += ".csv"

    if settings.OutputInYear == 1:
        cols = ','.join(['{}'.format(i) for i in range(settings.StartYear, settings.EndYear + 1, 1)])
    else:
        l = []
        for i in range(settings.StartYear, settings.EndYear + 1, 1):
            for m in range(1, 13):
                if m < 10:
                    mth = '0{}'.format(m)
                else:
                    mth = m
                l.append('{}{}'.format(i, mth))
        cols = ','.join(l)

    # set header
    hdr = 'id,{}'.format(cols)

    try:
        df.columns = hdr.split(',')
    except ValueError:
        raise

    df.to_csv(filename, index=False)


def SaveNetCDF(filename, data, Settings, varstr):
    filename = filename + ".nc"
    # open
    datagrp = spio.netcdf.netcdf_file(filename, 'w')
    (nrows, ncols) = data.shape

    # dimensions
    datagrp.createDimension('index', nrows)

    if Settings.OutputInYear:
        datagrp.createDimension('year', ncols)
        griddata = datagrp.createVariable('data', 'f4', ('index', 'year'))
    else:
        datagrp.createDimension('month', ncols)
        griddata = datagrp.createVariable('data', 'f4', ('index', 'month'))

    # variables
    unit = Settings.OutputUnitStr
    griddata.units = unit
    griddata.description = varstr + "_" + unit

    # data
    griddata[:, :] = data[:, :].copy()

    # close
    datagrp.close()


def writecsvMap(filename, data, Settings):
    years = map(str, range(Settings.StartYear, Settings.EndYear + 1))
    headerline = "id," + ",".join([year for year in years])

    with open(filename + '.csv', 'w') as outfile:
        np.savetxt(outfile, data, delimiter=',', header=headerline, fmt='%s', comments='')