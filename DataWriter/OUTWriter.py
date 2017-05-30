'''
Created on Oct 11, 2016
@author: lixi729
@Project: Xanthos V1.0

Output Settings:
OutputFormat:  = 0(default, netcdf file); = 1(csv file)
OutputUnit:    = 0(default, mm); = 1(km3)
OutputInYear:  = 0(default, per month); =1(per year, the output will combine 12-month results into annual result)

'''

import numpy as np
from scipy import io as spio
from netCDF4 import Dataset


def OUTWriter(Settings, area, PET, AET, Q, SAV, ChStorage, Avg_ChFlow):
    ChStorageNameStr = Settings.OutputNameStr
    SO = np.copy(SAV)

    flag = Settings.OutputFormat
    if flag == 0:
        print "Save in netcdf files"
    else:
        print "Save in csv files"

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

        Settings.OutputNameStr = "_".join(Settings.OutputNameStr.split("_")[0:-2]) \
                                 + "_" + str(Settings.StartYear) + "_" + str(Settings.EndYear)

        print "Output data annually"

    if Settings.OutputUnit == 1:  # convert the original unit mm/month to new unit km3/month
        conversion = area / 1e6  # mm -> km3

        for j in range(PET.shape[1]):
            PET[:, j] = PET[:, j] * conversion
            AET[:, j] = AET[:, j] * conversion
            Q[:, j] = Q[:, j] * conversion
            SAV[:, j] = SAV[:, j] * conversion
            Avg_ChFlow[:, j] = Avg_ChFlow[:, j] * conversion

        if Settings.OutputInYear == 1:
            Settings.OutputUnitStr = "km^3/year"
        else:
            Settings.OutputUnitStr = "km^3/month"
    else:
        if Settings.OutputInYear == 1:
            Settings.OutputUnitStr = "mm/year"
        else:
            Settings.OutputUnitStr = "mm/month"

    print "Unit is ", Settings.OutputUnitStr

    print "Output dimension is", PET.shape

    SaveData(Settings, 'pet', PET, flag)
    SaveData(Settings, 'aet', AET, flag)
    SaveData(Settings, 'q', Q, flag)
    SaveData(Settings, 's', SAV, flag)
    SaveData(Settings, 'Avg_ChFlow', Avg_ChFlow, flag)

    if Settings.HistFlag == 'True':
        print "The following two files are saved as initialization data sets (latest month) for future mode:"
        print "ChStorage: monthly output, unit is m^3/month, dimension is", ChStorage.shape
        Settings.OutputNameStr = ChStorageNameStr
        SaveData(Settings, 'ChStorage_hist', ChStorage, flag)

        print "Soil column moisture: monthly output, unit is mm/month, dimension is", SO.shape
        Settings.OutputNameStr = ChStorageNameStr
        SaveData(Settings, 'Sav_hist', SO, flag)

    return Q, Avg_ChFlow


def SaveData(Settings, varstr, data, flag):
    filename = Settings.OutputFolder + varstr + "_" + Settings.OutputNameStr
    if flag == 0:
        SaveNetCDF(filename, data, Settings, varstr)
    else:
        SaveCSV(filename, data)


def SaveMAT(filename, data, varstr):
    filename = filename + ".mat"
    spio.savemat(filename, {varstr: data})


def SaveCSV(filename, data):
    filename = filename + ".csv"
    with open(filename, 'w') as outfile:
        np.savetxt(outfile, data, delimiter=',')


def SaveNetCDF(filename, data, Settings, varstr):
    filename = filename + ".nc"
    # open
    datagrp = Dataset(filename, 'w', format='NETCDF4')
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
    griddata[:, :] = data[:, :]

    # close
    datagrp.close()


def writecsvMap(filename, data, Settings):
    years = map(str, range(Settings.StartYear, Settings.EndYear + 1))
    headerline = "ID," + ",".join([year for year in years]) + ", Unit (km^3/year)"

    with open(filename + '.csv', 'w') as outfile:
        np.savetxt(outfile, data, delimiter=',', header=headerline, fmt='%s')