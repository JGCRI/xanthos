"""
Created on Jan 5, 2017
@author: lixi729
@Project: Xanthos V1.0


License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute


Perform diagnostics by comparing the estimates of average total annual runoff (km^3/yr) of this study to other models.

# Estimates of average total annual runoff (km^3/yr)
# The comparison data file needs to be preprocessed.
# Dimension: (67420, 1)
# Unit: km3/year
#
# Runoff
# - VIC     The major comparison
# - WBM     Ref comparison: WBM (Fekete et al., 2000) and WBMc (Fekete et al., 2000) are also used as additional comparisons (2 column csv files)
# - UNH     Ref comparison: UNH-GRDC 1986-1995
"""

import numpy as np
import os
import matplotlib.pyplot as plt

from xanthos.data_reader.data_load import load_const_griddata as loadfile


def Diagnostics(settings, Q, Avg_ChFlow, ref):


    area = ref.area

    if settings.PerformDiagnostics:
        # Prepare the data
        ny = int(settings.EndYear - settings.StartYear + 1)
        # convert the original unit mm/month to new unit km3/year
        q = np.sum(Q[:, :], axis=1) / ny * area / 1e6
        # ac  = np.sum(Avg_ChFlow[:,:], axis=1)/ny * area/1e6

        VIC = loadfile(settings.VICDataFile, 0, "q")  # 67420*30
        VICyears = range(1971, 2001)
        try:
            si = VICyears.index(settings.StartYear)
        except:
            si = 0
        try:
            ei = VICyears.index(settings.EndYear) + 1
        except:
            ei = 30

        qq = np.sum(VIC[:, si:ei], axis=1) / (ei - si)
        plotname = 'VIC_' + str(VICyears[si]) + '-' + str(VICyears[ei - 1])

        UNH = loadfile(settings.UNHDataFile, 0, "q")  # 67420*1

        temp1 = loadfile(settings.WBMDataFile, 0, "q")
        temp2 = loadfile(settings.WBMCDataFile, 0, "q")
        wbm = np.zeros((settings.ncell), dtype=float)
        wbmc = np.zeros((settings.ncell), dtype=float)
        for i in range(temp1.shape[0]):
            wbm[int(temp1[i, 0]) - 1] = temp1[i, 1]

        for i in range(temp2.shape[0]):
            wbmc[int(temp2[i, 0]) - 1] = temp2[i, 1]

        # Only basins/countries/regions for which all four models have values are used to estimate the RMSE values

        if settings.DiagnosticScale == 0 or settings.DiagnosticScale == 1:
            # Basin Based
            qb = np.zeros((max(ref.basin_ids), 5), dtype=float)
            qb[:, 0] = Aggregation_Diagnostics(settings, ref.basin_ids, q)
            qb[:, 1] = Aggregation_Diagnostics(settings, ref.basin_ids, qq)
            qb[:, 2] = Aggregation_Diagnostics(settings, ref.basin_ids, wbm)
            qb[:, 3] = Aggregation_Diagnostics(settings, ref.basin_ids, wbmc)
            qb[:, 4] = Aggregation_Diagnostics(settings, ref.basin_ids, UNH)
            qb = np.insert(qb, 0, np.sum(qb, axis=0), axis=0)  # add global
            BasinNames = np.insert(ref.basin_names, 0, 'Global')

            writecsvDiagnostics(os.path.join(settings.OutputFolder, "Diagnostics_Runoff_Basin_Scale"), qb, plotname, BasinNames)

            outputname = os.path.join(settings.OutputFolder, "Diagnostics_Runoff_Basin_Scale")

            Plot_Diagnostics(qb[1:, :], outputname, 'Basin', plotname)

            for i in range(qb.shape[0]):
                if not (qb[i, 0] > 0 and qb[i, 1] > 0 and qb[i, 2] > 0 and qb[i, 3] > 0 and qb[i, 4] > 0):
                    qb[i, :] = 0

        if settings.DiagnosticScale == 0 or settings.DiagnosticScale == 2:
            # Country Based
            qc = np.zeros((max(ref.country_ids), 5), dtype=float)
            qc[:, 0] = Aggregation_Diagnostics(settings, ref.country_ids, q)
            qc[:, 1] = Aggregation_Diagnostics(settings, ref.country_ids, qq)
            qc[:, 2] = Aggregation_Diagnostics(settings, ref.country_ids, wbm)
            qc[:, 3] = Aggregation_Diagnostics(settings, ref.country_ids, wbmc)
            qc[:, 4] = Aggregation_Diagnostics(settings, ref.country_ids, UNH)
            qc = np.insert(qc, 0, np.sum(qc, axis=0), axis=0)  # add global
            CountryNames = np.insert(ref.country_names, 0, 'Global')

            writecsvDiagnostics(os.path.join(settings.OutputFolder, "Diagnostics_Runoff_Country_Scale"), qc, plotname, CountryNames)
            outputname = os.path.join(settings.OutputFolder, "Diagnostics_Runoff_Country_Scale")

            Plot_Diagnostics(qc[1:, :], outputname, 'Country', plotname)

            for i in range(qc.shape[0]):
                if not (qc[i, 0] > 0 and qc[i, 1] > 0 and qc[i, 2] > 0 and qc[i, 3] > 0 and qc[i, 4] > 0):
                    qc[i, :] = 0

                    # q1 = qc[np.nonzero(qc[:,0])[0][1:],0]
                    # q2 = qc[np.nonzero(qc[:,1])[0][1:],1]
                    # q3 = qc[np.nonzero(qc[:,2])[0][1:],2]
                    # q4 = qc[np.nonzero(qc[:,3])[0][1:],3]
                    # q5 = qc[np.nonzero(qc[:,4])[0][1:],4]

                    # print "RMSE at the country scale:          ", np.sqrt(((q1 - q2) ** 2).mean()), np.sqrt(((q1 - q3) ** 2).mean()), np.sqrt(((q1 - q4) ** 2).mean()), np.sqrt(((q1 - q5) ** 2).mean())

        if settings.DiagnosticScale == 0 or settings.DiagnosticScale == 3:
            # Region Based
            qr = np.zeros((max(ref.region_ids), 5), dtype=float)
            qr[:, 0] = Aggregation_Diagnostics(settings, ref.region_ids, q)
            qr[:, 1] = Aggregation_Diagnostics(settings, ref.region_ids, qq)
            qr[:, 2] = Aggregation_Diagnostics(settings, ref.region_ids, wbm)
            qr[:, 3] = Aggregation_Diagnostics(settings, ref.region_ids, wbmc)
            qr[:, 4] = Aggregation_Diagnostics(settings, ref.region_ids, UNH)
            qr = np.insert(qr, 0, np.sum(qr, axis=0), axis=0)  # add global
            RegionNames = np.insert(ref.region_names, 0, 'Global')

            writecsvDiagnostics(os.path.join(settings.OutputFolder, "Diagnostics_Runoff_Region_Scale"), qr, plotname, RegionNames)
            outputname = os.path.join(settings.OutputFolder, "Diagnostics_Runoff_Region_Scale")

            Plot_Diagnostics(qr[1:, :], outputname, 'Region', plotname)

            for i in range(qr.shape[0]):
                if not (qr[i, 0] > 0 and qr[i, 1] > 0 and qr[i, 2] > 0 and qr[i, 3] > 0 and qr[i, 4] > 0):
                    qr[i, :] = 0

    else:
        return


def Aggregation_Diagnostics(settings, Map, runoff):
    NB = max(Map)
    Map_runoff = np.zeros((NB,), dtype=float)

    for index in range(0, settings.ncell):
        if not np.isnan(runoff[index]) and Map[index] > 0:
            Map_runoff[Map[index] - 1] += runoff[index]

    return Map_runoff


def writecsvDiagnostics(filename, data, ComparisonDataName, Names):
    headerline = "Name,This Study," + ComparisonDataName + ",WBM,WBMc,UNH_1986-1995,Unit(km^3/year)"
    newdata = np.insert(data.astype(str), 0, Names, axis=1)

    with open(filename + '.csv', 'w') as outfile:
        np.savetxt(outfile, newdata, delimiter=',', header=headerline, fmt='%s')


def Plot_Diagnostics(data, outputname, titlestr, ComparisonDataName):
    fig = plt.figure()
    ax = plt.gca()
    ax.loglog([0.01, 100000], [0.01, 100000], 'grey')
    ax.scatter(data[:, 0], data[:, 1], c='black', alpha=0.5, edgecolors='none', label=ComparisonDataName)
    ax.scatter(data[:, 0], data[:, 2], c='Red', alpha=0.5, edgecolors='none', label='WBM')
    ax.scatter(data[:, 0], data[:, 3], c='Blue', alpha=0.5, edgecolors='none', label='WBMc')
    ax.scatter(data[:, 0], data[:, 4], c='green', alpha=0.5, edgecolors='none', label='UNH/GRDC_1986-1995')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.axis([0.01, 1e5, 0.01, 1e5])
    ax.legend(loc='lower right', bbox_to_anchor=(1, 0), fontsize=10)
    plt.title('Hydro Model Diagnostics at ' + titlestr + ' Scale', fontsize=12, fontweight='bold')
    plt.xlabel(r'This Study Estimated Averaged Annual Runoff ($km^3$/yr)', fontsize=12)
    plt.ylabel(r'Averaged Annual Runoff ($km^3$/yr)', fontsize=12)
    fig.savefig(outputname + '.png', dpi=300)
    plt.close(fig)