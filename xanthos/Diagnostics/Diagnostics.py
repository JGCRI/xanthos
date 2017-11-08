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

from xanthos.data_reader.DataLoad import load_const_griddata as loadfile


def Diagnostics(settings, Q, Avg_ChFlow, GridConstants):


    area = GridConstants['Area']

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
        # print "Runoff Deviation from Estimates of this study: " + plotname + "          WBM          WBMc          UNH_1986-1995"

        if settings.DiagnosticScale == 0 or settings.DiagnosticScale == 1:
            # Basin Based
            qb = np.zeros((max(GridConstants['BasinIDs']), 5), dtype=float)
            qb[:, 0] = Aggregation_Diagnostics(settings, GridConstants['BasinIDs'], q)
            qb[:, 1] = Aggregation_Diagnostics(settings, GridConstants['BasinIDs'], qq)
            qb[:, 2] = Aggregation_Diagnostics(settings, GridConstants['BasinIDs'], wbm)
            qb[:, 3] = Aggregation_Diagnostics(settings, GridConstants['BasinIDs'], wbmc)
            qb[:, 4] = Aggregation_Diagnostics(settings, GridConstants['BasinIDs'], UNH)
            qb = np.insert(qb, 0, np.sum(qb, axis=0), axis=0)  # add global
            BasinNames = np.insert(GridConstants['BasinNames'], 0, 'Global')

            writecsvDiagnostics(os.path.join(settings.OutputFolder, "Diagnostics_Runoff_Basin_Scale"), qb, plotname, BasinNames)

            outputname = os.path.join(settings.OutputFolder, "Diagnostics_Runoff_Basin_Scale")

            Plot_Diagnostics(qb[1:, :], outputname, 'Basin', plotname)

            for i in range(qb.shape[0]):
                if not (qb[i, 0] > 0 and qb[i, 1] > 0 and qb[i, 2] > 0 and qb[i, 3] > 0 and qb[i, 4] > 0):
                    qb[i, :] = 0

                    # q1 = qb[np.nonzero(qb[:,0])[0][1:],0]
                    # q2 = qb[np.nonzero(qb[:,1])[0][1:],1]
                    # q3 = qb[np.nonzero(qb[:,2])[0][1:],2]
                    # q4 = qb[np.nonzero(qb[:,3])[0][1:],3]
                    # q5 = qb[np.nonzero(qb[:,4])[0][1:],4]

                    # print "RMSE at the basin scale:            ", np.sqrt(((q1 - q2) ** 2).mean()), np.sqrt(((q1 - q3) ** 2).mean()), np.sqrt(((q1 - q4) ** 2).mean()), np.sqrt(((q1 - q5) ** 2).mean())

        if settings.DiagnosticScale == 0 or settings.DiagnosticScale == 2:
            # Country Based
            qc = np.zeros((max(GridConstants['CountryIDs']), 5), dtype=float)
            qc[:, 0] = Aggregation_Diagnostics(settings, GridConstants['CountryIDs'], q)
            qc[:, 1] = Aggregation_Diagnostics(settings, GridConstants['CountryIDs'], qq)
            qc[:, 2] = Aggregation_Diagnostics(settings, GridConstants['CountryIDs'], wbm)
            qc[:, 3] = Aggregation_Diagnostics(settings, GridConstants['CountryIDs'], wbmc)
            qc[:, 4] = Aggregation_Diagnostics(settings, GridConstants['CountryIDs'], UNH)
            qc = np.insert(qc, 0, np.sum(qc, axis=0), axis=0)  # add global
            CountryNames = np.insert(GridConstants['CountryNames'], 0, 'Global')

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
            qr = np.zeros((max(GridConstants['GCAMRegionIDs']), 5), dtype=float)
            qr[:, 0] = Aggregation_Diagnostics(settings, GridConstants['GCAMRegionIDs'], q)
            qr[:, 1] = Aggregation_Diagnostics(settings, GridConstants['GCAMRegionIDs'], qq)
            qr[:, 2] = Aggregation_Diagnostics(settings, GridConstants['GCAMRegionIDs'], wbm)
            qr[:, 3] = Aggregation_Diagnostics(settings, GridConstants['GCAMRegionIDs'], wbmc)
            qr[:, 4] = Aggregation_Diagnostics(settings, GridConstants['GCAMRegionIDs'], UNH)
            qr = np.insert(qr, 0, np.sum(qr, axis=0), axis=0)  # add global
            RegionNames = np.insert(GridConstants['GCAMRegionNames'], 0, 'Global')

            writecsvDiagnostics(os.path.join(settings.OutputFolder, "Diagnostics_Runoff_Region_Scale"), qr, plotname, RegionNames)
            outputname = os.path.join(settings.OutputFolder, "Diagnostics_Runoff_Region_Scale")

            Plot_Diagnostics(qr[1:, :], outputname, 'Region', plotname)

            for i in range(qr.shape[0]):
                if not (qr[i, 0] > 0 and qr[i, 1] > 0 and qr[i, 2] > 0 and qr[i, 3] > 0 and qr[i, 4] > 0):
                    qr[i, :] = 0

                    # q1 = qr[np.nonzero(qr[:,0])[0][1:],0]
                    # q2 = qr[np.nonzero(qr[:,1])[0][1:],1]
                    # q3 = qr[np.nonzero(qr[:,2])[0][1:],2]
                    # q4 = qr[np.nonzero(qr[:,3])[0][1:],3]
                    # q5 = qr[np.nonzero(qr[:,4])[0][1:],4]

                    # print "RMSE at the region scale:           ", np.sqrt(((q1 - q2) ** 2).mean()), np.sqrt(((q1 - q3) ** 2).mean()), np.sqrt(((q1 - q4) ** 2).mean()), np.sqrt(((q1 - q5) ** 2).mean())

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