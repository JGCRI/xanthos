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
import pandas as pd


# import matplotlib.pyplot as plt


def Diagnostics(settings, Q, ref):
    area = ref.area

    if not settings.PerformDiagnostics:
        return

    # Diagnostics requested, so prepare the data
    ny = int(settings.EndYear - settings.StartYear + 1)

    # convert the original unit mm/month to new unit km3/year
    q = np.sum(Q[:, :], axis=1) / ny * area / 1e6

    VIC = ref.vic

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

    UNH = ref.unh

    temp1 = ref.wbmd
    temp2 = ref.wbmc
    wbm = np.zeros((settings.ncell), dtype=float)
    wbmc = np.zeros((settings.ncell), dtype=float)
    for i in range(temp1.shape[0]):
        wbm[int(temp1[i, 0]) - 1] = temp1[i, 1]

    for i in range(temp2.shape[0]):
        wbmc[int(temp2[i, 0]) - 1] = temp2[i, 1]

    # Only basins/countries/regions for which all four models have values are used to estimate the RMSE values

    # Basin Based
    if settings.DiagnosticScale == 0 or settings.DiagnosticScale == 1:
        Write_Diagnostics(settings, ref, 'Basin', q, qq, wbm, wbmc, UNH, plotname)

    # Country Based
    if settings.DiagnosticScale == 0 or settings.DiagnosticScale == 2:
        Write_Diagnostics(settings, ref, 'Country', q, qq, wbm, wbmc, UNH, plotname)

    # Region Based
    if settings.DiagnosticScale == 0 or settings.DiagnosticScale == 3:
        Write_Diagnostics(settings, ref, 'Region', q, qq, wbm, wbmc, UNH, plotname)


def Aggregation_Diagnostics(settings, Map, runoff):
    NB = max(Map)
    Map_runoff = np.zeros((NB,), dtype=float)

    for index in range(0, settings.ncell):
        if not np.isnan(runoff[index]) and Map[index] > 0:
            Map_runoff[Map[index] - 1] += runoff[index]

    return Map_runoff


def Write_Diagnostics(settings, ref, scale, q, qq, wbm, wbmc, UNH, plotname):
    if scale == 'Region':
        ids = ref.region_ids
        names = ref.region_names
    elif scale == 'Basin':
        ids = ref.basin_ids
        names = ref.basin_names
    elif scale == 'Country':
        ids = ref.country_ids
        names = ref.country_names
    else:
        raise ValueError("Scale for diagnostics must be Region, Basin, or Country")

    qs = np.zeros((max(ids), 5), dtype=float)
    qs[:, 0] = Aggregation_Diagnostics(settings, ids, q)
    qs[:, 1] = Aggregation_Diagnostics(settings, ids, qq)
    qs[:, 2] = Aggregation_Diagnostics(settings, ids, wbm)
    qs[:, 3] = Aggregation_Diagnostics(settings, ids, wbmc)
    qs[:, 4] = Aggregation_Diagnostics(settings, ids, UNH)
    qs = np.insert(qs, 0, np.sum(qs, axis=0), axis=0)  # add global
    ScaleNames = np.insert(names, 0, 'Global')

    fname = "Diagnostics_Runoff_{}_Scale_km3peryr".format(scale)
    writecsvDiagnostics(os.path.join(settings.OutputFolder, fname), qs, plotname, ScaleNames)
    outputname = os.path.join(settings.OutputFolder, fname)

    # Plot_Diagnostics(qs[1:, :], outputname, scale, plotname)

    for i in range(qs.shape[0]):
        if not (qs[i, 0] > 0 and qs[i, 1] > 0 and qs[i, 2] > 0 and qs[i, 3] > 0 and qs[i, 4] > 0):
            qs[i, :] = 0


def writecsvDiagnostics(filename, data, ComparisonDataName, Names):
    hdr = "name,xanthos," + ComparisonDataName + ",WBM,WBMc,UNH_1986-1995"
    newdata = np.insert(data.astype(str), 0, Names, axis=1)

    df = pd.DataFrame(newdata)
    df.columns = hdr.split(',')

    if filename[-4:] != '.csv':
        filename = filename + '.csv'

    df.to_csv(filename, index=False)


#    with open(filename + '.csv', 'w') as outfile:
#        np.savetxt(outfile, newdata, delimiter=',', header=headerline, fmt='%s')


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
