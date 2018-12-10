"""
Perform diagnostics.

Perform diagnostics by comparing the estimates of average total annual
runoff (km^3/yr) of this study to other models.

Estimates of average total annual runoff (km^3/yr)
The comparison data file needs to be preprocessed.
Dimension: (67420, 1)
Unit: km3/year

Runoff
- VIC     The major comparison
- WBM     Ref comparison: WBM (Fekete et al., 2000) and WBMc (Fekete et al., 2000)
            are also used as additional comparisons (2 column csv files)
- UNH     Ref comparison: UNH-GRDC 1986-1995

Created on Jan 5, 2017
@author: lixi729
@Project: Xanthos V1.0

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute
"""

import numpy as np
import os
import pandas as pd
# import matplotlib.pyplot as plt


def Diagnostics(settings, Q, ref):
    """Aggregate and write results based on user settings."""
    if not settings.PerformDiagnostics:
        return

    # The name to use for plotting and comparison data output
    REF_DATA_NAME = 'VIC_1971-2000'

    area = ref.area

    # Diagnostics requested, so prepare the data
    nyear = int(settings.EndYear - settings.StartYear + 1)

    # convert the original unit mm/month to new unit km3/year
    q = np.sum(Q[:, :], axis=1) / nyear * area / 1e6

    VIC = ref.vic

    qq = np.mean(VIC, axis=1)

    UNH = ref.unh

    temp1 = ref.wbmd
    temp2 = ref.wbmc
    wbm = np.zeros((settings.ncell), dtype=float)
    wbmc = np.zeros((settings.ncell), dtype=float)
    for i in range(temp1.shape[0]):
        wbm[int(temp1[i, 0]) - 1] = temp1[i, 1]

    for i in range(temp2.shape[0]):
        wbmc[int(temp2[i, 0]) - 1] = temp2[i, 1]

    # Only basins/countries/regions for which all four models have
    # values are used to estimate the RMSE values

    # Basin Based
    if settings.DiagnosticScale == 0 or settings.DiagnosticScale == 1:
        Write_Diagnostics(settings, ref, 'Basin', q, qq, wbm, wbmc, UNH, REF_DATA_NAME)

    # Country Based
    if settings.DiagnosticScale == 0 or settings.DiagnosticScale == 2:
        Write_Diagnostics(settings, ref, 'Country', q, qq, wbm, wbmc, UNH, REF_DATA_NAME)

    # Region Based
    if settings.DiagnosticScale == 0 or settings.DiagnosticScale == 3:
        Write_Diagnostics(settings, ref, 'Region', q, qq, wbm, wbmc, UNH, REF_DATA_NAME)


def Aggregation_Diagnostics(settings, id_map, runoff):
    nregions = len(np.unique(id_map))
    Map_runoff = np.zeros((nregions,), dtype=float)

    for index in range(0, settings.ncell):
        if not np.isnan(runoff[index]) and Map[index] > 0:
            Map_runoff[Map[index] - 1] += runoff[index]

    return Map_runoff


def Write_Diagnostics(settings, ref, scale, q, qq, wbm, wbmc, UNH, plotname):
    if scale == 'Region':
        id_map = ref.region_ids
        name_map = ref.region_names
    elif scale == 'Basin':
        id_map = ref.basin_ids
        name_map = ref.basin_names
    elif scale == 'Country':
        id_map = ref.country_ids
        name_map = ref.country_names
    else:
        raise ValueError("Scale for diagnostics must be Region, Basin, or Country")

    runoff_df = pd.DataFrame({
        'xanthos': q,
        plotname: qq,
        'WBM': wbm,
        'WBMc': wbmc,
        'UNH_1986-1995': UNH
    })

    # Aggregate all grid cells using the id map
    runoff_df['id'] = id_map
    agg_df = runoff_df.groupby('id', as_index=False).sum()

    # Map on the region/basin/country names, keeping all names even where there are no values
    names_df = pd.DataFrame(name_map, columns=['name'])
    agg_df = names_df.merge(agg_df, 'left', left_index=True, right_on='id')
    agg_df.set_index('id', inplace=True)

    # Add global total as top row
    agg_df.loc[-1, 'name'] = 'Global'
    agg_df.loc[-1, 1:] = agg_df.sum(numeric_only=True)
    agg_df.index = agg_df.index + 1
    agg_df = agg_df.sort_index()

    file_name = "Diagnostics_Runoff_{}_Scale_km3peryr.csv".format(scale)
    output_name = os.path.join(settings.OutputFolder, file_name)

    # Write out as .csv, replacing nan values with zero
    agg_df.to_csv(output_name, na_rep=0, index=False)

    # Plot_Diagnostics(qs[1:, :], output_name, scale, plotname)


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
