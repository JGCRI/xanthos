"""
Module to create diagnostic files and charts.

Created on Jan 5, 2017
Modified on Dec 10, 2018

@author: lixi729, Caleb Braun
@Project: Xanthos V2.2

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2018, Battelle Memorial Institute
"""

import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt


class Diagnostics:
    """
    Perform diagnostics by comparing the estimates of average total annual
    runoff (km^3/yr) of xanthos to other models.

    Estimates of average total annual runoff (km^3/yr)
    The comparison data file needs to be preprocessed.
    Dimension: (67420, 1)
    Unit: km3/year

    Runoff
    - VIC     The major comparison
    - WBM     Ref comparison: WBM (Fekete et al., 2000) and WBMc (Fekete et al., 2000)
              are also used as additional comparisons (2 column csv files)
    - UNH     Ref comparison: UNH-GRDC 1986-1995
    """
    def __init__(self, settings, xanthos_q, ref):
        """
        Aggregate and write results based on user settings.

        :param settings:        parsed settings from input configuration file
        :param xanthos_q:       runoff from a xanthos run
        :param ref:             parsed reference data
        """
        if not settings.PerformDiagnostics:
            return

        # The name to use for plotting and comparison data output
        self.REF_DATA_NAME = 'VIC_1971-2000'

        self.output_folder = settings.OutputFolder

        area = ref.area

        # Diagnostics requested, so prepare the data
        nyear = int(settings.EndYear - settings.StartYear + 1)

        # convert the original unit mm/month to new unit km3/year
        q = np.sum(xanthos_q[:, :], axis=1) / nyear * area / 1e6

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
            self.write_diagnostics('Basin', ref.basin_ids, ref.basin_names, q, qq, wbm, wbmc, UNH, 1)

        # Country Based
        if settings.DiagnosticScale == 0 or settings.DiagnosticScale == 2:
            self.write_diagnostics('Country', ref.country_ids, ref.country_names, q, qq, wbm, wbmc, UNH)

        # Region Based
        if settings.DiagnosticScale == 0 or settings.DiagnosticScale == 3:
            self.write_diagnostics('Region', ref.region_ids, ref.region_names, q, qq, wbm, wbmc, UNH, 1)

    def write_diagnostics(self, scale, id_map, name_map, q, qq, wbm, wbmc, UNH, name_map_offset=0):
        """
        Combine reference data sets to write out diagnostic files.

        :param scale:           level of aggregation for diagnostics, one of 'Basin', 'Country', or 'Region'
        :param id_map:          map of grid cells to basin/country/region ids
        :param name_map:        map of aggregation region ids to names
        :param q:               xanthos runoff
        :param qq:              reference runoff
        :param wbm:             WBM runoff
        :param wbmc:            WBMc runoff
        :param UNH:             UNH runoff (1986-1995)
        :param name_map_offset: how much to offset the names id values from their index, default 0
        """
        if scale not in ['Basin', 'Country', 'Region']:
            raise ValueError("Scale for diagnostics must be Region, Basin, or Country")

        # Columns must be explicitly specified to keep non-alphabetical order for Python versions < 3.6
        runoff_df = pd.DataFrame({
            'xanthos': q,
            self.REF_DATA_NAME: qq,
            'WBM': wbm,
            'WBMc': wbmc,
            'UNH_1986-1995': UNH
        }, columns=['xanthos', self.REF_DATA_NAME, 'WBM', 'WBMc', 'UNH_1986-1995'])

        # Aggregate all grid cells using the id map
        runoff_df['id'] = id_map
        agg_df = runoff_df.groupby('id', as_index=False).sum()

        # Map on the region/basin/country names, keeping all names even where there are no values
        names_df = pd.DataFrame({'name': name_map})
        names_df.index += name_map_offset
        agg_df = names_df.merge(agg_df, 'left', left_index=True, right_on='id')
        agg_df.set_index('id', inplace=True)

        # Add global total as top row
        agg_df.loc[-1, 'name'] = 'Global'
        agg_df.loc[-1, 1:] = agg_df.sum(numeric_only=True)
        agg_df.index = agg_df.index + 1
        agg_df = agg_df.sort_index()

        file_name = 'Diagnostics_Runoff_{}_Scale_km3peryr'.format(scale)
        output_name = os.path.join(self.output_folder, file_name)

        # Write out as .csv, replacing nan values with zero
        agg_df.to_csv(output_name + '.csv', na_rep=0, index=False)

        # self.plot_diagnostics(qs[1:, :], output_name, scale)

    def plot_diagnostics(self, data, outputname, titlestr):
        """Plot diagnostics."""
        fig = plt.figure()
        ax = plt.gca()
        ax.loglog([0.01, 100000], [0.01, 100000], 'grey')
        ax.scatter(data[:, 0], data[:, 1], c='black', alpha=0.5, edgecolors='none', label=self.REF_DATA_NAME)
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
