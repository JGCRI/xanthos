"""
Module to write output data files.

Created on Oct 11, 2016
Modified on Dec 10, 2018

@author: lixi729, Caleb Braun
@Project: Xanthos V2.2

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2018, Battelle Memorial Institute
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
FORMAT_NPY = 4

UNIT_MM_MTH = 0
UNIT_KM3_MTH = 1

NMONTHS = 12


class OutWriter:
    """
    Write out main Xanthos output variables.

    Output settings:
        OutputFormat  =  0 (default, netcdf file); 1 (csv file); 2 (mat file); 3 (parquet file)
        OutputUnit    =  0 (default, mm); 1 (km3)
        OutputInYear  =  0 (default, per month); 1 (per year, 12-month results combine into annual result)
    """

    def __init__(self, settings, grid_areas, all_outputs):
        """
        Initialize necessary settings for writing output variables.

        :param settings:        parsed settings from input configuration file
        :param grid_areas:      map of basin indices to grid cell area, in km2 (numpy array)
        :param all_outputs:     dictionary mapping all output names (strings) to their values (numpy arrays)
        """
        self.output_names = [oname for oname in settings.output_vars if oname in all_outputs.keys()]
        self.outputs = [pd.DataFrame(all_outputs[out_name]) for out_name in self.output_names]

        # array for converting basin values from mm to km3
        self.conversion_mm_km3 = grid_areas / 1e6

        # output options
        self.proj_name = settings.ProjectName
        self.out_folder = settings.OutputFolder
        self.out_format = settings.OutputFormat
        self.out_unit = settings.OutputUnit
        self.out_unit_str = settings.OutputUnitStr
        self.output_in_year = settings.OutputInYear
        self.ChStorageNameStr = settings.OutputNameStr

        year_range = range(settings.StartYear, settings.EndYear + 1)
        if self.output_in_year:
            self.time_steps = [str(y) for y in year_range]
        else:
            self.time_steps = ['{}{:02}'.format(y, mth) for y in year_range for mth in range(1, NMONTHS + 1)]

        # parameter checks
        if self.out_format not in [FORMAT_NETCDF, FORMAT_CSV, FORMAT_MAT, FORMAT_PARQUET, FORMAT_NPY]:
            logging.warning("Output format {} is invalid; writing output as .csv".format(self.out_format))
            self.out_format = FORMAT_CSV

    def get(self, varstr):
        """Get an output variable from its name."""
        return self.outputs[self.output_names.index(varstr)]

    def write(self):
        """
        Format and call appropriate writer for output variables.

        Writes all data to disk, using the settings defined in the config file.
        Assumes all variables are DataFrames containing values in units of
        mm/mth for sequential months over all years.

        The output for average channel flow is a special case. It is expected
        in units of m3/sec and no conversions will be made. Aggregating to
        yearly time steps is done by taking the mean of all months over a year.
        """
        if not self.output_names:
            # output_names will be an empty list if no outputs were requested
            logging.debug("No valid output variables specified")
            return

        # output_in_year is 0 for monthly output and 1 for yearly
        logging.debug("Outputting data {}".format(("monthly", "annually")[self.output_in_year]))
        logging.debug("Unit is {}".format(self.out_unit_str))

        for i, var in enumerate(self.output_names):
            if var == 'avgchflow':
                agg_function = 'mean'
                unit = 'm3persec'
            else:
                agg_function = 'sum'
                unit = self.out_unit_str

            # Take the yearly mean for aggregating channel flow, otherwise take the total
            if self.output_in_year:
                self.outputs[i] = self.agg_to_year(self.outputs[i], agg_function)

            if self.out_unit == UNIT_KM3_MTH and var != 'avgchflow':
                self.outputs[i] = self.outputs[i].multiply(self.conversion_mm_km3, axis=0)

            logging.debug("{} output dimension is {}".format(var, self.outputs[i].shape))

            filename = '{}_{}_{}'.format(var, unit, self.proj_name)
            filename = os.path.join(self.out_folder, filename)

            # Outputs are by grid cell and grid cell ids start at 1
            self.outputs[i].index += 1

            self.write_data(filename, var, self.outputs[i], col_names=self.time_steps)

    def write_aggregates(self, ref, df, basin, country, region):
        """
        Spatially aggregate runoff and write out results.

        :param ref:         parsed reference data
        :param df:          pandas DataFrame of values to aggregate and output
        :param basin:       bool - aggregate by basin?
        :param country:     bool - aggregate by country?
        :param region:      bool - aggregate by GCAM region?
        """
        filename = '{}_{}_{}'.format('{}', self.out_unit_str, self.proj_name)
        filepath = os.path.join(self.out_folder, filename)

        if basin:
            logging.info("Aggregating by Basin")
            varstr = 'Basin_runoff'
            bsn_agg = self.agg_spatial(df, ref.basin_ids, ref.basin_names, inc_name_idx=True)
            self.write_data(filepath.format(varstr), varstr, bsn_agg)

        if country:
            logging.info("Aggregating by Country")
            varstr = 'Country_runoff'
            ctry_agg = self.agg_spatial(df, ref.country_ids, ref.country_names)
            self.write_data(filepath.format(varstr), varstr, ctry_agg)

        if region:
            logging.info("Aggregating by GCAM Region")
            varstr = 'GCAMRegion_runoff'
            rgn_agg = self.agg_spatial(df, ref.region_ids, ref.region_names, inc_name_idx=True)
            self.write_data(filepath.format(varstr), varstr, rgn_agg)

        logging.info("Aggregated unit is {}".format(self.out_unit_str))

    def write_data(self, filename, var, data, col_names=None):
        """Save output data as a NetCDF or .csv, .mat, parquet, or .npy file."""
        if self.out_format == FORMAT_NETCDF:
            self.save_netcdf(filename, data, var)

        elif self.out_format == FORMAT_CSV:
            filename += ".csv"
            self.save_csv(filename, data, col_names)

        elif self.out_format == FORMAT_MAT:
            filename += ".mat"
            self.save_mat(filename, data, var)

        elif self.out_format == FORMAT_PARQUET:
            self.save_parquet(filename, data, col_names)

        elif self.out_format == FORMAT_NPY:
            np.save(filename, data)

    def save_mat(self, filename, data, varstr):
        """Write output data in the .mat format."""
        spio.savemat(filename, {varstr: data})

    def save_csv(self, filename, df, col_names=None, add_id=True):
        """Write pandas DataFrame as a .csv file."""
        if col_names is None:
            col_names = list(df.columns)

        df.columns = col_names

        if add_id:
            # add in index as region, basin, country or grid cell id number
            df.to_csv(filename, index_label='id')
        else:
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

    def save_parquet(self, filename, df, col_names=None):
        """Write pandas DataFrame to parquet file."""
        from fastparquet import write as fp_write

        filename += ".parquet"
        append = os.path.exists(filename)

        if col_names is not None:
            df.columns = col_names

        fp_write(filename, df, row_group_offsets=len(df), compression="GZIP", file_scheme='hive', has_nulls=False, append=append)

    def agg_to_year(self, df, func='sum'):
        """
        Aggregate a DataFrame (cells x months) to (cells x years).

        The groups are defined by an array the length of the number of columns,
        where identical values are grouped. Integer division by 12 on the
        column index gives sequential groups of length 12.

        :param df:      DataFrame with months in columns
        :param func:    Function to use for aggregation [default sum]
        """
        return df.groupby(np.arange(len(df.columns)) // NMONTHS, axis=1).agg(func)

    def agg_spatial(self, df, id_map, name_map, inc_name_idx=False):
        """Aggregate a DataFrame (cells x time) to (geographic area x time)."""
        # Convert to DataFrame for joining
        names_df = pd.DataFrame({'name': name_map})
        if inc_name_idx:
            names_df.index = names_df.index + 1

        # Aggregate all grid cells using the id map
        df['id'] = id_map
        agg_df = df.groupby('id', as_index=False).sum()

        # Map on the region/basin/country names, keeping all names even where there are no values
        agg_df = names_df.merge(agg_df, 'left', left_index=True, right_on='id')
        agg_df.set_index('id', inplace=True)

        return agg_df
