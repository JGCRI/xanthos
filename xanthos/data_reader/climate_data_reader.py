import calendar
import pkg_resources

import numpy as np
import xarray as xr
import pandas as pd


class ClimateToXanthos:
    """Convert common climate data 3D formats to the 2D format required by Xanthos. Unit
    conversions are conducted per metric as well.

    USAGE:

    To save formatted outputs as datasets:
    ```python
    from xanthos import ClimateToXanthos
    ctx = ClimateToXanthos(save_outputs=True)

    # for precipitation ('pr')
    pr_arr = ctx.format_pr(ncdf="<path to NetCDF file>",
                            data_variable_name='pr',
                            start_yr=<int four digit year>,
                            through_yr=<int four digit year>,
                            yr_interval=1,
                            out_file="<path and name to save NPY file to>")

    # for `tas`
    tas_arr = ctx.format_tas(ncdf"<path to NetCDF file>",
                              data_variable_name='tas',
                              out_file="<path and name to save NPY file to>")
    ```

    :param time_dimension_name:     Name of the time dimension in the NetCDF file.

    :param save_outputs:            If True, outputs are saved as pickled NPY files; otherwise,
                                    False.

    """

    XANTHOS_GRID_NUM = 67420
    SECONDS_IN_DAY = 86400
    KELVIN_TO_CELSIUS = -273.15
    DY_PER_MTH_FORMATS = {'standard': [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
                          'leap': [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]}

    # data path for reference data
    REF_LANDCELL_FILE = pkg_resources.resource_filename('xanthos', 'data/xanthos_0p5deg_landcell_reference.csv')
    COL_LON_INDEX = 'longitude_index'
    COL_LAT_INDEX = 'latitude_index'
    COL_BASIN_ID = 'basin_id'

    def __init__(self, time_dimension_name='time', save_outputs=False):

        self.time_dimension_name = time_dimension_name
        self.save_outputs = save_outputs

        # read in required reference data
        ref_df = pd.read_csv(ClimateToXanthos.REF_LANDCELL_FILE)

        # generate landcell index array
        self.lat_lon_array = ref_df[[ClimateToXanthos.COL_LAT_INDEX, ClimateToXanthos.COL_LON_INDEX]].values

        # get basin by landcell array
        self.basin_id_array = ref_df[[ClimateToXanthos.COL_BASIN_ID]][ClimateToXanthos.COL_BASIN_ID].values

    @classmethod
    def construct_days_per_month(cls, start_yr, through_yr, yr_interval=1):
        """Construct an array of the number of days per month for all years when
        considering leap years.

        :param start_yr:                Four digit integer for year dataset starts in.
                                        Dataset must start on January 01.

        :param through_yr:              Four digit integer for year dataset covers through.
                                        Dataset must end covering December 31.

        :param yr_interval:             Integer for the year interval. Default: 1

        """

        l = []
        for yr in range(start_yr, through_yr + yr_interval, yr_interval):

            if calendar.isleap(yr):
                days_per_month = cls.DY_PER_MTH_FORMATS['leap']

            else:
                days_per_month = cls.DY_PER_MTH_FORMATS['standard']

            l.extend(days_per_month)

        return np.array(l)

    @staticmethod
    def check_month_match(data_arr, start_yr, through_yr, yr_interval):
        """Check to see if number of months in NetCDF file match those of the
        input years.

        :param data_arr:                3D array containing [day, longitude, latitude] for
                                        each value.

        :param start_yr:                Four digit integer for year dataset starts in.
                                        Dataset must start on January 01.

        :param through_yr:              Four digit integer for year dataset covers through.
                                        Dataset must end covering December 31.

        :param yr_interval:             Integer for the year interval. Default: 1

        """

        msg = "Year range from {} through {} represents {} months.  NetCDF file has {} months."

        # number of months from year entry
        n_months = len(range(start_yr, through_yr + yr_interval, yr_interval)) * 12

        if n_months != data_arr.shape[0]:
            raise AssertionError(msg.format(start_yr, through_yr, n_months, data_arr.shape[0]))

    def replace_nan(self, data_arr):
        """Replace any nan elements in the array with the mean value for the basin.

        :param data_arr:                Full path with file name and extension to the reference
                                        basin id lookup file.  File has header.

        :return:                        Corrected data array.

        """

        # get the indicies of all nan values for the first month
        data_nan = np.where(np.isnan(data_arr[:, 0]))[0]

        # for each basin with a nan calculate a mean
        for ix in self.basin_id_array[data_nan]:
            # target basin means for each month
            rep_mean = np.nanmean(data_arr[self.basin_id_array == ix], axis=0)

            # set basin indicies up for all months
            basin_idx_all = np.tile(self.basin_id_array, (data_arr.shape[1], 1)).T

            # replace nan values in array with mean values for basin
            data_arr = np.where((np.isnan(data_arr)) & (basin_idx_all == ix), rep_mean, data_arr)

        return data_arr

    def extract_xanthos_grids(self, data_arr):
        """Extract each Xanthos land cell to an ordered array.

        :param data_arr:                3D array containing [day, latitude, longitude] for
                                        each value.

        :return:                        2D array ordered by Xanthos id for [grid_cell, month].

        """

        # set zeros array for output
        arr_out = np.zeros(shape=(ClimateToXanthos.XANTHOS_GRID_NUM, data_arr.shape[0]))

        for index, coord_idx in enumerate(self.lat_lon_array):
            lat_idx, lon_idx = coord_idx

            # get Xanthos land cell value by corresponding index
            arr_out[index, :] = data_arr[:, lat_idx, lon_idx]

        return arr_out

    def format_pr(self, ncdf, data_variable_name, start_yr, through_yr, yr_interval=1,
                  out_file=None):
        """Convert daily precipitation in units kg m-2 s-1 from a 3D shape where
        [day, latitude, longitude] to monthly mean precipitation shape
        where [grid_cell, month] in units mm per month.

        :param ncdf:                    Full path with file name and extension to the NetCDF
                                        file.

        :param data_variable_name:      Name of the NetCDF variable for the data

        :param start_yr:                Four digit integer for year dataset starts in.
                                        Dataset must start on January 01.

        :param through_yr:              Four digit integer for year dataset covers through.
                                        Dataset must end covering December 31.

        :param yr_interval:             Integer for the year interval. Default: 1

        :param out_file:                Full path with file name and extension for the output
                                        NPY file. Default: None.

        :return:                        Formatted Xanthos monthly mean precipitation in units
                                        mm per month.  [landcell, month]

        """

        # read in NetCDF file
        ds = xr.open_dataset(ncdf)

        # convert daily to monthly mean
        dsr = ds.resample(time='m').mean(self.time_dimension_name)

        data_arr = dsr.variables[data_variable_name][:]

        # check to make sure the number of months match
        self.check_month_match(data_arr, start_yr, through_yr, yr_interval)

        # extract xanthos grids in ordered array
        arr_out = self.extract_xanthos_grids(data_arr)

        # construct an array of the number of days per month for all years
        days_per_month_array = self.construct_days_per_month(start_yr, through_yr)

        # tile to apply to each month
        dpm_array = np.tile(days_per_month_array, (arr_out.shape[0], 1))

        # convert kg m-2 s-1 to mm per month
        arr_out *= ClimateToXanthos.SECONDS_IN_DAY  # to mm per day
        arr_out *= dpm_array  # to mm per month

        # replace nan values with means from corresponding basin
        arr_out = self.replace_nan(arr_out)

        if self.save_outputs:
            np.save(out_file, arr_out)

        return arr_out

    def format_tas(self, ncdf, data_variable_name, out_file=None):
        """Convert daily near surface air temperature (tas, tasmin, tasmax) in
        Kelvin from a 3D shape where [day, latitude, longitude] to monthly mean
        temperature for each Xanthos land grid cell where [grid_cell, month] in
        degrees Celsius.

        :param ncdf:                    Full path with file name and extension to the NetCDF
                                        file. Temperature units for input NetCDF are expected
                                        to be in Kelvin.

        :param data_variable_name:      Name of the NetCDF variable for the data

        :param out_file:                Full path with file name and extension for the output
                                        NPY file. Default: None.

        :return:                        Formatted Xanthos monthly mean temperature in units
                                        Celsius.  [landcell, month]

        """

        # read in NetCDF file
        self.ds = xr.open_dataset(ncdf)

        # convert daily to monthly mean
        dsr = self.ds.resample(time='m').mean(self.time_dimension_name)

        # replace nan values with means from corresponding basin
        data_arr = dsr.variables[data_variable_name][:]

        # extract xanthos grids in ordered array
        arr_out = self.extract_xanthos_grids(data_arr)

        # convert K to C
        arr_out += ClimateToXanthos.KELVIN_TO_CELSIUS

        # replace nan values with means from corresponding basin
        arr_out = self.replace_nan(arr_out)

        if self.save_outputs:
            np.save(out_file, arr_out)

        return arr_out
