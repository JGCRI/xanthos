import numpy as np
from xanthos.data_reader.data_utils import DataUtils


class DataGwam(DataUtils):

    def __init__(self, config, area, region_ids, country_ids, basin_ids):

        self.config = config
        self.area = area
        self.region_ids = region_ids
        self.country_ids = country_ids
        self.basin_ids = basin_ids

        # monthly precipitation mm/mth
        self.precip = self.load_to_array(self.config.PrecipitationFile, var_name=self.config.PrecipVarName)

        # Max Soil Moisture Map (mm/month): 67420 x 1
        self.max_soil_moist = self.load_data(self.config.max_soil_moisture, 1)

        # Water Bodies: assign MSM = 999, 306 x 2, Col 1 is the cell number in 67420
        self.lakes_msm = self.load_data(self.config.lakes_msm).astype(int)
        self.lakes_msm[:, 0] -= 1

        # Additional water bodies: assign MSM = 999, 421 x 2,  Col 1 is the cell number in 67420
        self.addit_water_msm = self.load_data(self.config.addit_water_msm).astype(int)
        self.addit_water_msm[:, 0] -= 1

        # create a matrix (MSMC: Maximum Soil Moisture Capacity) with all data
        # 1: area; 2: region; 3: Max Soil Moisture (mm/month)
        msmc = self.load_soil_moisture()

        # harmonized grid area
        self.grid_area = np.copy(msmc[:, 0])

        # maximum soil moisture
        self.soil_moisture = np.copy(msmc[:, 2])

        # harmonized grid area
        self.grid_area = np.copy(msmc[:, 0])

        # load soil moisture file if running future
        if self.config.HistFlag.lower() == "true":
            self.sm_prev = 0.5 * self.soil_moisture
        else:
            self.sm_prev = self.load_soil_data()

    def load_soil_moisture(self, missing=-9999):
        """
        Load soil moisture data.

        Assign max soil moisture (mm/month) [2] to Sm.  For historic data use 0.5 * sm to an initial value to pass to
        runoff model. If future mode, read values from historical file.
        """
        data = np.zeros((self.config.ncell, 5), order='F')

        data[:, 0] = self.area
        data[:, 1] = self.region_ids
        data[:, 2] = self.max_soil_moist

        # add max value (999) where water is
        data[self.lakes_msm[:, 0], 2] = self.lakes_msm[:, 1]
        data[self.addit_water_msm[:, 0], 2] = self.addit_water_msm[:, 1]

        country = self.country_ids[:]
        basin = self.basin_ids[:]

        # Ignore all the cells in which we are missing an ID value for soil moisture, country, or basin.
        # Thus, country and basin coverage must be consistent.
        # Basins coverage is smaller, and GCAM region ignores Greenland.
        invalid = np.where((data[:, 2] == 0) | (country == 0) | (basin == 0))[0]

        # should this be 0:2
        data[invalid, 1:2] = 0

        # should these be returned?
        country[invalid] = missing
        basin[invalid] = missing

        return data

    def load_soil_data(self):
        """Load soil moisture file into array if in future mode, else stage zeros array."""
        try:
            # Initialize channel storage/soil moisture.
            if self.config.HistFlag.lower() == "true":
                return np.zeros((self.config.ncell,), dtype=float)

            # For future runs, initialize with the last value of the historical channel storage/soil moisture
            else:
                return self.load_data(self.config.SavFile, 0, self.config.SavVarName)[:, -1]

        # if not in use
        except AttributeError:
            return np.zeros((self.config.ncell,), dtype=float)
