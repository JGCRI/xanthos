import pkg_resources

import numpy as np
import pandas as pd

from xanthos.data_reader.data_utils import DataUtils


class DataReference(DataUtils):

    # field names for the reference landcell file
    COL_GRID_ID = 'grid_id'
    COL_LONGITUDE = 'longitude'
    COL_LATITUDE = 'latitude'
    COL_LON_INDEX = 'longitude_index'
    COL_LAT_INDEX = 'latitude_index'
    COL_BASIN_ID = 'basin_id'
    COL_COUNTRY_ID = 'country_id'
    COL_REGION_ID = 'region_id'
    COL_AREA_HA = 'area_hectares'

    # data paths for reference data
    REF_LANDCELL_FILE = pkg_resources.resource_filename('xanthos', 'data/xanthos_0p5deg_landcell_reference.csv')
    REF_BASIN_NAME_FILE = pkg_resources.resource_filename('xanthos', 'data/BasinNames235.txt')
    REF_REGION_NAME_FILE = pkg_resources.resource_filename('xanthos', 'data/Rgn32Names.csv')
    REF_COUNTRY_NAME_FILE = pkg_resources.resource_filename('xanthos', 'data/country-names.csv')
    REF_COORDS_FILE = pkg_resources.resource_filename('xanthos', 'data/coordinates.csv')

    def __init__(self, config=None, nmonths=None):

        if config is None:
            super().__init__(nmonths=nmonths)

        else:
            super().__init__(nmonths=config.nmonths)

        # load Xanthos landcell reference data as data frame
        ref_landcell_df = pd.read_csv(DataReference.REF_LANDCELL_FILE)

        # Area value for each land grid cell: 67420 x 1, convert from ha to km2
        self.area = ref_landcell_df[DataReference.COL_AREA_HA].values * 0.01

        # Coordinates for flattened grid:  67420 x 5, the columns are ID#, lon, lat, ilon, ilat
        self.coords = self.load_data(DataReference.REF_COORDS_FILE)

        # Basin ID Map: 67420 x 1, 235 Basins
        self.basin_ids = ref_landcell_df[DataReference.COL_BASIN_ID].values

        # Corresponding to Basin ID Map, 235 Basin Names: 1D String Array
        self.basin_names = self.load_data(DataReference.REF_BASIN_NAME_FILE)

        # GCAM Region ID Map :  67420 x 1 (The nonag region table will be the 'primary' region assignment)
        self.region_ids = ref_landcell_df[DataReference.COL_REGION_ID].values

        # Corresponding to GCAM Region ID Map
        self.region_names = self.get_region_names(DataReference.REF_REGION_NAME_FILE)

        # Country ID Map : 67420 x 1 (249 countries: 1-249)
        self.country_ids = ref_landcell_df[DataReference.COL_COUNTRY_ID].values

        # Corresponding to Country ID Map, 0-248 index number and 249 Country Names: 2D String Array
        self.country_names = self.get_country_names(DataReference.REF_COUNTRY_NAME_FILE)

        self.latitude = np.copy(self.coords[:, 2])

        self.lat_radians = np.radians(self.latitude)

    @staticmethod
    def get_country_names(f):
        """Get an array of country names corresponding to GCAM countries."""
        with open(f, 'r') as f:
            country = f.read().splitlines()
            return np.array([i.split(',') for i in country])[:, 1]

    @staticmethod
    def get_region_names(f):
        """Get an array of region names corresponding to the GCAM region id map."""
        with open(f, 'r') as f:
            f.readline()
            region = f.read().split('\n')
            return np.array([i.split(',') for i in region])[:, 0]
