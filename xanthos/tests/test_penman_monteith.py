import os
import pkg_resources
import unittest

import numpy as np
import pandas as pd

from xanthos.pet.penman_monteith import run_pmpet
from xanthos.data_reader.data_load import DataLoader


class TestPenmanMonteith(unittest.TestCase):
    """Test that the default outputs do not change."""

    INPUT_ELEV = pkg_resources.resource_filename('xanthos', 'tests/data/inputs/pet/penman_monteith/elev.npy')
    INPUT_ALBEDO = pkg_resources.resource_filename('xanthos', 'tests/data/inputs/pet/penman_monteith/gcam_albedo.csv')
    INPUT_ETPARS = pkg_resources.resource_filename('xanthos', 'tests/data/inputs/pet/penman_monteith/gcam_ET_para.csv')
    INPUT_LAI = pkg_resources.resource_filename('xanthos', 'tests/data/inputs/pet/penman_monteith/gcam_lai.csv')
    INPUT_LAIMAX = pkg_resources.resource_filename('xanthos', 'tests/data/inputs/pet/penman_monteith/gcam_laimax.csv')
    INPUT_LAIMIN = pkg_resources.resource_filename('xanthos', 'tests/data/inputs/pet/penman_monteith/gcam_laimin.csv')
    INPUT_LUCC = pkg_resources.resource_filename('xanthos', 'tests/data/inputs/pet/penman_monteith/lucc1901_2010_lump.npy')
    INPUT_RHS = pkg_resources.resource_filename('xanthos', 'tests/data/inputs/pet/penman_monteith/rhs_watch_monthly_percent_1971_1975.npy')
    INPUT_RLDS = pkg_resources.resource_filename('xanthos', 'tests/data/inputs/pet/penman_monteith/rlds_watch_monthly_wperm2_1971_1975.npy')
    INPUT_RSDS = pkg_resources.resource_filename('xanthos', 'tests/data/inputs/pet/penman_monteith/rsds_watch_monthly_wperm2_1971_1975.npy')
    INPUT_TAS = pkg_resources.resource_filename('xanthos', 'tests/data/inputs/pet/penman_monteith/tas_watch_monthly_degc_1971_1975.npy')
    INPUT_TASMIN = pkg_resources.resource_filename('xanthos', 'tests/data/inputs/pet/penman_monteith/tasmin_watch_monthly_degc_1971_1975.npy')
    INPUT_WIND = pkg_resources.resource_filename('xanthos', 'tests/data/inputs/pet/penman_monteith/wind_watch_monthly_mpers_1971_1975.npy')


    def test_outputs(self):
        """Ensure the outputs match what is expected."""

        data = LocalLoader()

        yr = 1971

        # pet_arr = run_pmpet(yr, data, )


class LocalLoader(DataLoader):
    """Data loading protocol from `data_load.py`

    TODO:  Replace the way the penman monteith code loads the data to be more accessible.
    TODO:  Make PM data load a method in data loader

    """
    
    INPUT_ELEV = pkg_resources.resource_filename('xanthos', 'tests/data/inputs/pet/penman_monteith/elev.npy')
    INPUT_ALBEDO = pkg_resources.resource_filename('xanthos', 'tests/data/inputs/pet/penman_monteith/gcam_albedo.csv')
    INPUT_ETPARS = pkg_resources.resource_filename('xanthos', 'tests/data/inputs/pet/penman_monteith/gcam_ET_para.csv')
    INPUT_LAI = pkg_resources.resource_filename('xanthos', 'tests/data/inputs/pet/penman_monteith/gcam_lai.csv')
    INPUT_LAIMAX = pkg_resources.resource_filename('xanthos', 'tests/data/inputs/pet/penman_monteith/gcam_laimax.csv')
    INPUT_LAIMIN = pkg_resources.resource_filename('xanthos', 'tests/data/inputs/pet/penman_monteith/gcam_laimin.csv')
    INPUT_LULC = pkg_resources.resource_filename('xanthos', 'tests/data/inputs/pet/penman_monteith/lucc1901_2010_lump.npy')
    INPUT_RHS = pkg_resources.resource_filename('xanthos', 'tests/data/inputs/pet/penman_monteith/rhs_watch_monthly_percent_1971_1975.npy')
    INPUT_RLDS = pkg_resources.resource_filename('xanthos', 'tests/data/inputs/pet/penman_monteith/rlds_watch_monthly_wperm2_1971_1975.npy')
    INPUT_RSDS = pkg_resources.resource_filename('xanthos', 'tests/data/inputs/pet/penman_monteith/rsds_watch_monthly_wperm2_1971_1975.npy')
    INPUT_TAS = pkg_resources.resource_filename('xanthos', 'tests/data/inputs/pet/penman_monteith/tas_watch_monthly_degc_1971_1975.npy')
    INPUT_TASMIN = pkg_resources.resource_filename('xanthos', 'tests/data/inputs/pet/penman_monteith/tasmin_watch_monthly_degc_1971_1975.npy')
    INPUT_WIND = pkg_resources.resource_filename('xanthos', 'tests/data/inputs/pet/penman_monteith/wind_watch_monthly_mpers_1971_1975.npy')

    def __init__(self):

        super().__init__('LocalLoader')
        
        # values from literature
        et_params = np.genfromtxt(LocalLoader.INPUT_ETPARS, delimiter=',')
    
        self.cL = et_params[:, 0]
        self.beta = et_params[:, 1]
        self.rslimit = et_params[:, 2]
    
        # correlation coefficient for calculating emissivity (range from 0.34 to 0.44), not the ABCD parameter a
        self.ae = et_params[:, 3]
    
        # correlation coefficient for calculating emissivity (range from -0.14 to -0.25), not the ABCD parameter b
        self.be = et_params[:, 4]
        self.Tminopen = et_params[:, 5]
        self.Tminclose = et_params[:, 6]
        self.VPDclose = et_params[:, 7]
        self.VPDopen = et_params[:, 8]
        self.RBLmin = et_params[:, 9]
        self.RBLmax = et_params[:, 10]
        self.rc = et_params[:, 11]
        self.emiss = et_params[:, 12]
    
        # 2-d, rows:11 land cover types, cols:12 months
        self.alpha = np.genfromtxt(LocalLoader.INPUT_ALBEDO, delimiter=',')
        self.lai = np.genfromtxt(LocalLoader.INPUT_LAI, delimiter=',')
        self.laimax = np.genfromtxt(LocalLoader.INPUT_LAIMAX, delimiter=',')
        self.laimin = np.genfromtxt(LocalLoader.INPUT_LAIMIN, delimiter=',')
    
        # convert missing values to 0
        self.tair_load = self.load_to_array(LocalLoader.INPUT_TAS, 'pm_tas', nan_to_num=True)
        self.TMIN_load = self.load_to_array(LocalLoader.INPUT_TASMIN, 'pm_tmin', nan_to_num=True)
        self.rhs_load = self.load_to_array(LocalLoader.INPUT_RHS, 'pm_rhs', nan_to_num=True)
        self.wind_load = self.load_to_array(LocalLoader.INPUT_WIND, 'pm_wind', nan_to_num=True)
        self.rsds_load = self.load_to_array(LocalLoader.INPUT_RSDS, 'pm_rsds', nan_to_num=True)
        self.rlds_load = self.load_to_array(LocalLoader.INPUT_RLDS, 'pm_rlds', nan_to_num=True)
    
        # set previous air temp value leaving the first value at 0
        self.tairprev_load = np.zeros_like(self.tair_load)
        self.tairprev_load[1:, :] = self.tair_load[:-1, :]
    
        # use land cover for each target year (nan to 0);  1-d:67420 cells, 2-d: land cover type, 3-d: years
        self.lct_load = np.nan_to_num(np.load(LocalLoader.INPUT_LULC))
    
        # static data
        self.elev = np.nan_to_num(np.load(LocalLoader.INPUT_ELEV))


if __name__ == '__main__':
    unittest.main()
