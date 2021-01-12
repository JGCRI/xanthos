import pkg_resources
import unittest

import numpy as np

from xanthos.data_reader.data_penman_monteith import DataPenmanMonteith
from xanthos.pet.penman_monteith import run_pmpet


class TestPenmanMonteith(unittest.TestCase):
    """Test that the default outputs do not change."""

    NCELLS = 67420
    NLCS = 8
    START_YR = 1971
    THROUGH_YR = 2001
    WATER_IDX = 0
    SNOW_IDX = 6
    LC_YEARS = [1900, 1910, 1920, 1930, 1940, 1950, 1960, 1970, 1980, 1990, 2000, 2005, 2010]

    INPUT_ELEV = pkg_resources.resource_filename('xanthos', 'tests/data/inputs/pet/penman_monteith/elev.npy')
    INPUT_ALBEDO = pkg_resources.resource_filename('xanthos', 'tests/data/inputs/pet/penman_monteith/gcam_albedo.csv')
    INPUT_ETPARS = pkg_resources.resource_filename('xanthos', 'tests/data/inputs/pet/penman_monteith/gcam_ET_para.csv')
    INPUT_LAI = pkg_resources.resource_filename('xanthos', 'tests/data/inputs/pet/penman_monteith/gcam_lai.csv')
    INPUT_LAIMAX = pkg_resources.resource_filename('xanthos', 'tests/data/inputs/pet/penman_monteith/gcam_laimax.csv')
    INPUT_LAIMIN = pkg_resources.resource_filename('xanthos', 'tests/data/inputs/pet/penman_monteith/gcam_laimin.csv')
    INPUT_LUCC = pkg_resources.resource_filename('xanthos', 'tests/data/inputs/pet/penman_monteith/lucc1901_2010_lump.npy')
    INPUT_RHS = pkg_resources.resource_filename('xanthos', f'tests/data/inputs/pet/penman_monteith/rhs_watch_monthly_percent_{START_YR}_{THROUGH_YR}.npy')
    INPUT_RLDS = pkg_resources.resource_filename('xanthos', f'tests/data/inputs/pet/penman_monteith/rlds_watch_monthly_wperm2_{START_YR}_{THROUGH_YR}.npy')
    INPUT_RSDS = pkg_resources.resource_filename('xanthos', f'tests/data/inputs/pet/penman_monteith/rsds_watch_monthly_wperm2_{START_YR}_{THROUGH_YR}.npy')
    INPUT_TAS = pkg_resources.resource_filename('xanthos', f'tests/data/inputs/pet/penman_monteith/tas_watch_monthly_degc_{START_YR}_{THROUGH_YR}.npy')
    INPUT_TASMIN = pkg_resources.resource_filename('xanthos', f'tests/data/inputs/pet/penman_monteith/tasmin_watch_monthly_degc_{START_YR}_{THROUGH_YR}.npy')
    INPUT_WIND = pkg_resources.resource_filename('xanthos', f'tests/data/inputs/pet/penman_monteith/wind_watch_monthly_mpers_{START_YR}_{THROUGH_YR}.npy')

    # TODO:  look up units for PET output
    COMP_PET = pkg_resources.resource_filename('xanthos', f'tests/data/comp_data/pet_penman_monteith_watch_mmpermth_{START_YR}_{THROUGH_YR}.npy')

    def test_outputs(self):
        """Ensure the outputs match what is expected."""

        # load module data
        data = DataPenmanMonteith(params_file=TestPenmanMonteith.INPUT_ETPARS,
                                  albedo_file=TestPenmanMonteith.INPUT_ALBEDO,
                                  lai_file=TestPenmanMonteith.INPUT_LAI,
                                  laimax_file=TestPenmanMonteith.INPUT_LAIMAX,
                                  laimin_file=TestPenmanMonteith.INPUT_LAIMIN,
                                  tas_file=TestPenmanMonteith.INPUT_TAS,
                                  tasmin_file=TestPenmanMonteith.INPUT_TASMIN,
                                  rhs_file=TestPenmanMonteith.INPUT_RHS,
                                  wind_file=TestPenmanMonteith.INPUT_WIND,
                                  rsds_file=TestPenmanMonteith.INPUT_RSDS,
                                  rlds_file=TestPenmanMonteith.INPUT_RLDS,
                                  lulc_file=TestPenmanMonteith.INPUT_LUCC,
                                  elev_file=TestPenmanMonteith.INPUT_ELEV,
                                  start_yr=TestPenmanMonteith.START_YR,
                                  through_yr=TestPenmanMonteith.THROUGH_YR)
        # run the module
        pet = run_pmpet(data,
                        ncells=TestPenmanMonteith.NCELLS,
                        nlcs=TestPenmanMonteith.NLCS,
                        start_yr=TestPenmanMonteith.START_YR,
                        end_yr=TestPenmanMonteith.THROUGH_YR,
                        water_idx=TestPenmanMonteith.WATER_IDX,
                        snow_idx=TestPenmanMonteith.SNOW_IDX,
                        land_cover_years=TestPenmanMonteith.LC_YEARS)

        np.save(TestPenmanMonteith.COMP_PET, pet)


if __name__ == '__main__':
    unittest.main()
