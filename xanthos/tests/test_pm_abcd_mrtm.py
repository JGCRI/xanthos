import os
import glob
import pkg_resources
import unittest

import pandas as pd

from xanthos import Xanthos


class TestPmAbcdMrtm(unittest.TestCase):

    """Test that the default outputs do not change."""

    DEFAULT_CONFIG_FILE = pkg_resources.resource_filename('xanthos', 'tests/data/inputs/pm_abcd_mrtm.ini')
    EXAMPLE_DIR = pkg_resources.resource_filename('xanthos', 'test/data')
    DEFAULT_OUTPUTS_DIR = pkg_resources.resource_filename('xanthos', 'test/data/outputs/pm_abcd_mrtm_watch_1971_2001')

    def test_outputs(self):
        """Test that Xanthos produces correct outputs for the default configuration."""

        # create output directory if not exist
        if not os.path.exists(TestPmAbcdMrtm.DEFAULT_OUTPUTS_DIR):
            os.makedirs(TestPmAbcdMrtm.DEFAULT_OUTPUTS_DIR)

        # Remember original (correct) example outputs
        old_csv_files = self.read_csv_outputs(data_dir='/Users/d3y010/repos/github/xanthos/example/output/pm_abcd_mrtm_watch_1971_2001')

        # Set up and run Xanthos
        xth = Xanthos(TestPmAbcdMrtm.DEFAULT_CONFIG_FILE)
        res = xth.execute()

        # Check result dimensions
        self.assertEqual(res.Q.shape, (67420, 372))

        # Test that new outputs equal old outputs.
        new_csv_files = self.read_csv_outputs()
        for k in new_csv_files.keys():
            pd.testing.assert_frame_equal(new_csv_files[k], old_csv_files[k])

    @classmethod
    def read_csv_outputs(cls, data_dir=None):
        """Read all .csv files in output directory."""

        if data_dir is None:
            out_file_names = glob.glob(f'{cls.DEFAULT_OUTPUTS_DIR}*.csv')
        else:
            out_file_names = glob.glob(f'{data_dir}*.csv')

        out_files = {}
        for f in out_file_names:
            df = pd.read_csv(f)
            out_files[f] = df

        return out_files


if __name__ == '__main__':
    unittest.main()
