#!/usr/bin/env python

import glob
import unittest
import pkg_resources
import pandas as pd

from xanthos import Xanthos
from xanthos.install_supplement import InstallSupplement


class TestEqualOutputs(unittest.TestCase):
    """Test that the default outputs do not change."""

    DEFAULT_CONFIG_FILE = pkg_resources.resource_filename('xanthos', 'test/configs/pm_abcd_mrtm.ini')

    # install example supplement for testing in Travis-CI
    InstallSupplement(pkg_resources.resource_filename('xanthos', 'test/data'))

    EXAMPLE_DIR = pkg_resources.resource_filename('xanthos', 'test/data/example')
    DEFAULT_OUTPUTS_DIR = pkg_resources.resource_filename('xanthos', 'test/data/example/output/pm_abcd_mrtm_watch_1971_2001')

    def setUp(self):
        pass

    def testOutputs(self):
        """Test that Xanthos produces correct outputs for the default configuration."""

        # Remember original (correct) example outputs
        old_files = self.read_outputs()

        # Set up and run Xanthos
        xth = Xanthos(TestEqualOutputs.DEFAULT_CONFIG_FILE)
        res = xth.execute()

        # Check result dimensions
        self.assertEqual(res.Q.shape, (67420, 372))

        # Test that new outputs equal old outputs.
        new_files = self.read_outputs()
        for k in new_files.keys():
            pd.testing.assert_frame_equal(new_files[k], old_files[k])

    @classmethod
    def read_outputs(cls):
        """Read all .csv files in output directory."""

        out_file_names = glob.glob('{}*.csv'.format(cls.DEFAULT_OUTPUTS_DIR))

        out_files = {}
        for f in out_file_names:
            df = pd.read_csv(f)
            out_files[f] = df

        return out_files


if __name__ == '__main__':
    unittest.main()
