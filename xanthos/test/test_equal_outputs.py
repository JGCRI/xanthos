#!/usr/bin/env python
"""
Test that the default outputs do not change.
"""

import os
import glob
import unittest
import pandas as pd

from xanthos import Xanthos


class TestEqualOutputs(unittest.TestCase):

    DEFAULT_CONFIG_FILE = os.path.join(os.path.dirname(__file__), 'configs', 'pm_abcd_mrtm.ini')

    # installed in setup.py as fetch from minted Zenodo supplementary data
    EXAMPLE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', 'example'))
    DEFAULT_OUTPUTS_DIR = os.path.join(EXAMPLE_DIR, 'output', 'pm_abcd_mrtm_watch_1971_2001')

    # def setUp(self):
    #     pass

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
