#!/usr/bin/env python
"""
Test that the default outputs do not change.
"""

from xanthos import Xanthos
import unittest
import pandas as pd
import glob


class TestEqualOutputs(unittest.TestCase):
    def setUp(self):
        self.ini = 'example/pm_abcd_mrtm.ini'
        self.xth = Xanthos(self.ini)

        # Allowable fraction difference between old and new data
        self.tolerance = 1e-8

        # Load original data
        self.old_files = self.read_outputs()

    def testRun(self):
        """Test that Xanthos runs."""
        res = self.xth.execute()
        self.assertEqual(res.Q.shape, (67420, 372))

    def testAllowableOutputs(self):
        """Test that new outputs equal old outputs."""
        new_files = self.read_outputs()

        for i in new_files:
            pd.testing.assert_frame_equal(new_files[i], self.old_files[i])

    def read_outputs(self):
        """Read all .csv files in output directory."""
        out_dir = 'xanthos/example/output/'
        out_file_names = glob.glob('{}*.csv'.format(out_dir))

        out_files = sorted([pd.read_csv(f) for f in out_file_names])
        return(out_files)


if __name__ == '__main__':
    unittest.main()
