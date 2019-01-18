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
        pass

    def testOutputs(self):
        """Test that Xanthos produces correct outputs."""
        # Remember original (correct) example outputs
        old_files = self.read_outputs()

        # Set up and run Xanthos
        ini = 'example/pm_abcd_mrtm.ini'
        xth = Xanthos(ini)
        res = xth.execute()

        # Check result dimensions
        self.assertEqual(res.Q.shape, (67420, 372))

        # Test that new outputs equal old outputs.
        new_files = self.read_outputs()
        for k in new_files.keys():
            pd.testing.assert_frame_equal(new_files[k], old_files[k])

    def read_outputs(self):
        """Read all .csv files in output directory."""
        out_dir = 'example/output/pm_abcd_mrtm_watch_1971_2001/'
        out_file_names = glob.glob('{}*.csv'.format(out_dir))

        out_files = {}
        for f in out_file_names:
            df = pd.read_csv(f)
            out_files[f] = df

        return(out_files)


if __name__ == '__main__':
    unittest.main()
