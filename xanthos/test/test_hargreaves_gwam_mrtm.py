#!/usr/bin/env python
"""
Test that the hargreaves_gwam_mrtm configuration runs.
"""

import os
import unittest
import pkg_resources
import numpy as np

from xanthos import Xanthos


class testHargreavesGwamMrtm(unittest.TestCase):
    """Test the hargreaves_gwam_mrtm configuration.

    Uses the config file test/configs/hargreaves_gwam_mrtm.ini which runs for
    one year with no outputs and no optional modules.

    """

    TEST_CONFIG_FILE = pkg_resources.resource_filename('xanthos', 'test/configs/hargreaves_gwam_mrtm.ini')

    NCELL = 67420
    NMONTH = 12

    def setUp(self):
        pass

    def testHargreavesGwamMrtm(self):
        """Test the hargreaves_gwam_mrtm configuration runs with reasonable input."""

        # Set up and run Xanthos
        xth = Xanthos(testHargreavesGwamMrtm.TEST_CONFIG_FILE)

        # Run with data that is within a reasonable range of actual inputs
        args = {
            "TemperatureFile": np.random.uniform(-5, 30, (testHargreavesGwamMrtm.NCELL, testHargreavesGwamMrtm.NMONTH)),
            "DailyTemperatureRangeFile": np.random.uniform(1, 10, (testHargreavesGwamMrtm.NCELL, testHargreavesGwamMrtm.NMONTH)),
            "PrecipitationFile": np.random.normal(50, 10, (testHargreavesGwamMrtm.NCELL, testHargreavesGwamMrtm.NMONTH))
        }

        res = xth.execute(args)

        self.assertEqual(res.Q.shape, (testHargreavesGwamMrtm.NCELL, testHargreavesGwamMrtm.NMONTH))
        self.assertFalse(np.any(np.isnan(res.Q)))
        self.assertFalse(np.any(res.Q < 0))


if __name__ == '__main__':
    unittest.main()
