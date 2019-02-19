#!/usr/bin/env python
"""
Test that the hargreaves_gwam_mrtm configuration runs.
"""

from xanthos import Xanthos
import unittest
import numpy as np


class testHargreavesGwamMrtm(unittest.TestCase):
    """
    Test the hargreaves_gwam_mrtm configuration.

    Uses the config file test/configs/hargreaves_gwam_mrtm.ini which runs for
    one year with no outputs and no optional modules.
    """
    NCELL = 67420
    NMONTH = 12

    def setUp(self):
        pass

    def testHargreavesGwamMrtm(cls):
        """Test the hargreaves_gwam_mrtm configuration runs with reasonable input."""
        # Set up and run Xanthos
        ini = 'xanthos/test/configs/hargreaves_gwam_mrtm.ini'
        xth = Xanthos(ini)

        # Run with data that is within a reasonable range of actual inputs
        args = {
            "TemperatureFile": np.random.uniform(-5, 30, (cls.NCELL, cls.NMONTH)),
            "DailyTemperatureRangeFile": np.random.uniform(1, 10, (cls.NCELL, cls.NMONTH)),
            "PrecipitationFile": np.random.normal(50, 10, (cls.NCELL, cls.NMONTH))
        }

        res = xth.execute(args)

        cls.assertEqual(res.Q.shape, (cls.NCELL, cls.NMONTH))
        cls.assertFalse(np.any(np.isnan(res.Q)))
        cls.assertFalse(np.any(res.Q < 0))


if __name__ == '__main__':
    unittest.main()
