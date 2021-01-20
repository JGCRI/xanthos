import os
import glob
import pkg_resources
import unittest

import pandas as pd

from xanthos import Xanthos


class TestCalibManaged(unittest.TestCase):

    """Test that the default outputs do not change."""

    DEFAULT_CONFIG_FILE = pkg_resources.resource_filename('xanthos', 'tests/data/inputs/pm_abcd_mrtm_managed.ini')
    EXAMPLE_DIR = pkg_resources.resource_filename('xanthos', 'test/data')
    DEFAULT_OUTPUTS_DIR = pkg_resources.resource_filename('xanthos', 'test/data/outputs/pm_abcd_mrtm_watch_1971_2001')

    def test_outputs(self):
        """Test that Xanthos produces correct outputs for the default configuration."""

        # run calibration
        xth = Xanthos(TestCalibManaged.DEFAULT_CONFIG_FILE)
        xth.execute()


if __name__ == '__main__':
    unittest.main()
