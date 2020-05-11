#!/usr/bin/env python

import unittest

from xanthos import ClimateToXanthos


class TestEqualOutputs(unittest.TestCase):
    """Test that the default outputs do not change."""

    XANTHOS_GRID_NUM = 67420
    SECONDS_IN_DAY = 86400
    KELVIN_TO_CELSIUS = -273.15
    DY_PER_MTH_FORMATS = {'standard': [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
                          'leap': [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]}

    def setUp(self):
        pass

    def testAttributes(self):
        """Test ClimateToXanthos attributes."""

        ctx = ClimateToXanthos()

        self.assertEqual(TestEqualOutputs.XANTHOS_GRID_NUM, ctx.XANTHOS_GRID_NUM)

        self.assertEqual(TestEqualOutputs.SECONDS_IN_DAY, ctx.SECONDS_IN_DAY)

        self.assertEqual(TestEqualOutputs.KELVIN_TO_CELSIUS, ctx.KELVIN_TO_CELSIUS)

        self.assertDictEqual(TestEqualOutputs.DY_PER_MTH_FORMATS, ctx.DY_PER_MTH_FORMATS)

        # check lat, lon index array shape
        self.assertEqual(ctx.lat_lon_array.shape, (67420, 2))

        # check basin id array shape
        self.assertEqual(ctx.basin_id_array.shape, (67420, ))


if __name__ == '__main__':
    unittest.main()
