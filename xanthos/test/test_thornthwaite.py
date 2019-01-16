#!/usr/bin/env python
"""
Test that the Thornthwaite PET module outputs correctly.
"""

from xanthos.pet import thornthwaite
import unittest
import numpy as np


class TestThornthwaite(unittest.TestCase):
    def testDaylightHours(self):
        """Test the daylight hour calculation method."""
        mth_days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

        equator = thornthwaite.calc_daylight_hours(mth_days, np.array([0]))
        north_pole = thornthwaite.calc_daylight_hours(mth_days, np.array([np.pi / 2]))
        south_pole = thornthwaite.calc_daylight_hours(mth_days, np.array([-np.pi / 2]))

        self.assertTrue(np.all(equator == 12))  # 12 hour days on equator
        self.assertTrue(np.any(north_pole[0] == 0.0))  # months with no daylight
        self.assertTrue(np.any(north_pole[0] == 24.0))  # months with all daylight
        self.assertTrue(np.all(24 - north_pole == south_pole))  # poles should be opposite

    def testThornthwaiteResults(self):
        """Make sure Thornthwaite results are as expected."""
        lat_radians = np.array([0.698132])   # 40 degrees N
        tas1 = np.array([[2,  5,  6,  8, 10, 12, 15, 12, 10,  8,  6,  5]])
        tas2 = np.array([[-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]])

        pet1_correct = np.array([[9.7, 22.9, 33.7, 47.6, 66.0, 78.8, 98.5, 74.4, 54.8, 40.9, 27.1, 22.3]])

        pet1 = thornthwaite.execute(tas1, lat_radians, 1999, 1999)  # test for one year
        pet2 = thornthwaite.execute(tas2, lat_radians, 1999, 1999)  # test for one year

        self.assertTrue(np.all(np.round(pet1, 1) == pet1_correct))
        self.assertTrue(np.all(pet2 == 0))


if __name__ == '__main__':
    unittest.main()
