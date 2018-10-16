"""
Calculate Monthly PET using the Hargreaves Method.

Rewritten:
@date: 10/07/2016
@author: lixi729
@Project: Xanthos V1.0

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute
"""

import numpy as np


def calculate_pet(temp, dtr, x, y, dr, m):
    """
    Calculate potential evapotranspiration for each month, mm/month.

    @:param temp:           temperature array (degrees Celsius)
    @:param dtr:            daily temperature range
    @:param x:              latitude in radians
    @:param y:              solar declination in radians
    @:param dr:             inverse relative distance Earth-Sun
    @:param m:              the number of days in the target month
    """
    ws = calc_daylight_hours(x, y)
    ra = calc_insolation(dr, ws, x, y)

    # treatment of negative values in D (temp range)
    dtr[np.where(dtr < 0)[0]] = 0.

    # factor of M converts this from daily to monthly volume
    evap = m * 0.0023 * ra * (temp + 17.8) * np.sqrt(dtr)
    _pet = np.maximum(evap, np.zeros_like(evap))
    pet = np.squeeze(_pet.T)

    return pet


def calc_insolation(dr, ws, X, Y):
    """Calculate Extraterrestrial solar radiation mm/day."""
    return 15.392 * dr * (ws * np.sin(X) * np.sin(Y) + np.cos(X) * np.cos(Y) * np.sin(ws))


def calc_daylight_hours(X, Y):
    """Calculate the number of daylight hours per day."""
    # sunset hour angle in radians
    return acos(-np.tan(X) * np.tan(Y))


def acos(x):
    """Inverse cosine of the elements of x.

    y = acos(x) returns the Inverse Cosine (cos-1) of the elements of x.
    For real values of x in the interval [-1,1], acos(x) returns real values in the interval [0,pi].
    For real values of x outside the interval [-1,1], acos(x) returns 0 for x > 1 and returns pi for x < -1
    """
    y = np.zeros_like(x)

    for i in range(x.shape[0]):

        if (x[i] <= 1) and (x[i] >= -1):
            y[i] = np.arccos(x[i])

        elif x[i] < -1:
            y[i] = np.arccos(-1)

        elif x[i] > 1:
            y[i] = np.arccos(1)

    return y
