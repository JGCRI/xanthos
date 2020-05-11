"""
Calculating Monthly PET using the Thornthwaite Method.

Rewritten:
@date: 9/7/18
@author: Caleb Braun (caleb.braun@pnnl.gov)
@Project: Xanthos V2.0

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2018, Battelle Memorial Institute
"""

import calendar
import numpy as np


def calc_daylight_hours(mth_days, lat_radians):
    """
    Calculate average day length (hours) for each month at each latitude.

    @:param mth_days:       Iterable containing number of days for each month
    @:param lat_radians:    Numpy array of latitudes in radians
    @:return:               Numpy array (lat_radians x months)
    """
    EARTH_OBLIQ = 0.409  # Earth's obliquitiy in radians
    HOURS_IN_DY = 24.0
    NDAYS_IN_YR = np.arange(sum(mth_days)) + 1

    # Solar declination
    solar_dec = EARTH_OBLIQ * np.sin(((2 * np.pi / 365.0) * NDAYS_IN_YR - 1.39))

    # Sunset hour angle (domain of arccos is -1 <= x <= 1 radians; in this case,
    # greater or smaller values represent 24 hours sunlight/darkness)
    sunset_hour_angle_cos = -np.tan(lat_radians[:, np.newaxis]) * np.tan(solar_dec[np.newaxis, :])
    sunset_hour_angle = np.arccos(np.clip(sunset_hour_angle_cos, -1, 1))

    # Sum daylight hours for each month, then divde by month length for mean
    daylight_hours = sunset_hour_angle * (HOURS_IN_DY / np.pi)
    month_indices = np.roll(np.cumsum(mth_days), 1)
    month_indices[0] = 0
    mthly_daylight_hours = np.add.reduceat(daylight_hours, month_indices, axis=1) / mth_days

    return mthly_daylight_hours


def execute(tas, lat_radians, start_yr, end_yr):
    """
    Estimate potential evapotranspiration (PET) using the Thornthwaite method.

    Thornthwaite equation:

        PET = 1.6 (L/12) (N/30) (10Ta / I)^a

    where, for each month:

    L  = the mean day length in hours of the month being calculated
    N  = the number of days in the month being calculated
    Ta = the mean daily air temperature in deg C (if negative use 0) of the
         month being calculated
    I  = a heat index which depends on the 12 monthly mean temperatures and
         is calculated as the sum of (Ta_i / 5)^1.514 for each month i, where
         Ta_i is the air temperature for each month in the year
    a  = (6.75 x 10-7) * I^3 - (7.71 x 10-5) * I^2 + (1.792 x 10-2) * I + 0.49239

    For more information, see Thornthwaite (1948), available here:
        https://doi.org/10.2307/210739

    :param tas:             Numpy array of mean daily air temperature in deg C
    :param lat_radians:     Vector of latitude in radians for each cell in tas
    :param start_yr:        Data start year
    :param end_yr:          Data end year

    :return:    Estimated monthly potential evaporation of each month of the
                year in mm/month
    """
    MONTHDAYS = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    LEAP_MONTHDAYS = (31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    NMONTHS = 12

    # All nans and negatives are treated as zero
    tas[np.logical_or(np.isnan(tas), tas < 0)] = 0

    # --- Heat Index and exponential calculations:
    i = np.copy(tas)

    # Monthly Thornthwaite Heat Index formula:
    i = np.power((i / 5.0), 1.514)

    # Aggregate to Annual Heat Index
    I = np.add.reduceat(i, np.arange(0, i.shape[1], NMONTHS), axis=1)

    # Calculate the exponential constants
    a = (.000000675 * I ** 3) - (.0000771 * I ** 2) + (.0179 * I) + .492

    # Spread yearly values back over all months
    I = np.repeat(I, NMONTHS, axis=1)
    a = np.repeat(a, NMONTHS, axis=1)

    # --- Pet calculations:
    # Calculate PET before adjusting for real month length and theoretical
    # sunshine hours.
    # If you don't pass `out` the indices where (I == 0) will be uninitialized.
    pet_unadj = np.divide(10 * tas, I, out=np.zeros_like(I), where=(I != 0))
    pet_unadj = 16 * np.power(pet_unadj, a)

    # --- Correction variables
    # Monthly mean daylight hours
    L = calc_daylight_hours(MONTHDAYS, lat_radians)
    L = np.repeat(L, (end_yr - start_yr + 1), axis=1)

    # Determine which years in the data are leap years
    leapyears = [calendar.isleap(y) for y in range(start_yr, end_yr + 1)]

    # Replace any leap years with leap year mean daylight hours
    if any(leapyears):
        L_leap = calc_daylight_hours(LEAP_MONTHDAYS, lat_radians)
        leap_months = np.repeat(leapyears, NMONTHS)
        L[:, leap_months] = np.tile(L_leap, sum(leapyears))

    # 1d array of the number of days in each month
    N = np.array([LEAP_MONTHDAYS if ly else MONTHDAYS for ly in leapyears]).flatten()

    # Final equation
    #   L = Monthly mean daylight hours  (shape = pet_unadj.shape)
    #   N = Number of days in each month (shape = nmonths)
    pet = pet_unadj * (L / NMONTHS) * (N / 30.0)

    # return numpy array of PET in mm/month [ncells x nmonths]
    return pet
