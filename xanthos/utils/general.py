"""
General helper functions.

@author: lixi729
@Project: Xanthos V1.0

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute
"""

import numpy as np


def set_month_arrays(n_months, start_year, end_year):
    """
    Construct array for [ [year, month_index, num_days], ... ] reference.

    :param n_months:                number of months
    :param start_year:              start year
    :param end_year:                end year
    :return:                        numpy array of arrays [ [year, month_index, num_days], ... ]
    """
    counter = 0

    # year, month, number of days in month
    M = np.zeros((n_months, 3), dtype=int)

    # regular year days in month
    M1 = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    # leap year days in month
    M2 = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    for i in range(start_year, end_year + 1):
        # leap year
        if np.mod(i, 4) == 0:
            M0 = M2[:]

        # regular year
        else:
            M0 = M1[:]

        for j in range(12):
            M[counter, 0] = i
            M[counter, 1] = j
            M[counter, 2] = M0[j]
            counter += 1

    return M


def calc_sinusoidal_factor(yr_imth_ndays, startmonth=1):
    """
    Calculate the solar declination and inverse relative distance Earth-Sun.

    @:param yr_imth_ndays:   array of [ [year, month_index, num_days], ... ]
    @:return dr:             inverse relative distance Earth-Sun
    @:return solar_dec:      solar declination in radians
    """
    solar_dec = np.zeros((yr_imth_ndays.shape[0],), dtype=float)
    dr = np.zeros((yr_imth_ndays.shape[0],), dtype=float)

    for i in range(yr_imth_ndays.shape[0]):

        # leap year
        if np.mod(yr_imth_ndays[i, 0], 4) == 0:
            j = np.arange(1, 367)
            m1 = [1, 32, 61, 92, 122, 153, 183, 214, 245, 275, 306, 336]
            m2 = [31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366]

        # regular year
        else:
            j = np.arange(1, 366)

            # Month start days and end days
            m1 = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]
            m2 = [31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]

        # correction for dat aets that do not start on a year boundary
        ph = (startmonth - 1.) / 12. * 2. * np.pi
        lambda_temp = 0.409 * np.sin(2 * np.pi * j / max(j) - 1.39 + ph)
        dr_temp = 1. + 0.033 * np.cos(2 * np.pi * j / max(j) + ph)

        # monthly averages of lambda and dr
        mth = yr_imth_ndays[i, 1]
        solar_dec[i] = np.mean(lambda_temp[m1[mth] - 1:m2[mth]])
        dr[i] = np.mean(dr_temp[m1[mth] - 1:m2[mth]])

    return solar_dec, dr
