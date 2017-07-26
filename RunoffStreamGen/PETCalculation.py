'''
Calculating Monthly PET (Hargreaves Method)

Rewritten:
@date: 10/07/2016
@author: lixi729
@Project: Xanthos V1.0


License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute

'''

import numpy as np


def calculate_pet(T, D, X, Y, dr, M):
    """ Calculate potential evapotranspiration for each month, mm/month

    T, D, X, Y, dr, M

    """
    Ws = calc_daylight_hours(X, Y)
    Ra = calc_insolation(dr, Ws, X, Y)

    D[np.where(D < 0)[0]] = 0.  # treatment of negative values in D (temp range)

    evap = M * 0.0023 * Ra * (T + 17.8) * np.sqrt(D)  # Factor of M converts this from daily to monthly volume
    PET = np.maximum(evap, np.zeros_like(evap))
    PET = np.squeeze(PET.T)

    return PET


def calc_insolation(dr, ws, X, Y):
    """Calculate Extraterrestrial solar radiation mm/day """

    # R = 24*60*0.408*0.0820/np.pi * dr * (ws * np.sin(X) * np.sin(Y) + np.cos(X) * np.cos(Y) * np.sin(ws))
    R = 15.392 * dr * (ws * np.sin(X) * np.sin(Y) + np.cos(X) * np.cos(Y) * np.sin(ws))

    return R


def calc_daylight_hours(X, Y):
    """Calculate the number of daylight hours per day."""

    ws = acos(-np.tan(X) * np.tan(Y))  # sunset hour angle in radians
    # N = 24. / np.pi * ws

    return ws


def calc_sinusoidal_factor(M, startmonth=1):
    """Calculate the sinusoidal factor
    dr:      inverse relative distance Earth-Sun
    lambdaT: solar declination in radians
    """

    lambdaT = np.zeros((M.shape[0],), dtype=float)
    dr = np.zeros((M.shape[0],), dtype=float)

    for i in range(M.shape[0]):

        if np.mod(M[i, 0], 4) == 0:  # leap year
            j = np.arange(1, 367)
            m1 = [1, 32, 61, 92, 122, 153, 183, 214, 245, 275, 306, 336]
            m2 = [31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366]
        else:  # regular year
            j = np.arange(1, 366)
            m1 = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]  # Month start days and end days
            m2 = [31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]

        ph = (startmonth - 1.) / 12. * 2. * np.pi  # Correction for datasets that don't start on a year boundary.
        lambda_temp = 0.409 * np.sin(2 * np.pi * j / max(j) - 1.39 + ph)  #
        dr_temp = 1. + 0.033 * np.cos(2 * np.pi * j / max(j) + ph)

        # monthly averages of lambda and dr
        m = M[i, 1]
        lambdaT[i] = np.mean(lambda_temp[m1[m] - 1:m2[m]])
        dr[i] = np.mean(dr_temp[m1[m] - 1:m2[m]])

    return lambdaT, dr


def acos(x):
    # y = acos(x) returns the Inverse Cosine (cos-1) of the elements of x.
    # For real values of x in the interval [-1,1], acos(x) returns real values in the interval [0,pi].
    # For real values of x outside the interval [-1,1], acos(x) returns 0 for x > 1 and returns pi for x < -1

    y = np.zeros_like(x)
    for i in range(x.shape[0]):
        if x[i] <= 1 and x[i] >= -1:
            y[i] = np.arccos(x[i])
        elif x[i] < -1:
            y[i] = np.arccos(-1)
        elif x[i] > 1:
            y[i] = np.arccos(1)

    return y


def set_month_arrays(settings):
    counter = 0
    nmonth = settings.nmonths
    StartYear = settings.StartYear
    EndYear = settings.EndYear

    M = np.zeros((nmonth, 3), dtype=int)  # year, month, number of days in month
    M1 = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    M2 = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    for i in range(StartYear, EndYear + 1):
        if np.mod(i, 4) == 0:  # leap year
            M0 = M2[:]
        else:  # regular year
            M0 = M1[:]
        for j in range(12):
            M[counter, 0] = i
            M[counter, 1] = j
            M[counter, 2] = M0[j]
            counter += 1
    return M