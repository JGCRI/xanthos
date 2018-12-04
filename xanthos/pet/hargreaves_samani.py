"""
Hargreaves-Samani PET.

@author   Yaling Liu
@email:   cauliuyaling@gmail.com
@Project: Xanthos 2.0

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute
"""

import logging
import calendar
import numpy as np


def days_per_month(start_year, end_year):
    """Get days per month as list. Account for leap years."""
    l = []
    yrs = list(range(start_year, end_year + 1, 1))

    for yr in yrs:

        for mth in range(1, 13, 1):
            l.append(calendar.monthrange(yr, mth)[1])

    return l


def pet(lat, mth, t, tmax, tmin):

    if t < 0:
        return 0

    else:

        # day for the middle of the month for each month in year
        j = [15, 45, 75, 105, 135, 165, 195, 225, 255, 285, 315, 345]

        # change index to month in sequence
        if mth == 0:
            dy = j[0]

        else:
            dy = j[mth % 12]

        delta = 0.4102 * np.sin(2 * (np.pi / 365) * (dy - 80))

        phi = (lat * np.pi / 180)

        tn = -np.tan(delta) * np.tan(phi)

        if (tn < -1.) or (tn > 1.):
            acs = 0
        else:
            acs = np.arccos(tn)

        # extraterrestrial radiation
        ra = 118 / np.pi * acs + np.cos(phi) * np.cos(delta) * np.sin(acs)

        # calculate pet
        pet = 0.408 * 0.0023 * ra * (t + 17.8) * np.sqrt(abs(tmax - tmin))

        return pet


def update_user(mth_idx, start_year, prev_yr=None):

    msg = "\t\tProcessing Year: {}"

    yr_val = (mth_idx + 1) / 12.0

    if mth_idx == 0:

        logging.info(msg.format(start_year))
        return start_year

    if yr_val.is_integer():

        yr = prev_yr + 1

        logging.info(msg.format(yr))
        return yr

    else:

        return prev_yr


def execute(config, data, out_file=None):

    # n_cells, n_months
    hgpet = np.zeros(shape=(config.ncell, config.nmonths)) * np.nan

    # get days per month for time range
    mon2 = days_per_month(config.StartYear, config.EndYear)

    py = None
    for mth_idx in range(0, config.nmonths, 1):

        py = update_user(mth_idx, config.StartYear, py)

        for cell in range(0, config.ncell, 1):

            lat = data.coords[:, 2][cell]
            t = data.hs_tas[cell, mth_idx]
            tmax = data.hs_tmax[cell, mth_idx]
            tmin = data.hs_tmin[cell, mth_idx]

            hgpet[cell, mth_idx] = pet(lat, mth_idx, t, tmax, tmin)

    # mm/day to mm/month
    hgpet *= mon2

    if out_file is not None:
        np.save(out_file, hgpet)

    return hgpet
