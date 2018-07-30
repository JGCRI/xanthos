import calendar
import numpy as np


def days_per_month(start_year, end_year):
    """
    Get days per month as list. Account for leap years.
    """
    l = []
    yrs = range(start_year, end_year + 1, 1)

    for yr in yrs:

        for mth in range(1, 13, 1):
            l.append(calendar.monthrange(yr, mth)[1])

    return l


def pet(lat, mth, t, tmax, tmin):

    # day for the middle of the month for each month in year
    j = [15, 45, 75, 105, 135, 165, 195, 225, 255, 285, 315, 345]

    # change index to month in sequence
    mth2 = (mth + 1) % 12
    if mth2 == 0:
        mth2 = 12

    j1 = j[mth2 - 1]

    delta = 0.4102 * np.sin(2 * (np.pi / 365) * (j1 - 80))

    phi = (lat * np.pi / 180)

    # extraterrestrial radiation
    ra = 118 / np.pi * (np.arccos(-np.tan(delta) * np.tan(phi))) + np.cos(phi) * np.cos(delta) * np.sin(np.arccos(-np.tan(delta) * np.tan(phi)))

    # calculate pet
    pet = 0.408 * 0.0023 * ra * (t + 17.8) * np.sqrt(tmax - tmin)

    if t < 0:
        pet = 0

    return pet


def execute(config, data, out_file=None):

    # n_cells, n_months
    hgpet = np.zeros(shape=(config.ncell, config.nmonths)) * np.nan

    # get days per month for time range
    mon2 = days_per_month(config.StartYear, config.EndYear)

    for mth_idx in range(0, config.nmonths, 1):

        yr_val = (mth_idx + 1) / 12

        if yr_val > 0:

            yr = config.StartYear + (yr_val - 1)

            print("\t\tProcessing Year: {}".format(yr))

        for cell in range(0, config.ncell, 1):
            hgpet[cell, mth_idx] = pet(data.coords[2][cell], mth_idx, data.hs_tas[cell, mth_idx], data.hs_tmax[cell, mth_idx], data.hs_tmin[cell, mth_idx])

    # mm/day to mm/month
    hgpet *= mon2

    if out_file is not None:
        np.save(out_file, hgpet)

    return hgpet