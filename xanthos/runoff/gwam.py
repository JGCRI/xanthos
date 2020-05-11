"""
Module for generating monthly runoff from climatological input data

@date   10/14/2016
@author: lixi729
@email: xinya.li@pnl.gov
@Project: Xanthos V1.0


License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute
"""

import numpy as np


def runoffgen(PET, P, settings, Sm, chstor, indexing=999):
    """
    Inputs:

    P:        monthly precip data from climate model
    T:        monthly temperature data from climate model
    D:        monthly average daily temperature range from climate model
    Sm:       (max) soil moisture
    X:        vector of gridpoint latitudes in radians
    Y:        monthly solar declination in radians
    dr:       monthly inverse relative distance Earth-Sun
    M:        days per month for each month in the dataset
    chstor:   channel storage (soil moisture), carried over from previous month

    Outputs: unit mm/month

    PET:               monthly potential evapo-transpiration
    AET:               monthly actual evapo-transpiration
    Q:                 monthly runoff volume
    Sav:               monthly soil column moisture
    """

    # Don't have an initial value for the soil moisture.
    # Initialize with half the maximum capacity instead.
    B = chstor + P - PET  # water in channel storage + precip net of PET.
    B = np.squeeze(B.T)

    # Initialize outputs
    Sav = np.zeros((settings.ncell,), dtype=float)
    Q = np.zeros((settings.ncell,), dtype=float)
    AET = np.zeros((settings.ncell,), dtype=float)

    # These are the branches of the conditional in the original
    c0 = np.where(Sm != 0)[0]
    c1 = np.where(Sm == indexing)[0]  # water bodies
    c1n = np.where(Sm != indexing)[0]
    tmp1 = np.where(B >= Sm)[0]  # B contains NaNs -> warning RuntimeWarning: invalid value encountered
    tmp2 = np.where(B < Sm)[0]  # B contains NaNs -> warning RuntimeWarning: invalid value encountered
    tmp3 = np.intersect1d(c0, c1n)
    c2 = np.intersect1d(tmp1, tmp3)
    c3 = np.intersect1d(tmp2, tmp3)

    # Condition 1: Lakes and such
    Q[c1] = np.maximum(0, P[c1] - PET[c1])  # NaNs in P-> warning RuntimeWarning: invalid value encountered
    Q[np.isnan(Q)] = 0.0  # treat NaNs from P
    AET[c1] = np.minimum(P[c1], PET[c1])  # NaNs in P -> warning RuntimeWarning: invalid value encountered
    AET[np.isnan(AET)] = PET[np.isnan(AET)]  # treat NaNs from P

    # Condition 2: Net precipitation is greater than max soil moisture
    Q[c2] = B[c2] - Sm[c2]
    Sav[c2] = Sm[c2]
    AET[c2] = PET[c2]

    # Condition 3: Net precipitation is less than max soil moisture
    alpha = 1
    tmp3 = chstor[c3] + P[c3]
    tmp5 = (5. * chstor[c3] / Sm[c3] - 2. * (chstor[c3] / Sm[c3]) ** 2.) / 3.
    tmp6 = np.minimum(np.ones_like(tmp5), tmp5)
    tmp7 = PET[c3] * np.maximum(0.1 * np.ones_like(tmp6), tmp6)
    AET[c3] = np.minimum(tmp3, tmp7)

    tmp8 = chstor[c3] * (1 - np.exp(-alpha * chstor[c3] / Sm[c3])) / (1 - np.exp(-alpha)) + (P[c3] - AET[c3])
    Sav[c3] = np.minimum(Sm[c3], tmp8)

    # Need special treatment if Sav turns out to be zero.
    subcond3a = np.intersect1d(c3, np.where(Sav <= 0)[0])
    Sav[subcond3a] = 0
    AET[subcond3a] = P[subcond3a] + chstor[subcond3a]
    Q[c3] = np.maximum(np.zeros_like(Q[c3]), chstor[c3] + P[c3] - AET[c3] - Sav[c3])

    return [PET, AET, Q, Sav]
