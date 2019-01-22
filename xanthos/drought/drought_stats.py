"""
Functions for computing drought statistics from hydrological model
output.
"""

import numpy as np


def droughtstats(hydroout, threshvals):
    """Compute Severity, Intensity, and Duration statistics for a matrix of hydrological output.

    :param hydroout: matrix[ntime x ngrid] of outputs for a hydrological
                     variable, such as runoff or soil moisture.
    :param threshvals: matrix[K x ngrid] of drought thresholds.  The K different
                     sets of thresholds will be recycled as necessary over all
                     the ntime time values.  Typically, either K==1 (single set
                     of thresholds for all times), or K==12 (different set of
                     thresholds for each month of the year).
    :return: (S, I, D) matrices of severity, intensity, and duration, each with
                     the same dimensions as hydroout.

    Duration is defined as the number of consecutive months that the
    hydrological variable has been under the drought threshold.  By definition,
    the duration of a grid cell that is not currently under drought conditions
    is zero.

    Severity is the cumulative sum since the start of the current drought
    episode of the amount by which the hydrological variable fell short of the
    drought threshold.  If q_0(t) is the (possibly time variable) threshold, and
    t_0 is the start of the current drought episode, then
    $$
       S(t) = \sum_{t'=t_0}^t q_0(t') - q(t').
    $$
    Severity is also defined to be zero for a grid cell that is not currently
    under drought conditions.

    Intensity is the average severity over a drought period.  I = S/D.  It is
    defined to be (you guessed it) zero for a grid cell not currently under
    drought conditions.
    """

    S = np.empty_like(hydroout)
    I = np.empty_like(hydroout)
    D = np.empty_like(hydroout)

    # Handle the first case specially, since there won't be a previous time
    # step.

    isdrought = hydroout[0, :] < threshvals[0, :]
    D[0, :] = np.where(isdrought, 1.0, 0.0)
    I[0, :] = S[0, :] = np.where(isdrought, (threshvals[0, :] - hydroout[0, :])/threshvals[0, :], 0.0)

    ntime = hydroout.shape[0]
    nthresh = threshvals.shape[0]

    for t in range(1, ntime):
        tm1 = t-1
        m = t % nthresh         # "month"
        thresh = threshvals[m, :]
        hydro = hydroout[t, :]

        dtm1 = D[tm1, :]
        stm1 = S[tm1, :]

        dt = D[t, :]
        st = S[t, :]
        it = I[t, :]

        isdrought = hydro < thresh
        dt = np.where(isdrought, dtm1+1, 0.0)
        st = np.where(isdrought, stm1 + (thresh - hydro)/thresh, 0.0)
        it = np.where(isdrought, st / dt, 0.0)

    return (S, I, D)


def getthresh(histout, nper, quantile=0.1):
    """
    Compute quantile-based thresholds using data from a reference period.

    :param histout: Matrix[ntime x ngrid] of output variable during the reference period.
    :param nper: Number of periods to calculate thresholds for.  Generally either
                 1 (single threshold for all periods), or 12 (thresholds by month)
    :param quantile: Quantile to use for determining the threshold.
    :return: Matrix[nper x ngrid] of threshold values.

    The array passed as histout should contain only values from the reference
    period, from which drought thresholds will be calculated.

    """

    (ntime, ngrid) = histout.shape
    # Check that ntime is a multiple of nper
    nyear = int(ntime / nper)

    histout = np.reshape(histout, (nyear, nper, ngrid))

    return np.quantile(histout, quantile, axis=0)
