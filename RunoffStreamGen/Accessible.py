'''
Created on March 1, 2017

Original:
@author: Chris R. Vernon's version (chris.vernon@pnnl.gov)

Rewritten:
@author: lixi729
@email: xinya.li@pnl.gov
@Project: Xanthos V1.0

'''

import numpy as np
import pandas as pd


def AccessibleWater(settings, GridConstants, runoff):
    bdf = GridConstants['BasinNames']
    # read in reservoir capacity at basin level
    rdf = pd.read_csv(settings.ResCapacityFile, header=None, names=['res_capacity'])
    # read in baseflow index (BFI) file
    bfi = pd.read_csv(settings.BfiFile)['bfi_avg']

    # read in hydrology model runoff data for future and historical combined
    # runoff: mm/month -> km3/year
    ny = int(settings.nmonths / 12)
    q = np.zeros((settings.ncell, ny), dtype=float)
    conversion = GridConstants['Area'] / 1e6  # mm -> km3
    for i in range(ny):
        q[:, i] = np.sum(runoff[:, i * 12:(i + 1) * 12], axis=1) * conversion

    # Basin Aggregation
    NM = settings.ncell
    Map = GridConstants['BasinIDs']
    NB = max(Map)
    Map_runoff = np.zeros((NB, ny), dtype=float)
    for y in range(0, ny):
        for index in range(0, NM):
            if not np.isnan(q[index, y]) and Map[index] > 0:
                Map_runoff[Map[index] - 1, y] += q[index, y]

                # calculate rolling mean for q data
    qs = RollingWindowFilter(Map_runoff, settings.MovingMeanWindow)
    q_gcam = QInGCAMYears(qs, settings)

    # calculate baseflow using qtot multiplied by BFI per basin
    bflow = np.transpose(np.transpose(q_gcam) * np.array(bfi))
    # calculate Environmental Flow Requirements (EFR) per basin using 10% of historical mean
    if settings.StartYear > settings.HistEndYear:
        print 'Warning! No historical data used in calculating Environmental Flow Requirements (EFR) per basin for Accessible Water'
        edf = settings.Env_FlowPercent * np.mean(Map_runoff, axis=1)
    elif settings.EndYear <= settings.HistEndYear:
        edf = settings.Env_FlowPercent * np.mean(Map_runoff, axis=1)
    else:
        ValidYears = range(settings.StartYear, settings.EndYear + 1)
        hey = ValidYears.index(settings.HistEndYear)
        edf = settings.Env_FlowPercent * np.mean(Map_runoff[:, :(hey + 1)], axis=1)

    # calculate available water for each basin and year
    ac = accessible_water(q_gcam, bflow, edf, rdf.as_matrix())

    # create accessible water output dataframe at GCAM year time step
    filename = settings.OutputFolder + "Accessible_Water_" + settings.OutputNameStr
    genGCAMOutput(filename, ac, bdf, settings)


def RollingWindowFilter(data, window, Dimension=0):
    # Obtain the moving average, the result is set to the center of the window (interval)
    # Dimension = 0 , column; = 1 , row
    # Window: odd integer, size of the interval
    # data: 1D or 2D data array

    weights = np.repeat(1.0, window) / window
    it = int((window - 1) / 2) + 1
    if data.ndim == 1:
        sma = np.convolve(data, weights, 'same')
    else:
        sma = np.zeros(data.shape, dtype=float)
        if Dimension == 1:  # moving average for each column
            for i in range(data.shape[1]):
                sma[:, i] = np.convolve(data[:, i], weights, 'same')
                sma[0, i] = np.mean(data[:it, i])
                sma[data.shape[0] - 1, i] = np.mean(data[data.shape[0] - it:, i])
        elif Dimension == 0:  # moving average for each row
            for i in range(data.shape[0]):
                sma[i, :] = np.convolve(data[i, :], weights, 'same')
                sma[i, 0] = np.mean(data[i, :it])
                sma[i, data.shape[1] - 1] = np.mean(data[i, data.shape[1] - it:])

    return sma


def QInGCAMYears(qs, settings):
    """
    Create data frame with only target GCAM years.
    """
    ValidYears = range(settings.StartYear, settings.EndYear + 1)
    GCAMYears = range(settings.GCAM_StartYear, settings.GCAM_EndYear + 1, settings.GCAM_YearStep)
    q_gcam = np.zeros((qs.shape[0], len(GCAMYears)), dtype=float)

    for i in range(len(GCAMYears)):
        q_gcam[:, i] = qs[:, ValidYears.index(GCAMYears[i])]

    return q_gcam


def accessible_water(qtot, base, efr, res):
    """
    Calculate accessible water.
    """
    ac = np.zeros(qtot.shape, dtype=float)
    for i in range(qtot.shape[1]):
        a = qtot[:, i] - efr
        b = base[:, i] - efr + res
        c = np.min(np.vstack((a, b)), axis=0)
        ac[:, i] = np.where(c < 0, 0, c)

    return ac


def genGCAMOutput(filename, data, bdf, settings):
    # Create data frame containing basin_id, basin_name, and accessible water by year.
    years = map(str, range(settings.GCAM_StartYear, settings.GCAM_EndYear + 1, settings.GCAM_YearStep))
    headerline = "ID,Name," + ",".join([year for year in years])

    maxID = len(bdf)
    MapId = np.arange(1, maxID + 1, 1, dtype=int).astype(str)
    newdata = np.insert(data.astype(str), 0, bdf, axis=1)
    Result = np.insert(newdata.astype(str), 0, MapId, axis=1)

    with open(filename + '.csv', 'w') as outfile:
        np.savetxt(outfile, Result, delimiter=',', header=headerline, fmt='%s')