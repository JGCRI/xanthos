"""
Generate time series plots.

Created on Feb 10, 2017
@author: lixi729
@Project: Xanthos V1.0

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute
"""

import os
import logging
import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates


def TimeSeriesPlot(settings, Q, Avg_ChFlow, ref):
    """Create time series plots of Xanthos output data."""
    if not settings.CreateTimeSeriesPlot:
        return

    Folder = os.path.join(settings.OutputFolder, 'TimeSeriesPlot')

    if not os.path.exists(Folder):
        os.makedirs(Folder)

    x = {}
    # Prepare the data
    if settings.OutputInYear == 1:  # convert the original unit month to new unit year
        TimeUnit = 'year'
        x['data'] = np.array([datetime.datetime(i, 1, 1) for i in range(settings.StartYear, settings.EndYear + 1)])
        x['xmin'] = datetime.datetime(settings.StartYear - 1, 1, 1)
        x['xmax'] = datetime.datetime(settings.EndYear + 1, 1, 1)
    else:
        TimeUnit = 'month'
        x['data'] = np.array(
            [datetime.datetime(i, j, 1) for i in range(settings.StartYear, settings.EndYear + 1) for j in
             range(1, 13)])
        x['xmin'] = datetime.datetime(settings.StartYear - 1, 12, 1)
        x['xmax'] = datetime.datetime(settings.EndYear + 1, 1, 1)

    if settings.OutputUnit == 1:
        LengthUnit = 'km^3'
    else:
        LengthUnit = 'mm'

    if settings.TimeSeriesScale == 1:
        scalestr = 'Basin'
        CreateData_TimeSeriesScale(settings, Q, Avg_ChFlow, ref, scalestr, TimeUnit, LengthUnit, x)

    elif settings.TimeSeriesScale == 2:
        scalestr = 'Country'
        CreateData_TimeSeriesScale(settings, Q, Avg_ChFlow, ref, scalestr, TimeUnit, LengthUnit, x)

    elif settings.TimeSeriesScale == 3:
        scalestr = 'GCAMRegion'
        CreateData_TimeSeriesScale(settings, Q, Avg_ChFlow, ref, scalestr, TimeUnit, LengthUnit, x)
    else:
        scalestr = 'Basin'
        CreateData_TimeSeriesScale(settings, Q, Avg_ChFlow, ref, scalestr, TimeUnit, LengthUnit, x)
        scalestr = 'Country'
        CreateData_TimeSeriesScale(settings, Q, Avg_ChFlow, ref, scalestr, TimeUnit, LengthUnit, x)
        scalestr = 'GCAMRegion'
        CreateData_TimeSeriesScale(settings, Q, Avg_ChFlow, ref, scalestr, TimeUnit, LengthUnit, x)


def CreateData_TimeSeriesScale(settings, Q, Avg_ChFlow, ref, scalestr, TimeUnit, LengthUnit, x):

    if scalestr == 'Basin':
        id_fle = ref.basin_ids
        nm_fle = ref.basin_names
    elif scalestr == 'Country':
        id_fle = ref.country_ids
        nm_fle = ref.country_names
    elif scalestr == 'GCAMRegion':
        id_fle = ref.region_ids
        nm_fle = ref.region_names

    q = Aggregation_Map(id_fle, Q)
    ac = Aggregation_Map(id_fle, Avg_ChFlow)
    Names = np.insert(nm_fle, 0, 'Global')

    Folder = os.path.join(settings.OutputFolder, 'TimeSeriesPlot', '{}'.format(scalestr))

    if not os.path.exists(Folder):
        os.makedirs(Folder)

    q = np.insert(q, 0, np.sum(q, axis=0), axis=0)  # add global
    ac = np.insert(ac, 0, np.sum(ac, axis=0), axis=0)  # add global

    try:
        for i in settings.TimeSeriesMapID:

            outputname = os.path.join(Folder, '{0}{1}_{2}'.format(scalestr, i, Names[i]))

            Plot_TS(q[i, :], outputname, 'runoff', TimeUnit, LengthUnit, x)
            Plot_TS(ac[i, :], outputname, 'streamflow', 'sec', 'm^3', x)
            logging.info("Scale: {}, Create plots for {}_{}".format(scalestr, i, Names[i]))

    except:

        if settings.TimeSeriesMapID == 999:
            for i in range(0, q.shape[0]):

                outputname = os.path.join(Folder, '{0}{1}_{2}'.format(scalestr, i, Names[i]))

                Plot_TS(q[i, :], outputname, 'runoff', TimeUnit, LengthUnit, x)
                Plot_TS(ac[i, :], outputname, 'streamflow', 'sec', 'm^3', x)
            logging.info("Scale: {}, Create plots for all".format(scalestr))
        else:
            i = settings.TimeSeriesMapID

            outputname = os.path.join(Folder, '{0}{1}_{2}'.format(scalestr, i, Names[i]))

            Plot_TS(q[i, :], outputname, 'runoff', TimeUnit, LengthUnit, x)
            Plot_TS(ac[i, :], outputname, 'streamflow', 'sec', 'm^3', x)
            logging.info("Scale: {}, Create plots for {}_{}".format(scalestr, i, Names[i]))

    return


def Aggregation_Map(Map, runoff):

    NY = runoff.shape[1]
    NM = runoff.shape[0]
    NB = max(Map)
    Map_runoff = np.zeros((NB, NY), dtype=float)

    for y in range(0, NY):
        for index in range(0, NM):
            if not np.isnan(runoff[index, y]) and Map[index] > 0:
                Map_runoff[Map[index] - 1, y] += runoff[index, y]

    return Map_runoff


def Plot_TS(data, outputname, qstr, TimeUnit, LengthUnit, X):

    years = mdates.YearLocator()  # every year
    months = mdates.MonthLocator()  # every month
    yearsFmt = mdates.DateFormatter('%Y')

    fig = plt.figure()
    ax = plt.gca()
    ax.plot(X['data'], data)
    ax.xaxis.set_major_locator(years)
    ax.xaxis.set_major_formatter(yearsFmt)

    if TimeUnit == 'month':
        ax.xaxis.set_minor_locator(months)

    ax.set_xlim(X['xmin'], X['xmax'])
    ax.grid(True)
    plt.xlabel('Time (' + TimeUnit + ')', fontsize=12)
    plt.ylabel(qstr + ' ($' + LengthUnit + '$/' + TimeUnit + ')', fontsize=12)
    fig.autofmt_xdate()
    fig.savefig('{0}_{1}.png'.format(outputname, qstr), dpi=300)
    plt.close(fig)
