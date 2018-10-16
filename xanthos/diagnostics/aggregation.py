"""
Aggregation functions.

Created on Oct 26, 2016
@author: lixi729
@Project: Xanthos V1.0

1. The option to aggregate the gridded results by the river basins / country / region
2. Output in csv files, total runoff by basin/country/region... (row) and month/year (column)
   unit is determined by user settings

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute
"""

import numpy as np
import os
import logging
import pandas as pd


def Aggregation(settings, ref, q):
    """Aggregate and write results based on user settings."""
    Aggregation = {}

    if (settings.AggregateRunoffBasin > 0 or
            settings.AggregateRunoffCountry > 0 or
            settings.AggregateRunoffGCAMRegion > 0):

        if settings.AggregateRunoffBasin > 0:  # Basin
            logging.info("Aggregating by Basin")
            Aggregation['Basin_runoff'] = Aggregation_Map(settings, ref.basin_ids,
                                                          ref.basin_names, q, "Basin_runoff")
            logging.info("Basin_runoff: unit is {}".format(settings.OutputUnitStr))

        if settings.AggregateRunoffCountry > 0:  # Country
            logging.info("Aggregating by Country")
            Aggregation['Country_runoff'] = Aggregation_Map(settings, ref.country_ids, ref.country_names, q,
                                                            "Country_runoff")
            logging.info("Country_runoff: unit is {}".format(settings.OutputUnitStr))

        if settings.AggregateRunoffGCAMRegion > 0:  # GCAMRegion
            logging.info("Aggregating by GCAM Region")
            Aggregation['Region_runoff'] = Aggregation_Map(settings, ref.region_ids, ref.region_names, q,
                                                           "GCAMRegion_runoff")
            logging.info("Country_runoff: unit is ", settings.OutputUnitStr)

    return Aggregation


def Aggregation_Map(settings, Map, Names, runoff, varstr):
    """
    Aggregate runoff by basin/country/region.

    Aggregate monthly values to annual and output as csv.

    :return:    table of total runoff by basin/country/region (row) and year (column)
    """
    if settings.OutputInYear == 1:
        NT = int(settings.nmonths / 12)
    else:
        NT = int(settings.nmonths)

    NM = settings.ncell
    NB = max(Map)
    Map_runoff = np.zeros((NB, NT), dtype=float)

    for y in range(0, NT):
        for index in range(0, NM):
            if not np.isnan(runoff[index, y]) and Map[index] > 0:
                Map_runoff[Map[index] - 1, y] += runoff[index, y]

    # output
    maxID = max(Map)
    MapId = np.arange(1, maxID + 1, 1, dtype=int).astype(str)
    newdata = np.insert(Map_runoff.astype(str), 0, Names, axis=1)
    Result = np.insert(newdata.astype(str), 0, MapId, axis=1)

    writecsvAggregation(Result, settings, varstr)

    return Map_runoff


def writecsvAggregation(data, settings, var):
    """Save aggregate data as csv file."""
    filename = os.path.join(settings.OutputFolder,
                            '{}_{}_{}'.format(var, settings.OutputUnitStr, '_'.join(settings.ProjectName.split(' '))))

    # convert to data frame to set header and basin number in file
    df = pd.DataFrame(data)

    filename += ".csv"

    if settings.OutputInYear == 1:
        cols = ','.join(['{}'.format(i) for i in range(settings.StartYear, settings.EndYear + 1, 1)])
    else:
        l = []
        for i in range(settings.StartYear, settings.EndYear + 1, 1):
            for m in range(1, 13):
                if m < 10:
                    mth = '0{}'.format(m)
                else:
                    mth = m
                l.append('{}{}'.format(i, mth))
        cols = ','.join(l)

    # set header
    hdr = 'id,name,{}'.format(cols)
    df.columns = hdr.split(',')

    df.to_csv(filename, index=False)
