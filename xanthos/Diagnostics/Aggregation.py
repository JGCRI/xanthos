'''
Created on Oct 26, 2016
@author: lixi729
@Project: Xanthos V1.0

1. The option to aggregate the gridded results by the river basins / country / region
2. Output in csv files, total runoff by basin/country/region... (row) and month/year (column), unit is determined by user settings


License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute

'''

import numpy as np


def Aggregation(settings, GridConstants, q):
    Aggregation = {}

    if settings.AggregateRunoffBasin > 0 or settings.AggregateRunoffCountry > 0 or settings.AggregateRunoffGCAMRegion > 0:

        if settings.AggregateRunoffBasin > 0:  # Basin
            print("Aggregation by Basin")
            Aggregation['Basin_runoff'] = Aggregation_Map(settings, GridConstants['BasinIDs'],
                                                          GridConstants['BasinNames'], q, "Basin_runoff")        
            print "Basin_runoff: unit is ", settings.OutputUnitStr

        if settings.AggregateRunoffCountry > 0:  # Country
            print("Aggregation by Country")
            Aggregation['Country_runoff'] = Aggregation_Map(settings, GridConstants['CountryIDs'],
                                                            GridConstants['CountryNames'], q, "Country_runoff")
            print "Country_runoff: unit is ", settings.OutputUnitStr

        if settings.AggregateRunoffGCAMRegion > 0:  # GCAMRegion
            print("Aggregation by GCAM Region")
            Aggregation['Region_runoff'] = Aggregation_Map(settings, GridConstants['GCAMRegionIDs'],
                                                           GridConstants['GCAMRegionNames'], q, "GCAMRegion_runoff")
            print "Country_runoff: unit is ", settings.OutputUnitStr

    return Aggregation


def Aggregation_Map(settings, Map, Names, runoff, varstr):
    '''aggregate runoff by basin/country/region..., aggregate monthly values to annual.
    Return table of total runoff by basin/country/region... (row) and year (column) and output in csv'''

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
    filename = settings.OutputFolder + varstr + "_" + str(maxID) + "_" + settings.OutputNameStr
    writecsvAggregation(filename, Result, settings)

    return Map_runoff


def writecsvAggregation(filename, data, Settings):
    years = map(str, range(Settings.StartYear, Settings.EndYear + 1))
    if Settings.OutputInYear == 1:
        headerline = "ID, Name," + ",".join([year for year in years]) + ", Unit (" + Settings.OutputUnitStr + ")"
    else:
        MonthStr = np.chararray((len(years) * 12,), itemsize=6)
        for y in years:
            N = years.index(y)
            MonthStr[N * 12:(N + 1) * 12] = [str(y) + str(i).zfill(2) for i in range(1, 13)]
        headerline = "ID, Name," + ",".join([k for k in MonthStr]) + ", Unit (" + Settings.OutputUnitStr + ")"

    with open(filename + '.csv', 'w') as outfile:
        np.savetxt(outfile, data, delimiter=',', header=headerline, fmt='%s')