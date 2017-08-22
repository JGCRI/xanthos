"""
@date     10/14/2016
@author   lixi729
@email:   xinya.li@pnl.gov
@Project: Xanthos V1.0
@Function: Run model to create runoff and streamrouting results

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute
"""
import numpy as np
import time

import DataReader.DataLoad as DataLoad
from DataWriter.OUTWriter import OUTWriter
from RunoffStreamGen import PETCalculation
from RunoffStreamGen import RunoffGen as RG
from RunoffStreamGen import StreamRouting as SR
from Diagnostics.Aggregation import Aggregation
from Diagnostics.Diagnostics import Diagnostics
from Diagnostics.TimeSeries import TimeSeriesPlot
from RunoffStreamGen.Accessible import AccessibleWater


def Hydro(settings):
    '''
    Inputs:
    precip:      monthly precipitation averages [ncell x nmonth]
    temp:        monthly temperature averages [ncell x nmonth]
    dtr:         monthly average daily temperature range [ncell x nmonth]
    initstorage: initial channel storage
    unitflag:    1= output in km^3/yr.  0= output in mm/yr
    histflag:    1= running historical hydro; 0= running future hydro
    startmonth:  month of year for the first month in the dataset.  1=Jan, 2=Feb, etc.

    Outputs:
    q:           monthly runoff [ncell x nmonth]
    pet:         monthly potential evapo-transpiration
    aet:         monthly actual evapo-transpiration
    s:           monthly channel storage
    Avg_ChFlow:  monthly average stream flow [ncell x nmonth]
    ChStorage:   month-end channel storage [ncell x nmonth]
    '''

    t0 = time.time()
    print "# Start Hydrologic Model:"
    # Load Climate Data
    t1 = time.time()
    #
    Precip, Temp, DTR, sav_prev, chs_prev = DataLoad.load_gcm_data(settings)  # original mat file data
    #
    t2 = time.time()
    print("---Load Climate Data: %s seconds ---" % (t2 - t1))

    # Loads grid constants (gridded maps) as dictionary
    #
    GridConstants = DataLoad.load_map_data(settings)
    # create a matrix (MSMC: Maximum Soil Moisture Capacity) with all data
    # 1: area; 2: region# ; 3: Max Soil Moisture (mm/month)
    #
    MSMC = DataLoad.get_MaxSoilMoisture_matrix(GridConstants, settings.ncell)
    #
    t1 = time.time()
    print("---Load Gridded Map Data: %s seconds ---" % (t1 - t2))

    print "---Simulation of the Global Water Balance Model (GWBM):"

    # PET calculation related
    latitude = np.copy(GridConstants['Coord'][:, 2])
    X = latitude * np.pi / 180.  # Latitude in radians
    M = PETCalculation.set_month_arrays(settings)  # month-of-year and days per month for each month
    lambdaT, dr = PETCalculation.calc_sinusoidal_factor(M)

    # MSMC related data
    Sm = np.copy(MSMC[:, 2])
    if settings.HistFlag == "True":
        sav_prev = 0.5 * Sm

    # Flow related data
    FlowDis_1D = np.copy(GridConstants['FlowDis'])
    FlowDis_1D[np.where(FlowDis_1D < 1000)[0]] = 1000
    Ins_ChFlow = np.zeros((settings.ncell,), dtype=float)
    ChVeloc_1D = np.ones((settings.ncell,), dtype=float)
    dsid = SR.downstream(GridConstants['Coord'], GridConstants['FlowDir'], settings)  # -1, 1-67420
    upid = SR.upstream(GridConstants['Coord'], dsid, settings)  # 1-67420
    UM = SR.upstream_genmatrix(upid)

    # Preallocation
    PET = np.zeros_like(Precip)
    AET = np.zeros_like(Precip)
    Q = np.zeros_like(Precip)
    Sav = np.zeros_like(Precip)
    ChStorage = np.zeros_like(Precip)
    Avg_ChFlow = np.zeros_like(Precip)

    t2 = time.time()
    print("------Prepare simulation Data: %s seconds ------" % (t2 - t1))

    if settings.SpinUp > 0:
        print "------Start the spin-up initialization:"
        # start the spin-up loop for initialization
        for nm in range(int(settings.SpinUp * 12)):
            if np.mod(nm, 12) == 0: print "Year", M[nm, 0]
            P = np.copy(Precip[:, nm])  # keep nan in P
            T = np.nan_to_num(Temp[:, nm])  # nan to zero
            D = np.nan_to_num(DTR[:, nm])  # nan to zero
            Y = np.copy(lambdaT[nm])
            DR = np.copy(dr[nm])
            mm = np.copy(M[nm, 2])
            # monthly potential ET, actual ET, runoff, soil moisture
            PET[:, nm], AET[:, nm], Q[:, nm], Sav[:, nm] = RG.runoffgen(P, T, D, settings, Sm, X, Y, mm, DR, sav_prev)
            sav_prev = np.copy(Sav[:, nm])

            deltat_routing = 3 * 3600  # 3 hour time step for simplified solver
            # channel storage, avg. channel flow (m^3/sec), instantaneous channel flow (m^3/sec)
            ChStorage[:, nm], Avg_ChFlow[:, nm], Ins_ChFlow = SR.streamrouting(FlowDis_1D, chs_prev, Ins_ChFlow,
                                                                               ChVeloc_1D, Q[:, nm], MSMC[:, 0], mm,
                                                                               deltat_routing, UM);
            chs_prev = np.copy(ChStorage[:, nm])

        t1 = time.time()
        print("------Spin-up has finished successfully: %s seconds ------" % (t1 - t2))

    print "------Start the simulation for all months: "
    # Start the main loop for all months
    for nm in range(settings.nmonths):
        if np.mod(nm, 12) == 0: print "Year", M[nm, 0]

        # P      = np.nan_to_num(Precip[:,nm])
        P = np.copy(Precip[:, nm])  # keep nan in P
        T = np.nan_to_num(Temp[:, nm])  # nan to zero
        D = np.nan_to_num(DTR[:, nm])  # nan to zero
        Y = np.copy(lambdaT[nm])
        DR = np.copy(dr[nm])
        mm = np.copy(M[nm, 2])
        # monthly potential ET, actual ET, runoff, soil moisture
        PET[:, nm], AET[:, nm], Q[:, nm], Sav[:, nm] = RG.runoffgen(P, T, D, settings, Sm, X, Y, mm, DR, sav_prev)
        sav_prev = np.copy(Sav[:, nm])

        deltat_routing = 3 * 3600  # 3 hour time step for simplified solver
        # channel storage, avg. channel flow (m^3/sec), instantaneous channel flow (m^3/sec)
        ChStorage[:, nm], Avg_ChFlow[:, nm], Ins_ChFlow = SR.streamrouting(FlowDis_1D, chs_prev, Ins_ChFlow, ChVeloc_1D,
                                                                           Q[:, nm], MSMC[:, 0], mm, deltat_routing,
                                                                           UM);
        chs_prev = np.copy(ChStorage[:, nm])

    t1 = time.time()
    print("------Simulation has finished successfully: %s seconds ------" % (t1 - t2))

    # AssessibleWater
    if settings.CalculateAccessibleWater:
        t3 = time.time()
        print("---Start Accessible Water:")
        AccessibleWater(settings, GridConstants, Q)
        t4 = time.time()
        print("---Accessible Water has finished successfully: %s seconds ------" % (t4 - t3))


    # Diagnostics
    if settings.PerformDiagnostics:
        t3 = time.time()
        print("---Start Diagnostics:")
        Diagnostics(settings, Q, Avg_ChFlow, GridConstants)
        t4 = time.time()
        print("---Diagnostics has finished successfully: %s seconds ------" % (t4 - t3))


    # Output simulation results
    print "---Output simulation results:"
    q, ac = OUTWriter(settings, GridConstants['Area'], PET, AET, Q, Sav, ChStorage, Avg_ChFlow)
    t2 = time.time()
    print("---Output finished: %s seconds ---" % (t2 - t1))

    # Aggregation by Basin/Country/GCAMRegion...
    if settings.AggregateRunoffBasin > 0 or settings.AggregateRunoffCountry > 0 or settings.AggregateRunoffGCAMRegion > 0:
        t3 = time.time()
        print("---Start Aggregation:")
        Aggregation(settings, GridConstants, q)
        t4 = time.time()
        print("---Aggregation has finished successfully: %s seconds ------" % (t4 - t3))

    # Create time series plots
    if settings.CreateTimeSeriesPlot:
        t3 = time.time()
        print("---Creating Time Series Plots:")
        TimeSeriesPlot(settings, q, ac, GridConstants)
        t4 = time.time()
        print("---Plots has finished successfully: %s seconds ------" % (t4 - t3))
        print("# Hydrologic Model total time: %s seconds" % (t2 - t0 + t4 - t3))
    else:
        print("# Hydrologic Model total time: %s seconds" % (t2 - t0))





