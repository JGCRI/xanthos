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

import xanthos.DataReader.DataLoad as DataLoad
from xanthos.DataWriter.OUTWriter import OUTWriter

import xanthos.pet.hargreaves as PETCalculation
import xanthos.runoff.hejazi as RG
import xanthos.routing.simple as SR

from xanthos.Diagnostics.Aggregation import Aggregation
from xanthos.Diagnostics.Diagnostics import Diagnostics
from xanthos.Diagnostics.TimeSeries import TimeSeriesPlot

from xanthos.accessible.accessible import AccessibleWater
from xanthos.hydropower.potential import HydropowerPotential
from xanthos.hydropower.actual import HydropowerActual


class Hydro:

    def __init__(self, settings):

        self.s = settings

        # climate data
        self.Precip = None
        self.Temp = None
        self.DTR = None
        self.sav_prev = None
        self.chs_prev = None

        # reference data
        self.GridConstants = None
        self.MSMC = None

        # pet
        self.latitude = None
        self.X = None
        self.M = None
        self.lambdaT = None
        self.dr = None

        # msmc
        self.Sm = None
        self.sav_prev = None

        # flow
        self.FlowDis_1D = None
        self.Ins_ChFlow = None
        self.ChVeloc_1D = None
        self.dsid = None
        self.upid = None
        self.UM = None

        # allocation
        self.PET = None
        self.AET = None
        self.Q = None
        self.Sav = None
        self.ChStorage = None
        self.Avg_ChFlow = None

        # simulation
        self.P = None
        self.T = None
        self.D = None
        self.Y = None
        self.DR = None
        self.mm = None
        self.deltat_routing = None

        # outputs
        self.q = None
        self.ac = None

    def load_climate(self):
        """
        Load climate data.
        @:return    Precipitation, Temperature, Daily Temperature Range arrays and settings.
        """
        # Load Climate Data
        t0 = time.time()

        self.Precip, self.Temp, self.DTR, self.sav_prev, self.chs_prev = DataLoad.load_gcm_data(self.s)

        t1 = time.time()
        print("---Load Climate Data: %s seconds ---" % (t1 - t0))

    def load_reference(self):
        """
        Load reference data.
        """
        t0 = time.time()

        # Loads grid constants (gridded maps) as dictionary
        self.GridConstants = DataLoad.load_map_data(self.s)

        # create a matrix (MSMC: Maximum Soil Moisture Capacity) with all data
        # 1: area; 2: region; 3: Max Soil Moisture (mm/month)
        self.MSMC = DataLoad.get_MaxSoilMoisture_matrix(self.GridConstants, self.s.ncell)

        t1 = time.time()
        print("---Load Gridded Map Data: %s seconds ---" % (t1 - t0))

    def setup_pet(self):

        # PET calculation related
        self.latitude = np.copy(self.GridConstants['Coord'][:, 2])

        # Latitude in radians
        self.X = self.latitude * np.pi / 180.

        # month-of-year and days per month for each month
        self.M = PETCalculation.set_month_arrays(self.s)

        self.lambdaT, self.dr = PETCalculation.calc_sinusoidal_factor(self.M)

    def setup_msmc(self):

        self.Sm = np.copy(self.MSMC[:, 2])

        if self.s.HistFlag == "True":
            self.sav_prev = 0.5 * self.Sm

    def setup_flow(self):

        self.FlowDis_1D = np.copy(self.GridConstants['FlowDis'])
        self.FlowDis_1D[np.where(self.FlowDis_1D < 1000)[0]] = 1000
        self.Ins_ChFlow = np.zeros((self.s.ncell,), dtype=float)
        self.ChVeloc_1D = np.ones((self.s.ncell,), dtype=float)

        # -1, 1-67420
        self.dsid = SR.downstream(self.GridConstants['Coord'], self.GridConstants['FlowDir'], self.s)

        # 1-67420
        self.upid = SR.upstream(self.GridConstants['Coord'], self.dsid, self.s)
        self.UM = SR.upstream_genmatrix(self.upid)

    def setup_alloc(self):

        self.PET = np.zeros_like(self.Precip)
        self.AET = np.zeros_like(self.Precip)
        self.Q = np.zeros_like(self.Precip)
        self.Sav = np.zeros_like(self.Precip)
        self.ChStorage = np.zeros_like(self.Precip)
        self.Avg_ChFlow = np.zeros_like(self.Precip)

    def prep_arrays(self, nm):
        """
        Prepare arrays
        :param nm:
        :param Precip:
        :param Temp:
        :param DTR:
        :param lambdaT:
        :param dr:
        :param M:
        :return:
        """
        self.P = np.copy(self.Precip[:, nm])  # keep nan in P
        self.T = np.nan_to_num(self.Temp[:, nm])  # nan to zero
        self.D = np.nan_to_num(self.DTR[:, nm])  # nan to zero
        self.Y = np.copy(self.lambdaT[nm])
        self.DR = np.copy(self.dr[nm])
        self.mm = np.copy(self.M[nm, 2])

    def staging(self):
        """
        Stage data for processing.
        """
        t0 = time.time()

        # load climate data
        self.load_climate()

        # load reference data
        self.load_reference()

        print("---Simulation of the Global Water Balance Model (GWBM): ")

        # PET calculation related
        self.setup_pet()

        # MSMC related data
        self.setup_msmc()

        # Flow related data
        self.setup_flow()

        # allocation
        self.setup_alloc()

        t1 = time.time()
        print("------Prepare simulation Data: %s seconds ------" % (t1 - t0))

    def calculate_pet(self):
        """
        Calculate monthly potential evapo-transpiration.
        """
        if self.s.pet.lower() == 'hargreaves':
            self.pet_t = PETCalculation.calculate_pet(self.T, self.D, self.X, self.Y, self.DR, self.mm)

        else:
            # use Hargreaves as default
            self.pet_t = PETCalculation.calculate_pet(self.T, self.D, self.X, self.Y, self.DR, self.mm)

    def calculate_runoff(self, nm):
        """
        Calculate runoff.
        """
        if self.s.runoff.lower() == 'hejazi':
            rg = RG.runoffgen(self.pet_t, self.P, self.T, self.D, self.s, self.Sm, self.X, self.Y, self.mm, self.DR, self.sav_prev)
            self.PET[:, nm], self.AET[:, nm], self.Q[:, nm], self.Sav[:, nm] = rg

        else:
            # use hejazi as default
            rg = RG.runoffgen(self.PET, self.P, self.T, self.D, self.s, self.Sm, self.X, self.Y, self.mm, self.DR, self.sav_prev)
            self.PET[:, nm], self.AET[:, nm], self.Q[:, nm], self.Sav[:, nm] = rg

    def calculate_routing(self, nm):
        """
        Calculate routing.
        """
        if self.s.routing.lower() == 'simple':
            # 3 hour time step for simplified solver
            self.deltat_routing = 3 * 3600

            sr = SR.streamrouting(self.FlowDis_1D, self.chs_prev, self.Ins_ChFlow, self.ChVeloc_1D,
                                  self.Q[:, nm], self.MSMC[:, 0], self.mm, self.deltat_routing, self.UM)

            self.ChStorage[:, nm], self.Avg_ChFlow[:, nm], self.Ins_ChFlow = sr

        else:
            # use simple as default

            # 3 hour time step for simplified solver
            self.deltat_routing = 3 * 3600

            sr = SR.streamrouting(self.FlowDis_1D, self.chs_prev, self.Ins_ChFlow, self.ChVeloc_1D,
                                  self.Q[:, nm], self.MSMC[:, 0], self.mm, self.deltat_routing, self.UM)

            self.ChStorage[:, nm], self.Avg_ChFlow[:, nm], self.Ins_ChFlow = sr

    def simulation(self, nm):
        """
        Run model simulation for a single time step.
        """

        if np.mod(nm, 12) == 0:
            print("Year: {}".format(self.M[nm, 0]))

        # set up data for processing
        self.prep_arrays(nm)

        # calculate pet
        self.calculate_pet()

        # calculate runoff and generate monthly potential ET, actual ET, runoff, and soil moisture
        self.calculate_runoff(nm)

        # channel storage, avg. channel flow (m^3/sec), instantaneous channel flow (m^3/sec)
        self.calculate_routing(nm)

        # update soil moisture (sav) and channel storage (chs) arrays for next step
        self.sav_prev = np.copy(self.Sav[:, nm])
        self.chs_prev = np.copy(self.ChStorage[:, nm])

    def spin_up(self):
        """
        Run spin-up for a user-defined number of years. If 0, then pass.
        """
        if self.s.SpinUp > 0:

            print("------Start the spin-up initialization:")
            t0 = time.time()

            for nm in range(int(self.s.SpinUp * 12)):
                self.simulation(nm)

            print("------Spin-up has finished successfully: %s seconds ------" % (time.time() - t0))

    def execute(self):
        """
        Run simulation for defined months.
        """
        print("------Start the simulation for all months: ")
        t0 = time.time()

        for nm in range(self.s.nmonths):
            self.simulation(nm)

        print("------Simulation has finished successfully: %s seconds ------" % (time.time() - t0))

    def accessible_water(self):
        """
        Run accessible water module
        """
        if self.s.CalculateAccessibleWater:
            print("---Start Accessible Water:")
            t0 = time.time()

            AccessibleWater(self.s, self.GridConstants, self.Q)

            print("---Accessible Water has finished successfully: %s seconds ------" % (time.time() - t0))

    def hydropower_potential(self):
        """
        Run hydropower potential module.
        """
        if self.s.CalculateHydropowerPotential:
            print("---Start Hydropower Potential:")
            t0 = time.time()

            HydropowerPotential(self.s, self.Avg_ChFlow)

            print("---Hydropower Potential has finished successfully: %s seconds ------" % (time.time() - t0))

    def hydropower_actual(self):
        """
        Run hydropower actual module.
        """
        if self.s.CalculateHydropowerActual:
            print("---Start Hydropower Actual:")
            t0 = time.time()

            HydropowerActual(self.s, self.Avg_ChFlow)

            print("---Hydropower Actual has finished successfully: %s seconds ------" % (time.time() - t0))

    def diagnostics(self):
        """
        Run diagnostics.
        """
        if self.s.PerformDiagnostics:
            print("---Start Diagnostics:")
            t0 = time.time()

            Diagnostics(self.s, self.Q, self.Avg_ChFlow, self.GridConstants)

            print("---Diagnostics has finished successfully: %s seconds ------" % (time.time() - t0))

    def output_simulation(self):
        """
        Output simulation results.
        """
        print "---Output simulation results:"
        t0 = time.time()

        self.q, self.ac = OUTWriter(self.s, self.GridConstants['Area'], self.PET, self.AET, self.Q, self.Sav, self.ChStorage, self.Avg_ChFlow)

        print("---Output finished: %s seconds ---" % (time.time() - t0))

    def aggregate_outputs(self):
        """
        Aggregation by Basin, Country, and/or Region.
        """
        if self.s.AggregateRunoffBasin > 0 or self.s.AggregateRunoffCountry > 0 or self.s.AggregateRunoffGCAMRegion > 0:
            print("---Start Aggregation:")
            t0 = time.time()

            Aggregation(self.s, self.GridConstants, self.q)

            print("---Aggregation has finished successfully: %s seconds ------" % (time.time() - t0))

    def plots(self):
        """
        Create time series plots.
        """
        if self.s.CreateTimeSeriesPlot:
            print("---Creating Time Series Plots:")
            t0 = time.time()

            TimeSeriesPlot(self.s, self.q, self.ac, self.GridConstants)

            print("---Plots has finished successfully: %s seconds ------" % (time.time() - t0))

    def process(self):
        """
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
        """
        t0 = time.time()
        print("# Start Hydrologic Model:")

        # set up data for run
        self.staging()

        # run spin-up
        self.spin_up()

        # execute simulation
        self.execute()

        # accessible water module
        self.accessible_water()

        # hydropower potential
        self.hydropower_potential()

        # hydropower actual
        self.hydropower_actual()

        # diagnostics
        self.diagnostics()

        # output simulation data
        self.output_simulation()

        # aggregate outputs
        self.aggregate_outputs()

        # create time series plots
        self.plots()

        print("# Hydrologic Model total time: %s seconds" % (time.time() - t0))
