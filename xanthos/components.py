"""
Components for use in model configurations.

@author   Chris R. Vernon , lixi729
@email:   chris.vernon@pnnl.gov; xinya.li@pnl.gov
@Project: Xanthos 2.0

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute
"""

import numpy as np
import time

import xanthos.data_reader.DataLoad as DataLoad
from xanthos.DataWriter.OUTWriter import OUTWriter
from xanthos.Diagnostics.Aggregation import Aggregation
from xanthos.Diagnostics.Diagnostics import Diagnostics
from xanthos.Diagnostics.TimeSeries import TimeSeriesPlot
from xanthos.accessible.accessible import AccessibleWater
from xanthos.hydropower.potential import HydropowerPotential
from xanthos.hydropower.actual import HydropowerActual

import xanthos.pet.hargreaves as pet_mod
import xanthos.routing.simple as routing_mod


class Components:
    """
    Components for use in model configurations.

    @author   Chris R. Vernon , lixi729
    @email:   chris.vernon@pnnl.gov; xinya.li@pnl.gov
    @Project: Xanthos 2.0

    License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

    Copyright (c) 2017, Battelle Memorial Institute
    """

    def __init__(self, config):

        self.s = config

        # climate data
        self.Precip = None
        self.Temp = None
        self.DTR = None

        # reference data
        self.GridConstants = None

        # pet
        self.latitude = None
        self.X = None
        self.M = None
        self.lambdaT = None
        self.dr = None
        self.pet_t = None
        self.pet_out = None

        # soil moisture content
        self.MSMC = None
        self.Sm = None
        self.grid_area = None

        # flow
        self.FlowDis_1D = None
        self.Ins_ChFlow = None
        self.ChVeloc_1D = None
        self.dsid = None
        self.upid = None
        self.UM = None

        # runoff feedback for soil moisture; initial values from historical file if future mode
        self.sav_prev = None

        # routing feedback for channel storage; initial values from historical file if future mode
        self.chs_prev = None

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

        self.Precip, self.Temp, self.DTR = DataLoad.load_gcm_data(self.s)

        t1 = time.time()
        print("---Load Climate Data: %s seconds ---" % (t1 - t0))

    def load_reference(self):
        """
        Load reference data.
        """
        t0 = time.time()

        # Loads grid constants (gridded maps) as dictionary
        self.GridConstants = DataLoad.load_map_data(self.s)

        t1 = time.time()
        print("---Load Gridded Map Data: %s seconds ---" % (t1 - t0))

    def setup_pet(self):

        # PET calculation related
        self.latitude = np.copy(self.GridConstants['Coord'][:, 2])

        # Latitude in radians
        self.X = self.latitude * np.pi / 180.

        # month-of-year and days per month for each month
        self.M = pet_mod.set_month_arrays(self.s.nmonths, self.s.StartYear, self.s.EndYear)

        self.lambdaT, self.dr = pet_mod.calc_sinusoidal_factor(self.M)

    def setup_msmc(self):
        """
        Assign max soil moisture (mm/month) [2] to Sm.  For historic data use 0.5 * sm to an initial value to pass to
        runoff model. If future mode, read values from historical file.
        """

        # create a matrix (MSMC: Maximum Soil Moisture Capacity) with all data
        # 1: area; 2: region; 3: Max Soil Moisture (mm/month)
        self.MSMC = DataLoad.get_MaxSoilMoisture_matrix(self.GridConstants, self.s.ncell)

        # extract grid area
        self.grid_area = np.copy(self.MSMC[:, 0])

        # soil moisture
        self.Sm = np.copy(self.MSMC[:, 2])

        # load historic soil moisture data if future or create zeros array if historic
        self.sav_prev = DataLoad.load_soil_data(self.s)
        self.chs_prev = DataLoad.load_chs_data(self.s)

        if self.s.HistFlag == "True":
            self.sav_prev = 0.5 * self.Sm

    def setup_flow(self):

        self.FlowDis_1D = np.copy(self.GridConstants['FlowDis'])
        self.FlowDis_1D[np.where(self.FlowDis_1D < 1000)[0]] = 1000

        self.Ins_ChFlow = np.zeros((self.s.ncell,), dtype=float)

        self.ChVeloc_1D = np.copy(self.GridConstants['ChVeloc'])
        self.ChVeloc_1D[np.where(self.ChVeloc_1D < 0)[0]] = 0

        # -1, 1-67420
        self.dsid = routing_mod.downstream(self.GridConstants['Coord'], self.GridConstants['FlowDir'], self.s)

        # 1-67420
        self.upid = routing_mod.upstream(self.GridConstants['Coord'], self.dsid, self.s)
        self.UM = routing_mod.upstream_genmatrix(self.upid)

    def setup_alloc(self):

        self.PET = np.zeros(shape=(self.s.ncell, self.s.nmonths))
        self.AET = np.zeros(shape=(self.s.ncell, self.s.nmonths))
        self.Q = np.zeros(shape=(self.s.ncell, self.s.nmonths))
        self.Sav = np.zeros(shape=(self.s.ncell, self.s.nmonths))
        self.ChStorage = np.zeros(shape=(self.s.ncell, self.s.nmonths))
        self.Avg_ChFlow = np.zeros(shape=(self.s.ncell, self.s.nmonths))

    def prep_arrays(self, nm=None):
        """
        Prepare arrays
        """
        if nm is None:
            self.P = np.copy(self.Precip)  # keep nan in P
            self.T = np.nan_to_num(self.Temp)  # nan to zero
            self.D = np.nan_to_num(self.DTR)  # nan to zero
            self.pet_out = np.zeros_like(self.Precip)

        else:
            self.P = np.copy(self.Precip[:, nm])  # keep nan in P
            self.T = np.nan_to_num(self.Temp[:, nm])  # nan to zero
            self.D = np.nan_to_num(self.DTR[:, nm])  # nan to zero
            self.pet_out = None

    def prep_pet(self, nm):
        """
        Prepare PET for processing.
        """
        self.Y = np.copy(self.lambdaT[nm])
        self.DR = np.copy(self.dr[nm])
        self.mm = np.copy(self.M[nm, 2])
        self.temp_pet = np.nan_to_num(self.Temp[:, nm])  # nan to zero
        self.dtr_pet = np.nan_to_num(self.DTR[:, nm])  # nan to zero

    def calculate_pet(self, archive=False):
        """
        Calculate monthly potential evapo-transpiration.
        """
        if self.s.pet_module == 'hargreaves':

            self.pet_t = pet_mod.calculate_pet(self.temp_pet, self.dtr_pet, self.X, self.Y, self.DR, self.mm)

    def calculate_runoff(self, nm=None):
        """
        Calculate runoff.

        Runoff takes a feedback of soil moisture content (sav_prev).  This is updated for each iteration of a month.
        """
        if self.s.runoff_module == 'hejazi':
            import xanthos.runoff.hejazi as runoff_mod

            rg = runoff_mod.runoffgen(self.pet_t, self.P, self.T, self.D, self.s, self.Sm, self.X, self.Y, self.mm, self.DR, self.sav_prev)
            self.PET[:, nm], self.AET[:, nm], self.Q[:, nm], self.Sav[:, nm] = rg

        elif self.s.runoff_module == 'abcd':
            import xanthos.runoff.abcd as runoff_mod

            # run all basins at once in parallel
            rg = runoff_mod.abcd_execute(n_basins=235, basin_ids=self.GridConstants['BasinIDs'],
                                         pet=self.pet_out, precip=self.P, tmin=self.T, calib_file=self.s.calib_file,
                                         out_dir=self.s.ro_out_dir, n_months=self.s.nmonths,
                                         spinup_factor=self.s.SpinUp, jobs=self.s.ro_jobs)

            self.PET, self.AET, self.Q, self.Sav = rg

    def calculate_routing(self, nm):
        """
        Calculate routing.  Routing takes a simulated runoff (Q) from the runoff output and
        previous channel storage (chs_prev) from previous channel storage.
        """
        if self.s.routing_module == 'simple':
            # 3 hour time step for simplified solver
            self.deltat_routing = 3 * 3600

            sr = routing_mod.streamrouting(self.FlowDis_1D, self.chs_prev, self.Ins_ChFlow, self.ChVeloc_1D,
                                           self.Q[:, nm], self.grid_area, self.mm, self.deltat_routing, self.UM)

            self.ChStorage[:, nm], self.Avg_ChFlow[:, nm], self.Ins_ChFlow = sr

    def simulation(self, num_steps=None, pet=True, runoff=True, runoff_step='month',
                   routing=True, routing_step='month', notify='simulation'):
        """
        Run model simulation for a defined configuration.

        :param num_steps:           The number of time steps to process (INT)
        :param pet:                 True if running PET, False if embedded in runoff model
        :param runoff:              True if running Runoff, False if not
        :param runoff_step:         The time unit as a string; if None the runoff model iterates internally, else 'month'
        :param routing:             True if running Routing, False if not
        :param routing_step:        The time unit as a string; if None the routing model iterates internally, else 'month'
        :param notify:              A string that is used to add to log print that describes whether the simulation is
                                    spin-up or regular
        """
        # pass simulation if there are no steps to process
        if num_steps == 0:
            pass

        else:
            if (runoff_step == 'month') and (routing_step == 'month'):

                print("---{} in progress...".format(notify))
                t0 = time.time()

                for nm in range(num_steps):

                    if np.mod(nm, 12) == 0:
                        print("Processing Year: {}".format(self.M[nm, 0]))

                    # set up climate data for processing
                    self.prep_arrays(nm)

                    # set up PET data for processing
                    self.prep_pet(nm)

                    # calculate pet
                    if pet:
                        self.calculate_pet()

                    # calculate runoff and generate monthly potential ET, actual ET, runoff, and soil moisture
                    if runoff:
                        self.calculate_runoff(nm)

                        # update soil moisture (sav) array for next step
                        self.sav_prev = np.copy(self.Sav[:, nm])

                    # channel storage, avg. channel flow (m^3/sec), instantaneous channel flow (m^3/sec)
                    if routing:

                        self.calculate_routing(nm)

                        # update channel storage (chs) array for next step
                        self.chs_prev = np.copy(self.ChStorage[:, nm])

                print("---{0} has finished successfully: {1} seconds ---".format(notify, time.time() - t0))

            elif (runoff_step is None) and (routing_step == 'month'):

                print("---{} in progress... ".format(notify))
                t0 = time.time()

                # set up climate data for processing
                self.prep_arrays()

                # calculate PET
                if pet:
                    print("---Processing PET...")
                    t = time.time()

                    for nm in range(num_steps):

                        # set up PET data for processing
                        self.prep_pet(nm)

                        # calculate pet
                        self.calculate_pet()

                        # archive pet month in array
                        self.pet_out[:, nm] = self.pet_t

                    print("---PET processed in {} seconds---".format(time.time() - t))

                # calculate runoff for all basins all months
                if runoff:
                    print("---Processing Runoff...")
                    t = time.time()

                    self.calculate_runoff()

                    print("---Runoff processed in {} seconds---".format(time.time() - t))

                # update soil moisture (sav) array for use in future mode
                self.sav_prev = np.copy(self.Sav)

                # process routing
                if routing:

                    print("---Processing Routing...")
                    t = time.time()

                    for nm in range(self.s.nmonths):

                        # calculate the number of days in the month
                        #self.mm = np.copy(self.M[nm, 2])

                        # channel storage, avg. channel flow (m^3/sec), instantaneous channel flow (m^3/sec)
                        self.calculate_routing(nm)

                        # update channel storage (chs) arrays for next step
                        self.chs_prev = np.copy(self.ChStorage[:, nm])

                    print("---Routing processed in {} seconds---".format(time.time() - t))

                print("---{0} has finished successfully: {1} seconds ---".format(notify, time.time() - t0))

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
        pass