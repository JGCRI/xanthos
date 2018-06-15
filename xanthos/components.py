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

import data_reader.data_load as fetch
import utils.general as helper
import utils.math as umth
import calibrate.calibrate_abcd as calib_mod
from data_writer.out_writer import OUTWriter
from diagnostics.aggregation import Aggregation
from diagnostics.diagnostics import Diagnostics
from diagnostics.time_series import TimeSeriesPlot
from accessible.accessible import AccessibleWater
from hydropower.potential import HydropowerPotential
from hydropower.actual import HydropowerActual
from data_reader.data_load import LoadData


class Components:
    """
    Components for use in model configurations.

    @author   Chris R. Vernon , Xinya Li
    @email:   chris.vernon@pnnl.gov; xinya.li@pnl.gov
    @Project: Xanthos 2.0

    License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

    Copyright (c) 2017, Battelle Memorial Institute
    """

    def __init__(self, config):

        self.s = config

        # import desired modules for PET, Runoff, and Routing
        self.import_core()

        # load data
        self.data = LoadData(config)

        # index arrays
        self.yr_imth_dys = helper.set_month_arrays(self.s.nmonths, self.s.StartYear, self.s.EndYear)
        self.map_index = umth.sub2ind([self.s.ngridrow, self.s.ngridcol],
                                      self.data.coords[:, 4].astype(int) - 1,
                                      self.data.coords[:, 3].astype(int) - 1)

        # pet
        if self.s.pet_module == 'hargreaves':
            self.mth_temp_pet = None
            self.mth_dtr_pet = None
            self.pet_t = None
            self.pet_out = None
            sft = helper.calc_sinusoidal_factor(self.yr_imth_dys)
            self.solar_dec = sft[0]
            self.dr = sft[1]

        elif self.s.pet_module == 'pm':
            pass

        # runoff
        if self.s.runoff_module == 'gwam':
            self.soil_moisture = None
            self.grid_area = None
            self.sm_prev = None

        elif self.s.runoff_module == 'abcd':
            pass

        # routing
        if self.s.routing_module == 'mrtm':
            self.flow_dist = None
            self.flow_dir = None
            self.instream_flow = None
            self.str_velocity = None
            self.dsid = None
            self.upid = None
            self.um = None
            self.routing_timestep_hours = 3 * 3600
            self.chs_prev = None

        # outputs
        self.PET = np.zeros(shape=(self.s.ncell, self.s.nmonths))
        self.AET = np.zeros(shape=(self.s.ncell, self.s.nmonths))
        self.Q = np.zeros(shape=(self.s.ncell, self.s.nmonths))
        self.Sav = np.zeros(shape=(self.s.ncell, self.s.nmonths))
        self.ChStorage = np.zeros(shape=(self.s.ncell, self.s.nmonths))
        self.Avg_ChFlow = np.zeros(shape=(self.s.ncell, self.s.nmonths))

        # simulation
        self.P = None
        self.T = None
        self.D = None
        self.mth_solar_dec = None
        self.mth_dr = None
        self.mth_days = None

        # outputs
        self.q = None
        self.ac = None

    def import_core(self):
        """
        Import desired core modules.
        """
        global pet_mod
        global runoff_mod
        global routing_mod

        # import desired module for PET
        if self.s.pet_module == 'hargreaves':
            import pet.hargreaves as pet_mod

        elif self.s.pet_module == 'pm':
            import pet.penman_monteith as pet_mod

        # import desired module for Runoff
        if self.s.runoff_module == 'gwam':
            import runoff.gwam as runoff_mod

        elif self.s.runoff_module == 'abcd':
            import runoff.abcd as runoff_mod

        # import desired module for Routing
        if self.s.routing_module == 'mrtm':
            import routing.mrtm as routing_mod

    def prep_arrays(self, nm=None):
        """
        Prepare arrays.
        """
        if nm is None:
            self.P = np.copy(self.data.precip)  # keep nan in P
            self.T = np.nan_to_num(self.data.temp)  # nan to zero
            self.D = np.nan_to_num(self.data.dtr)  # nan to zero

        else:
            self.P = np.copy(self.data.precip[:, nm])  # keep nan in P
            self.T = np.nan_to_num(self.data.temp[:, nm])  # nan to zero
            self.D = np.nan_to_num(self.data.dtr[:, nm])  # nan to zero

    def prep_pet(self, nm):
        """
        Prepare PET for processing.

        @:param nm:                     time-step number

        @:return mth_solar_dec:         solar declination for the target step
        @:return mth_dr:                inverse relative distance Earth-Sun for the target step
        @:return mth_days:              number of days in the target month
        @:return temp_pet:              PET for the target month
        @:return dtr_pet:               Daily temperature range for the target month
        """
        if self.s.pet_module == 'hargreaves':

            self.mth_solar_dec = np.copy(self.solar_dec[nm])
            self.mth_dr = np.copy(self.dr[nm])
            self.mth_days = np.copy(self.yr_imth_dys[nm, 2])
            self.mth_temp_pet = np.nan_to_num(self.data.temp[:, nm])
            self.mth_dtr_pet = np.nan_to_num(self.data.dtr[:, nm])

        elif self.s.pet_module == 'pm':
            pass

    def calculate_pet(self, nm=None):
        """
        Calculate monthly potential evapo-transpiration.
        """
        if self.s.pet_module == 'hargreaves':

            return pet_mod.calculate_pet(self.mth_temp_pet, self.mth_dtr_pet, self.data.lat_radians, self.mth_solar_dec,
                                         self.mth_dr, self.mth_days)

        elif self.s.pet_module == 'pm':

            return pet_mod.run_pmpet(self.data, self.s.ncell, self.s.pm_nlcs, self.s.StartYear, self.s.EndYear,
                                     self.s.pm_water_idx, self.s.pm_snow_idx, self.s.pm_lc_years)

        elif self.s.pet_module == 'none':
            return self.data.pet_out

    def calculate_runoff(self, nm=None, pet=None):
        """
        Calculate runoff.

        Runoff takes a feedback of soil moisture content (sm_prev).  This is updated for each iteration of a month.

        :returns:               PET : Potential Evapotranspiration (mm/month)
                                AET : Actual Evapotranspiration (mm/month)
                                Q   : Runoff (mm/month)
                                Sav : Soil Moisture content (mm/month)
        """
        if self.s.runoff_module == 'gwam':

            rg = runoff_mod.runoffgen(self.pet_t, self.P, self.s, self.soil_moisture, self.sm_prev)

            self.PET[:, nm], self.AET[:, nm], self.Q[:, nm], self.Sav[:, nm] = rg

        elif self.s.runoff_module == 'abcd':

            rg = runoff_mod.abcd_execute(n_basins=self.s.n_basins, basin_ids=self.data.basin_ids,
                                         pet=pet, precip=self.data.precip, tmin=np.nan_to_num(self.data.tmin),
                                         calib_file=self.s.calib_file, n_months=self.s.nmonths,
                                         spinup_steps=self.s.runoff_spinup, jobs=self.s.ro_jobs)

            self.PET, self.AET, self.Q, self.Sav = rg

        else:

            # if user is providing a custom runoff file
            if self.s.alt_runoff is not None:
                self.Q = np.load(self.s.alt_runoff)

    def calculate_routing(self, runoff):
        """
        Calculate routing.  Routing takes a simulated runoff (Q) from the runoff output and
        previous channel storage (chs_prev) from previous channel storage.

        :returns:                   ChStorage     : Channel storage (m3)
                                    Avg_ChFlow    : Average streamflow (m3/s)
                                    instream_flow : Streamflow (m3/s)
        """
        if self.s.routing_module == 'mrtm':

            # load routing data
            self.flow_dist = fetch.load_routing_data(self.s.FlowDis, self.s.ngridrow, self.s.ngridcol,
                                                     self.map_index, rep_val=1000)
            self.flow_dir = fetch.load_routing_data(self.s.FlowDir, self.s.ngridrow, self.s.ngridcol, self.map_index)
            self.instream_flow = np.zeros((self.s.ncell,), dtype=float)
            self.str_velocity = fetch.load_routing_data(self.s.strm_veloc, self.s.ngridrow, self.s.ngridcol,
                                                        self.map_index, rep_val=0)
            self.dsid = routing_mod.downstream(self.data.coords, self.flow_dir, self.s)
            self.upid = routing_mod.upstream(self.data.coords, self.dsid, self.s)
            self.um = routing_mod.upstream_genmatrix(self.upid)
            self.chs_prev = fetch.load_chs_data(self.s)

            # process spin up for channel storage from historic period
            for nm in range(0, self.s.routing_spinup, 1):

                sr = routing_mod.streamrouting(self.flow_dist, self.chs_prev, self.instream_flow, self.str_velocity,
                                               runoff[:, nm], self.data.area, self.yr_imth_dys[nm, 2],
                                               self.routing_timestep_hours, self.um)

                self.ChStorage[:, nm], self.Avg_ChFlow[:, nm], self.instream_flow = sr

                # update channel storage (chs) arrays for next step
                self.chs_prev = np.copy(self.ChStorage[:, nm])

            # run routing simumlation
            for nm in range(self.s.nmonths):
                # channel storage, avg. channel flow (m^3/sec), instantaneous channel flow (m^3/sec)
                sr = routing_mod.streamrouting(self.flow_dist, self.chs_prev, self.instream_flow, self.str_velocity,
                                               runoff[:, nm], self.data.area, self.yr_imth_dys[nm, 2],
                                               self.routing_timestep_hours, self.um)

                self.ChStorage[:, nm], self.Avg_ChFlow[:, nm], self.instream_flow = sr

                # update channel storage (chs) arrays for next step
                self.chs_prev = np.copy(self.ChStorage[:, nm])

            return self.Avg_ChFlow

    def simulation(self, pet=True, pet_num_steps=0, pet_step='month',
                   runoff=True, runoff_num_steps=0, runoff_step='month',
                   routing=True, routing_num_steps=0, routing_step='month',
                   notify='simulation'):
        """
        Run model simulation for a defined configuration.

        :param num_steps:           The number of time steps to process (INT)
        :param pet:                 True if running PET, False if embedded in runoff model
        :param runoff:              True if running Runoff, False if not
        :param runoff_step:         The time unit as a string; if None the runoff model iterates internally, else 'month'
        :param routing_num_steps:   The number of steps for to run the routing module; different on spin-up, same as runoff
                                    when running normal.  Specified in months in the config file.
        :param routing:             True if running Routing, False if not
        :param routing_step:        The time unit as a string; if None the routing model iterates internally, else 'month'
        :param notify:              A string that is used to add to log print that describes whether the simulation is
                                    spin-up or regular
        """

        # default to calibration if selected
        if self.s.calibrate == 1:
            self.calibrate()

        else:
            # pass simulation if there are no steps to process
            if (pet_num_steps + runoff_num_steps + routing_num_steps) == 0:
                pass

            else:

                # --------------------------------------------------
                # USED FOR THE FOLLOWING CONFIGURATIONS:
                #
                # hargreaves-gwam-mrtm
                # --------------------------------------------------
                if (pet_step == 'month') and (runoff_step == 'month') and (routing_step == 'month'):

                    print("---{} in progress...".format(notify))
                    t0 = time.time()

                    if pet:

                        print("\tProcessing PET...")
                        t = time.time()

                        for nm in range(pet_num_steps):
                            # set up climate data for processing
                            self.prep_arrays(nm)

                            # set up PET data for processing
                            self.prep_pet(nm)

                            # calculate pet
                            self.calculate_pet()

                        print("\tPET processed in {} seconds---".format(time.time() - t))

                    # for the case where the user provides a PET dataset
                    else:
                        # load user provided data
                        self.calculate_pet()

                    if runoff:

                        print("\tProcessing Runoff...")
                        t = time.time()

                        for nm in range(runoff_num_steps):

                            # calculate runoff and generate monthly potential ET, actual ET, runoff, and soil moisture
                            if runoff:
                                self.calculate_runoff(nm)

                                # update soil moisture (sav) array for next step
                                self.sm_prev = np.copy(self.Sav[:, nm])

                        print("\tRunoff processed in {} seconds---".format(time.time() - t))

                    # channel storage, avg. channel flow (m^3/sec), instantaneous channel flow (m^3/sec)
                    if routing:

                        print("\tProcessing Routing...")
                        t = time.time()

                        # channel storage, avg. channel flow (m^3/sec), instantaneous channel flow (m^3/sec)
                        self.calculate_routing(self.Q)

                        print("\tRouting processed in {} seconds---".format(time.time() - t))

                    print("---{0} has finished successfully: {1} seconds ---".format(notify, time.time() - t0))

                # --------------------------------------------------
                # USED FOR THE FOLLOWING CONFIGURATIONS:
                #
                # hargreaves-abcd-mrtm
                # --------------------------------------------------
                elif (pet_step == 'month') and (runoff_step is None) and (routing_step == 'month'):

                    print("---{} in progress... ".format(notify))
                    t0 = time.time()

                    # calculate PET
                    if pet:
                        print("\tProcessing PET...")
                        t = time.time()
                        pet_out = np.zeros_like(self.data.precip)

                        for nm in range(pet_num_steps):
                            # set up PET data for processing
                            self.prep_pet(nm)

                            # calculate pet
                            pet_out[:, nm] = self.calculate_pet()

                        print("\tPET processed in {} seconds---".format(time.time() - t))

                    # for the case where the user provides a PET dataset
                    else:
                        # load user provided data
                        pet_out = self.calculate_pet()

                    # calculate runoff for all basins all months
                    if runoff:
                        print("\tProcessing Runoff...")
                        t = time.time()

                        self.calculate_runoff(pet=pet_out)

                        print("\tRunoff processed in {} seconds---".format(time.time() - t))

                    # process routing
                    if routing:

                        print("\tProcessing Routing...")
                        t = time.time()

                        # channel storage, avg. channel flow (m^3/sec), instantaneous channel flow (m^3/sec)
                        self.calculate_routing(self.Q)

                        print("\tRouting processed in {} seconds---".format(time.time() - t))

                    print("---{0} has finished successfully: {1} seconds ---".format(notify, time.time() - t0))

                # --------------------------------------------------
                # USED FOR THE FOLLOWING CONFIGURATIONS:
                #
                # pm-abcd-mrtm
                # --------------------------------------------------
                elif (pet_step is None) and (runoff_step is None) and (routing_step == 'month'):

                    print("---{} in progress... ".format(notify))
                    t0 = time.time()

                    # calculate PET
                    if pet:
                        print("\tProcessing PET...")
                        t = time.time()

                        # calculate pet
                        pet_out = self.calculate_pet()

                        print("\tPET processed in {} seconds---".format(time.time() - t))

                    # for the case where the user provides a PET dataset
                    else:
                        # load user provided data
                        pet_out = self.calculate_pet()

                    # calculate runoff for all basins all months
                    if runoff:
                        print("\tProcessing Runoff...")
                        t = time.time()

                        self.calculate_runoff(pet=pet_out)

                        print("\tRunoff processed in {} seconds---".format(time.time() - t))

                    # process routing
                    if routing:

                        print("\tProcessing Routing...")
                        t = time.time()

                        # process spin up for channel storage from historic period
                        self.calculate_routing(self.Q)

                        print("\tRouting processed in {} seconds---".format(time.time() - t))

                    print("---{0} has finished successfully: {1} seconds ---".format(notify, time.time() - t0))

    def accessible_water(self):
        """
        Run accessible water module
        """
        if self.s.CalculateAccessibleWater:
            print("---Start Accessible Water:")
            t0 = time.time()

            AccessibleWater(self.s, self.data, self.Q)

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

            Diagnostics(self.s, self.Q, self.data)

            print("---Diagnostics has finished successfully: %s seconds ------" % (time.time() - t0))

    def output_simulation(self):
        """
        Output simulation results.  This step both converts the data to the user specified format and
        """
        print("---Output simulation results:")
        t0 = time.time()

        self.q, self.ac = OUTWriter(self.s, self.data.area, self.PET, self.AET, self.Q, self.Sav, self.ChStorage,
                                    self.Avg_ChFlow)

        print("---Output finished: %s seconds ---" % (time.time() - t0))

    def aggregate_outputs(self):
        """
        Aggregation by Basin, Country, and/or Region.
        """
        if self.s.AggregateRunoffBasin > 0 or self.s.AggregateRunoffCountry > 0 or self.s.AggregateRunoffGCAMRegion > 0:
            print("---Start Aggregation:")
            t0 = time.time()

            Aggregation(self.s, self.data, self.q)

            print("---Aggregation has finished successfully: %s seconds ------" % (time.time() - t0))

    def plots(self):
        """
        Create time series plots.
        """
        if self.s.CreateTimeSeriesPlot:
            print("---Creating Time Series Plots:")
            t0 = time.time()

            TimeSeriesPlot(self.s, self.q, self.ac, self.data)

            print("---Plots has finished successfully: %s seconds ------" % (time.time() - t0))

    def calibrate(self):
        """
        Run calibration to generate parameters for the ABCD model
        """
        print("---Processing PET...")
        t = time.time()

        pet_out = self.calculate_pet()

        print("---PET processed in {} seconds---".format(time.time() - t))

        print("---Running calibration:")

        calib_mod.calibrate_all(settings=self.s, data=self.data, pet=pet_out, router_function=self.calculate_routing)
