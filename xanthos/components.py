"""
Components for use in model configurations.

@author   Chris R. Vernon, lixi729
@email:   chris.vernon@pnnl.gov; xinya.li@pnl.gov
@Project: Xanthos 2.0

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute
"""

import numpy as np
import time
import logging

import xanthos.utils.general as helper
import xanthos.calibrate.calibrate_abcd as calib_mod
from xanthos.data_writer.out_writer import OutWriter
from xanthos.diagnostics.diagnostics import Diagnostics
from xanthos.diagnostics.time_series import TimeSeriesPlot
from xanthos.accessible.accessible import AccessibleWater
from xanthos.hydropower.potential import HydropowerPotential
from xanthos.hydropower.actual import HydropowerActual
from xanthos.data_reader.data_load import DataLoader
from xanthos.drought.drought_stats import DroughtStats


class Components:
    """
    Components for use in model configurations.

    @author   Chris R. Vernon, Xinya Li
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
        self.data = DataLoader(config)

        # index arrays
        self.yr_imth_dys = helper.set_month_arrays(self.s.nmonths, self.s.StartYear, self.s.EndYear)

        # pet
        if self.s.pet_module == 'hargreaves':
            self.mth_temp_pet = None
            self.mth_dtr_pet = None
            self.pet_t = None
            self.pet_out = None
            sft = helper.calc_sinusoidal_factor(self.yr_imth_dys)
            self.solar_dec = sft[0]
            self.dr = sft[1]

        elif self.s.pet_module == 'hs':
            pass

        elif self.s.pet_module == 'pm':
            pass

        elif self.s.pet_module == 'thornthwaite':
            pass

        # runoff
        if self.s.runoff_module == 'gwam':
            self.soil_moisture = self.data.soil_moisture
            self.sm_prev = self.data.sm_prev

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

        # aggregated outputs
        self.q = None
        self.ac = None

    def import_core(self):
        """Import desired core modules."""
        global pet_mod
        global runoff_mod
        global routing_mod

        # import desired module for PET
        if self.s.pet_module == 'hargreaves':
            import xanthos.pet.hargreaves as pet_mod

        elif self.s.pet_module == 'hs':
            import xanthos.pet.hargreaves_samani as pet_mod

        elif self.s.pet_module == 'pm':
            import xanthos.pet.penman_monteith as pet_mod

        elif self.s.pet_module == 'thornthwaite':
            import xanthos.pet.thornthwaite as pet_mod

        # import desired module for Runoff
        if self.s.runoff_module == 'gwam':
            import xanthos.runoff.gwam as runoff_mod

        elif self.s.runoff_module == 'abcd':
            import xanthos.runoff.abcd as runoff_mod

        # import desired module for Routing
        if self.s.routing_module == 'mrtm':
            import xanthos.routing.mrtm as routing_mod

    def prep_arrays(self, nm=None):
        """
        Prepare arrays.

        @:param nm:     time-step number
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
            self.mth_temp_pet = np.nan_to_num(self.T)
            self.mth_dtr_pet = np.nan_to_num(self.D)

        elif self.s.pet_module == 'hs':
            pass

        elif self.s.pet_module == 'pm':
            pass

        elif self.s.pet_module == 'thornthwaite':
            pass

    def calculate_pet(self):
        """Calculate monthly potential evapo-transpiration."""
        if self.s.pet_module == 'hargreaves':

            return pet_mod.calculate_pet(self.mth_temp_pet, self.mth_dtr_pet, self.data.lat_radians,
                                         self.mth_solar_dec, self.mth_dr, self.mth_days)
        elif self.s.pet_module == 'hs':

            return pet_mod.execute(self.s, self.data)

        elif self.s.pet_module == 'pm':

            return pet_mod.run_pmpet(self.data, self.s.ncell, self.s.pm_nlcs, self.s.StartYear, self.s.EndYear,
                                     self.s.pm_water_idx, self.s.pm_snow_idx, self.s.pm_lc_years)

        elif self.s.pet_module == 'thornthwaite':

            return pet_mod.execute(self.data.tair, self.data.lat_radians,
                                   self.s.StartYear, self.s.EndYear)

        elif self.s.pet_module == 'none':
            return self.data.pet_out

    def calculate_runoff(self, step_num=None, pet=None):
        """
        Calculate runoff.

        GWAM iterates monthly and takes a feedback of soil moisture content
        (sm_prev).  This is updated for each iteration of a month.

        ABCD iterates internally.

        :param step_num:        Integer for iteration time step (required for GWAM)
        :param pet:             Array of PET (required for ABCD)

        :returns:               PET : Potential Evapotranspiration (mm/month)
                                AET : Actual Evapotranspiration (mm/month)
                                Q   : Runoff (mm/month)
                                Sav : Soil Moisture content (mm/month)
        """
        if self.s.runoff_module == 'gwam':
            rg = runoff_mod.runoffgen(self.pet_t, self.P, self.s, self.soil_moisture, self.sm_prev)

            self.PET[:, step_num], self.AET[:, step_num], self.Q[:, step_num], self.Sav[:, step_num] = rg

        elif self.s.runoff_module == 'abcd':

            rg = runoff_mod.abcd_execute(n_basins=self.s.n_basins, basin_ids=self.data.basin_ids,
                                         pet=pet, precip=self.data.precip, tmin=self.data.tmin,
                                         calib_file=self.s.calib_file, n_months=self.s.nmonths,
                                         spinup_steps=self.s.runoff_spinup, jobs=self.s.ro_jobs)

            self.PET, self.AET, self.Q, self.Sav = rg

        else:

            # if user is providing a custom runoff file
            if self.s.alt_runoff is not None:
                self.Q = np.load(self.s.alt_runoff)

    def calculate_routing(self, runoff):
        """
        Calculate routing.

        Routing takes a simulated runoff (Q) from the runoff output and
        previous channel storage (chs_prev) from previous channel storage.

        :returns:                   ChStorage     : Channel storage (m3)
                                    Avg_ChFlow    : Average streamflow (m3/s)
                                    instream_flow : Streamflow (m3/s)
        """
        if self.s.routing_module == 'mrtm':
            # initialize routing data
            self.flow_dist = self.data.flow_dist
            self.flow_dir = self.data.flow_dir
            self.instream_flow = self.data.instream_flow
            self.str_velocity = self.data.str_velocity
            self.chs_prev = self.data.chs_prev

            self.dsid = routing_mod.downstream(self.data.coords, self.flow_dir, self.s)
            self.upid = routing_mod.upstream(self.data.coords, self.dsid, self.s)
            self.um = routing_mod.upstream_genmatrix(self.upid)

            # process spin up for channel storage from historic period
            for nm in range(0, self.s.routing_spinup, 1):

                sr = routing_mod.streamrouting(self.flow_dist, self.chs_prev, self.instream_flow, self.str_velocity,
                                               runoff[:, nm], self.data.area, self.yr_imth_dys[nm, 2],
                                               self.routing_timestep_hours, self.um)

                self.ChStorage[:, nm], self.Avg_ChFlow[:, nm], self.instream_flow = sr

                # update channel storage (chs) arrays for next step
                self.chs_prev = np.copy(self.ChStorage[:, nm])

            # run routing simulation
            for nm in range(self.s.nmonths):
                # channel storage, avg. channel flow (m^3/sec), instantaneous channel flow (m^3/sec)
                sr = routing_mod.streamrouting(self.flow_dist, self.chs_prev, self.instream_flow, self.str_velocity,
                                               runoff[:, nm], self.data.area, self.yr_imth_dys[nm, 2],
                                               self.routing_timestep_hours, self.um)

                self.ChStorage[:, nm], self.Avg_ChFlow[:, nm], self.instream_flow = sr

                # update channel storage (chs) arrays for next step
                self.chs_prev = np.copy(self.ChStorage[:, nm])

            return self.Avg_ChFlow

    def simulation(self, run_pet, run_runoff, run_routing, pet_num_steps=0, runoff_num_steps=0, routing_num_steps=0,
                   notify='simulation'):
        """
        Run model simulation for a defined configuration.

        :param run_pet:             Boolean - True to run PET module, False to load PET from file
        :param run_runoff:          Boolean - True to run runoff module
        :param run_routing:         Boolean - True to run routing module
        :param pet_num_steps:       The number of time steps for running the PET module; if 0, the
                                    PET module should iterate internally
        :param runoff_num_steps:    The number of time steps for running the runoff module; if 0,
                                    the runoff module should iterate internally
        :param routing_num_steps:   The number of time steps for running the routing module; if 0,
                                    the routing module should iterate internally
        :param notify:              A string that is used to add to log print that describes
                                    whether the simulation is spin-up or regular
        """
        # default to calibration if selected
        if self.s.calibrate:
            self.calibrate()
            return

        logging.info("---{} in progress...".format(notify))
        t0 = time.time()

        # Run PET
        if run_pet:
            logging.info("\tProcessing PET...")
            t = time.time()

            # Calculate PET step by step
            if pet_num_steps > 0:
                pet_out = np.zeros_like(self.data.precip)

                for nm in range(pet_num_steps):
                    # set up climate data for processing (used by hargreaves)
                    self.prep_arrays(nm)

                    # set up PET data for processing
                    self.prep_pet(nm)

                    # calculate pet
                    pet_out[:, nm] = self.calculate_pet()

            # Calculate PET all at once
            else:
                pet_out = self.calculate_pet()

            logging.info("\tPET processed in {} seconds---".format(time.time() - t))

        # Otherwise calculate_pet() will load user-provided PET dataset
        else:
            pet_out = self.calculate_pet()

        # Process runoff
        if run_runoff:
            logging.info("\tProcessing Runoff...")
            t = time.time()

            # Calculate runoff step by step (GWAM)
            if runoff_num_steps > 0:
                for nm in range(runoff_num_steps):
                    self.pet_t = pet_out[:, nm]

                    # calculate runoff and generate monthly potential ET, actual ET, runoff, and soil moisture
                    self.calculate_runoff(step_num=nm)

                    # update soil moisture (sav) array for next step
                    self.sm_prev = np.copy(self.Sav[:, nm])

            # Calculate runoff all at once (ABCD)
            else:
                self.calculate_runoff(pet=pet_out)

            logging.info("\tRunoff processed in {} seconds---".format(time.time() - t))

        # Process routing (only MRTM currently implemented -- don't allow timestep iteration)
        if run_routing:
            logging.info("\tProcessing Routing...")
            t = time.time()

            # channel storage, avg. channel flow (m^3/sec), instantaneous channel flow (m^3/sec)
            self.calculate_routing(self.Q)

            logging.info("\tRouting processed in {} seconds---".format(time.time() - t))

        logging.info("---{0} has finished successfully: {1} seconds ---".format(notify, time.time() - t0))

    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    # OPTIONAL POST-PROCESSING MODULE METHODS
    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    def drought(self):
        """Run drought module."""
        if self.s.CalculateDroughtStats:
            logging.info("---Start Drought Statistics:")
            t0 = time.time()

            DroughtStats(self.s, self.Q, self.Sav)

            logging.info("---Drought Statistics has finished successfully: %s seconds ------" % (time.time() - t0))

    def accessible_water(self):
        """Run accessible water module."""
        if self.s.CalculateAccessibleWater:
            logging.info("---Start Accessible Water:")
            t0 = time.time()

            AccessibleWater(self.s, self.data, self.Q)

            logging.info("---Accessible Water has finished successfully: %s seconds ------" % (time.time() - t0))

    def hydropower_potential(self):
        """Run hydropower potential module."""
        if self.s.CalculateHydropowerPotential:
            logging.info("---Start Hydropower Potential:")
            t0 = time.time()

            HydropowerPotential(self.s, self.Avg_ChFlow)

            logging.info("---Hydropower Potential has finished successfully: %s seconds ------" % (time.time() - t0))

    def hydropower_actual(self):
        """Run hydropower actual module."""
        if self.s.CalculateHydropowerActual:
            logging.info("---Start Hydropower Actual:")
            t0 = time.time()

            HydropowerActual(self.s, self.Avg_ChFlow)

            logging.info("---Hydropower Actual has finished successfully: %s seconds ------" % (time.time() - t0))

    def diagnostics(self):
        """Run diagnostics."""
        if self.s.PerformDiagnostics:
            logging.info("---Start Diagnostics:")
            t0 = time.time()

            Diagnostics(self.s, self.Q, self.data)

            logging.info("---Diagnostics has finished successfully: %s seconds ------" % (time.time() - t0))

    def output_simulation(self):
        """
        Output simulation results.

        This step converts the data to a pandas DataFrame in the user-specified format.
        """
        logging.info("---Output simulation results:")
        t0 = time.time()

        all_outputs = {
            'pet': self.PET,
            'aet': self.AET,
            'q': self.Q,
            'soilmoisture': self.Sav,
            'avgchflow': self.Avg_ChFlow
        }

        output_writer = OutWriter(self.s, self.data.area, all_outputs)
        output_writer.write()

        try:
            self.q = output_writer.get('q')
        except ValueError:
            self.q = self.Q

        try:
            self.ac = output_writer.get('avgchflow')
        except ValueError:
            self.ac = self.Avg_ChFlow

        output_writer.write_aggregates(self.data, self.q, self.s.AggregateRunoffBasin, self.s.AggregateRunoffCountry,
                                       self.s.AggregateRunoffGCAMRegion)

        logging.info("---Output finished: %s seconds ---" % (time.time() - t0))

    def plots(self):
        """Create time series plots."""
        if self.s.CreateTimeSeriesPlot:
            logging.info("---Creating Time Series Plots:")
            t0 = time.time()

            TimeSeriesPlot(self.s, self.q, self.ac, self.data)

            logging.info("---Plots has finished successfully: %s seconds ------" % (time.time() - t0))

    def calibrate(self):
        """Run calibration to generate parameters for the ABCD model."""
        logging.info("---Processing PET...")
        t = time.time()

        pet_out = self.calculate_pet()

        logging.info("---PET processed in {} seconds---".format(time.time() - t))

        logging.info("---Running calibration:")

        calib_mod.calibrate_all(settings=self.s, data=self.data, pet=pet_out, router_function=self.calculate_routing)
