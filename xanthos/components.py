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
import xanthos.calibrate.calibrate_managed as calib_managed
from xanthos.data_writer.out_writer import OutWriter
from xanthos.diagnostics.diagnostics import Diagnostics
from xanthos.diagnostics.time_series import TimeSeriesPlot
from xanthos.accessible.accessible import AccessibleWater
from xanthos.hydropower.potential import HydropowerPotential
from xanthos.hydropower.actual import HydropowerActual
from xanthos.drought.drought_stats import DroughtStats
from xanthos.data_reader.data_penman_monteith import DataPenmanMonteith
from xanthos.data_reader.data_hargreaves import DataHargreaves
from xanthos.data_reader.data_hargreaves_semani import DataHargreavesSemani
from xanthos.data_reader.data_thornthwaite import DataThornthwaite
from xanthos.data_reader.data_reference import DataReference
from xanthos.data_reader.data_gwam import DataGwam
from xanthos.data_reader.data_mrtm import DataMrtm
from xanthos.data_reader.data_mrtm_managed import DataMrtmManaged
from xanthos.data_reader.data_mrtm import DataUtils
from xanthos.data_reader.data_abcd import DataAbcd
from xanthos.data_reader.data_diagnostics import DataDiagnostics
from xanthos.data_reader.data_calibration import DataCalibration
from xanthos.data_reader.data_calibration_managed import DataCalibrationManaged


class Components(DataUtils):
    """
    Components for use in model configurations.

    @author   Chris R. Vernon, Xinya Li
    @email:   chris.vernon@pnnl.gov; xinya.li@pnl.gov
    @Project: Xanthos 2.0

    License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

    Copyright (c) 2017, Battelle Memorial Institute
    """

    def __init__(self, config):

        super().__init__(nmonths=config.nmonths)

        self.s = config

        # import desired modules for PET, Runoff, and Routing
        self.import_core()

        # load reference data
        self.data_reference = DataReference(config)

        # index arrays
        self.yr_imth_dys = helper.set_month_arrays(self.s.nmonths, self.s.StartYear, self.s.EndYear)

        # setup PET data
        if self.s.pet_module == 'hargreaves':
            self.mth_temp_pet = None
            self.mth_dtr_pet = None
            self.pet_t = None
            self.pet_out = None
            sft = helper.calc_sinusoidal_factor(self.yr_imth_dys)
            self.solar_dec = sft[0]
            self.dr = sft[1]

            self.data_hargreaves = DataHargreaves(config)

        elif self.s.pet_module == 'hs':
            self.data_hargreaves_semani = DataHargreavesSemani(config)

        elif self.s.pet_module == 'pm':
            self.penman_monteith_data = DataPenmanMonteith(config)

        elif self.s.pet_module == 'thornthwaite':
            self.data_thornthwaite = DataThornthwaite(config)

        # runoff
        if self.s.runoff_module == 'gwam':
            self.data_gwam = DataGwam(config, self.data_reference.area, self.data_reference.region_ids,
                                      self.data_reference.country_ids, self.data_reference.basin_ids)

        elif self.s.runoff_module in ('abcd', 'abcd_managed'):
            self.data_abcd = DataAbcd(config)

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

            self.data_mrtm = DataMrtm(config_obj=config)

        elif self.s.routing_module == 'mrtm_managed':
            self.flow_dist = None
            self.flow_dir = None
            self.instream_flow = None
            self.str_velocity = None
            self.dsid = None
            self.upid = None
            self.um = None
            self.routing_timestep_hours = 3 * 3600
            self.chs_prev = None

            self.data_mrtm_managed = DataMrtmManaged(config_obj=config)
            self.calib_data = DataCalibrationManaged(config_obj=self.s)

        # outputs
        self.PET = np.zeros(shape=(self.s.ncell, self.s.nmonths))
        self.AET = np.zeros(shape=(self.s.ncell, self.s.nmonths))
        self.Q = np.zeros(shape=(self.s.ncell, self.s.nmonths))
        self.Sav = np.zeros(shape=(self.s.ncell, self.s.nmonths))
        self.ChStorage = np.zeros(shape=(self.s.ncell, self.s.nmonths))
        self.Avg_ChFlow = np.zeros(shape=(self.s.ncell, self.s.nmonths))

        # simulation
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

        elif self.s.runoff_module == 'abcd_managed':
            import xanthos.runoff.abcd_managed as runoff_mod

        # import desired module for Routing
        if self.s.routing_module == 'mrtm':
            import xanthos.routing.mrtm as routing_mod

        elif self.s.routing_module == 'mrtm_managed':
            import xanthos.routing.mrtm_managed as routing_mod

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
            self.mth_temp_pet = np.nan_to_num(self.data_hargreaves.temp)
            self.mth_dtr_pet = np.nan_to_num(self.data_hargreaves.dtr)

        elif self.s.pet_module == 'hs':
            pass

        elif self.s.pet_module == 'pm':
            pass

        elif self.s.pet_module == 'thornthwaite':
            pass

    def calculate_pet(self):
        """Calculate monthly potential evapo-transpiration."""
        if self.s.pet_module == 'hargreaves':

            return pet_mod.calculate_pet(self.mth_temp_pet, self.mth_dtr_pet, self.data_reference.lat_radians,
                                         self.mth_solar_dec, self.mth_dr, self.mth_days)
        elif self.s.pet_module == 'hs':

            return pet_mod.execute(self.s, self.data_hargreaves_semani)

        elif self.s.pet_module == 'pm':

            return pet_mod.run_pmpet(self.penman_monteith_data, self.s.ncell, self.s.pm_nlcs, self.s.StartYear, self.s.EndYear,
                                     self.s.pm_water_idx, self.s.pm_snow_idx, self.s.pm_lc_years)

        elif self.s.pet_module == 'thornthwaite':

            return pet_mod.execute(self.data_thornthwaite.tair, self.data_reference.lat_radians,
                                   self.s.StartYear, self.s.EndYear)

        elif self.s.pet_module == 'none':
            return self.load_to_array(self.s.pet_file)

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
            rg = runoff_mod.runoffgen(self.pet_t, self.data_gwam.precip, self.s, self.data_gwam.soil_moisture,
                                      self.data_gwam.sm_prev)

            self.PET[:, step_num], self.AET[:, step_num], self.Q[:, step_num], self.Sav[:, step_num] = rg

        elif self.s.runoff_module == 'abcd':

            rg = runoff_mod.abcd_execute(n_basins=self.s.n_basins, basin_ids=self.data_reference.basin_ids,
                                         pet=pet, precip=self.data_abcd.precip, tmin=self.data_abcd.tmin,
                                         calib_file=self.s.calib_file, n_months=self.s.nmonths,
                                         spinup_steps=self.s.runoff_spinup, jobs=self.s.ro_jobs)

            self.PET, self.AET, self.Q, self.Sav = rg

        elif self.s.runoff_module == 'abcd_managed':

            self.params = [self.s.a_param, self.s.b_param, self.s.c_param, self.s.d_param, self.s.m_param]

            ncell_in_basin = len(self.calib_data.basin_ids)

            self.target_basin_ids = np.zeros(ncell_in_basin)

            he = runoff_mod.AbcdManaged(pars=self.params,
                                        soil_water_initial=np.squeeze(self.calib_data.sm),
                                        pet=pet,
                                        precip=self.data_abcd.precip,
                                        tmin=self.data_abcd.tmin,
                                        basin_ids=self.target_basin_ids,
                                        process_steps=self.nmonths,
                                        spinup_steps=self.s.runoff_spinup,
                                        method="dist")

            he.emulate()

            self.PET = he.pet.T
            self.AET = he.actual_et.T
            self.Q = he.rsim.T
            self.Sav = he.soil_water_storage.T

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
            self.flow_dist = self.data_mrtm.flow_dist
            self.flow_dir = self.data_mrtm.flow_dir
            self.instream_flow = self.data_mrtm.instream_flow
            self.str_velocity = self.data_mrtm.str_velocity
            self.chs_prev = self.data_mrtm.chs_prev

            self.dsid = routing_mod.downstream(self.data_reference.coords, self.flow_dir, self.s)
            self.upid = routing_mod.upstream(self.data_reference.coords, self.dsid, self.s)
            self.um = routing_mod.upstream_genmatrix(self.upid)

            # process spin up for channel storage from historic period
            for nm in range(0, self.s.routing_spinup, 1):

                sr = routing_mod.streamrouting(self.flow_dist, self.chs_prev, self.instream_flow, self.str_velocity,
                                               runoff[:, nm], self.data_reference.area, self.yr_imth_dys[nm, 2],
                                               self.routing_timestep_hours, self.um)

                self.ChStorage[:, nm], self.Avg_ChFlow[:, nm], self.instream_flow = sr

                # update channel storage (chs) arrays for next step
                self.chs_prev = np.copy(self.ChStorage[:, nm])

            # run routing simulation
            for nm in range(self.s.nmonths):

                # channel storage, avg. channel flow (m^3/sec), instantaneous channel flow (m^3/sec)
                sr = routing_mod.streamrouting(self.flow_dist, self.chs_prev, self.instream_flow, self.str_velocity,
                                               runoff[:, nm], self.data_reference.area, self.yr_imth_dys[nm, 2],
                                               self.routing_timestep_hours, self.um)

                self.ChStorage[:, nm], self.Avg_ChFlow[:, nm], self.instream_flow = sr

                # update channel storage (chs) arrays for next step
                self.chs_prev = np.copy(self.ChStorage[:, nm])

            return self.Avg_ChFlow

        elif self.s.routing_module == 'mrtm_managed':

            # index for gauge station locations
            self.grdcData_info = np.copy(self.calib_data.grdc_coord_index_file)
            self.grdc_xanthosID = self.grdcData_info[np.where(self.grdcData_info[:, 0] == self.s.basin_num), 1][0][0] - 1

            # load additional data
            self.wdirr = np.copy(self.calib_data.total_demand_cumecs)
            self.irrmean = np.mean(self.wdirr, axis=1)  # mean demand
            self.ppose = np.copy(self.calib_data.purpose)
            self.cpa = np.copy(self.calib_data.capacity) * 10 ** 6  # m3
            self.HP_Release = np.copy(self.calib_data.hp_release)
            self.WConsumption = np.copy(self.calib_data.water_consumption)
            self.chs_ini = np.copy(self.calib_data.ini_channel_storage)
            self.Initial_instream_flow = np.copy(self.calib_data.instream_flow_natural)
            self.SM = np.squeeze(np.copy(self.calib_data.sm))
            self.mtifl_natural = np.copy(self.calib_data.mtif_natural)
            self.maxtifl_natural = np.copy(self.calib_data.maxtif_natural)

            self.flow_dist = self.data_mrtm_managed.flow_dist
            grid_size = np.sqrt(self.data_reference.area) * 1000
            nn_grids = np.where(self.flow_dist < grid_size)[0]
            self.flow_dist[nn_grids] = grid_size[nn_grids]

            self.flow_dir = self.data_mrtm_managed.flow_dir

            self.str_v = self.data_mrtm_managed.str_velocity
            self.str_v[self.str_v < 0.01] = 0.01
            self.str_velocity = self.str_v * self.s.beta_param

            self.chs_prev = self.data_mrtm_managed.chs_prev

            self.dsid = routing_mod.downstream(self.data_reference.coords, self.flow_dir, self.data_mrtm_managed.NGRIDROW,
                                               self.data_mrtm_managed.NGRIDCOL)
            self.upid = routing_mod.upstream(self.data_reference.coords, self.dsid, self.data_mrtm_managed.NGRIDROW,
                                             self.data_mrtm_managed.NGRIDCOL)
            self.um, self.up = routing_mod.upstream_genmatrix(self.upid)

            # routing initializations
            self.chs_prev = np.squeeze(self.chs_ini)  # initial channel storage in m3
            self.instream_flow = np.squeeze(self.Initial_instream_flow)  # initial channel flow in m3 /s
            self.Sini = 0.5 * np.squeeze(self.cpa)  # initial reservoir storage in m3
            self.res_prev = np.squeeze(self.cpa)  # reservoir storage from t-1 in m3

            # Preallocation for variables to be returned from routing module
            self.Avg_ChFlow = np.zeros_like(self.data_abcd.precip)  # average channel flow m3/s
            self.ChStorage = np.zeros_like(self.data_abcd.precip)  # channel storage
            self.ResStorage = np.zeros_like(self.data_abcd.precip)  # reservoir storage
            self.Qout_res_avg = np.zeros_like(self.data_abcd.precip)  # reservoir storage
            self.Qin_res_avg = np.zeros_like(self.data_abcd.precip)  # reservoir storage

            # Reservoir flag
            self.res_flag = 1  # 1 if with reservoir, 0 without reservoir
            # routing time step
            self.routing_timestep = 3 * 3600  # seconds

            # spin up run
            for nm in range(self.s.routing_spinup):

                sr = routing_mod.streamrouting(self.flow_dist,
                                               self.chs_prev,
                                               self.instream_flow,
                                               self.str_velocity,
                                               runoff[:, nm],
                                               self.calib_data.area,
                                               self.yr_imth_dys[nm, 2],
                                               self.routing_timestep,
                                               self.um,
                                               self.up,
                                               self.Sini,
                                               self.wdirr[:, nm],
                                               self.irrmean,
                                               self.mtifl_natural,
                                               self.ppose,
                                               self.cpa,
                                               self.HP_Release[:, :, np.mod(nm, 12)],
                                               self.maxtifl_natural,
                                               self.WConsumption[:, nm],
                                               self.s.alpha_param,
                                               self.res_prev,
                                               self.res_flag)

                (self.ChStorage[:, nm],
                 self.Avg_ChFlow[:, nm],
                 self.instream_flow,
                 self.Qin_Channel_avg,
                 self.Qout_channel_avg,
                 self.Qin_res_avg[:, nm],
                 self.Qout_res_avg[:, nm],
                 self.ResStorage[:, nm]) = sr

                # update data
                self.chs_prev = np.copy(self.ChStorage[:, nm])
                self.res_prev = np.copy(self.ResStorage[:, nm])

                # update the reservoir storage at beginning of year
                if np.mod(nm, 12) == 11:
                    self.Sini = self.ResStorage[:, nm]

            # simulation run
            for nm in range(self.nmonths):

                sr = routing_mod.streamrouting(self.flow_dist,
                                               self.chs_prev,
                                               self.instream_flow,
                                               self.str_velocity,
                                               runoff[:, nm],
                                               self.calib_data.area,
                                               self.yr_imth_dys[nm, 2],
                                               self.routing_timestep,
                                               self.um,
                                               self.up,
                                               self.Sini,
                                               self.wdirr[:, nm],
                                               self.irrmean,
                                               self.mtifl_natural,
                                               self.ppose,
                                               self.cpa,
                                               self.HP_Release[:, :, np.mod(nm, 12)],
                                               self.maxtifl_natural,
                                               self.WConsumption[:, nm],
                                               self.s.alpha_param,
                                               self.res_prev,
                                               self.res_flag)

                (
                    self.ChStorage[:, nm],
                    self.Avg_ChFlow[:, nm],
                    self.instream_flow,
                    self.Qin_Channel_avg,
                    self.Qout_channel_avg,
                    self.Qin_res_avg[:, nm],
                    self.Qout_res_avg[:, nm],
                    self.ResStorage[:, nm],
                ) = sr

                # update channel storage (chs) arrays for next step
                # update data
                self.chs_prev = np.copy(self.ChStorage[:, nm])
                self.res_prev = np.copy(self.ResStorage[:, nm])

                # update the reservoir storage at beginning of year
                if np.mod(nm, 12) == 11:
                    self.Sini = self.ResStorage[:, nm]

        return self.Avg_ChFlow[self.grdc_xanthosID, :]

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
                pet_out = np.zeros_like(self.data_abcd.precip)

                for nm in range(pet_num_steps):

                    # set up PET data for processing; hargreaves
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

            # TODO: check on this data reference
            AccessibleWater(self.s, self.data_reference, self.Q)

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

            Diagnostics(self.s, self.Q, DataDiagnostics())

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

        output_writer = OutWriter(self.s, self.data_reference.area, all_outputs)
        output_writer.write()

        try:
            self.q = output_writer.get('q')
        except ValueError:
            self.q = self.Q

        try:
            self.ac = output_writer.get('avgchflow')
        except ValueError:
            self.ac = self.Avg_ChFlow

        # TODO: check data reference
        output_writer.write_aggregates(self.data_reference, self.q, self.s.AggregateRunoffBasin, self.s.AggregateRunoffCountry,
                                       self.s.AggregateRunoffGCAMRegion)

        logging.info("---Output finished: %s seconds ---" % (time.time() - t0))

    def plots(self):
        """Create time series plots."""
        if self.s.CreateTimeSeriesPlot:
            logging.info("---Creating Time Series Plots:")
            t0 = time.time()

            # check data reference
            TimeSeriesPlot(self.s, self.q, self.ac, self.data_reference)

            logging.info("---Plots has finished successfully: %s seconds ------" % (time.time() - t0))

    def calibrate(self):
        """Run calibration to generate parameters for the ABCD model."""
        logging.info("---Processing PET...")
        t = time.time()

        pet_out = self.calculate_pet()

        logging.info("---PET processed in {} seconds---".format(time.time() - t))

        logging.info("---Loading calibration data:")
        calibration_data = DataCalibrationManaged(config_obj=self.s)

        logging.info("---Running calibration:")
        calib_managed.calibrate_all(settings=self.s,
                                    calibration_data=calibration_data,
                                    pet=pet_out)

        # calib_mod.calibrate_all(settings=self.s,
        #                         data=calibration_data,
        #                         pet=pet_out,
        #                         router_function=self.calculate_routing)
