"""
Calibrate the ABCD model.

@author   Caleb Braun, Chris R. Vernon
@email:   caleb.braun@pnnl.gov
@Project: Xanthos 2.0

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2018, Battelle Memorial Institute
"""

import os
import numpy as np
import logging
import time
from scipy.optimize import differential_evolution
import spotpy
from xanthos.runoff.abcd import ABCD

from xanthos.data_reader.data_calibration import DataCalibration
from xanthos.data_reader.data_reference import DataReference

import xanthos.routing.mrtm as routing_mod


import xanthos.utils.general as helper
import xanthos.utils.math as umth


class Calibrate:
    """Calibrate the ABCD runoff module."""

    NCELL = 67420
    NGRIDROW = 360
    NGRIDCOL = 720

    LB = 1e-4
    UB = 1 - LB
    LB1 = 0.75
    LB2 = 1
    LBA = 9.9
    UB2 = 10.0

    def __init__(self, basin_num, basin_ids, basin_areas, precip, pet, obs,
                 tmin, runoff_spinup, set_calibrate, obs_unit,
                 out_dir, router_func=None, config_obj=None, cal_observed=None, purpose_file=None,
                 capacity_file=None, hp_release_file=None,
                 water_consumption_file=None, instream_natural_flow_file=None,
                 initial_channel_storage_natural_file=None, sm_file=None, mtif_natural_file=None,
                 maxtif_natural_file=None, total_demand_cumecs_file=None,
                 grdc_xanthos_coord_index_file=None, start_year=None, end_year=None, gen=100, nchains=10):
        """Initialize calibration data and parameters.

        :param basin_num:      basin number as an integer
        :param basin_ids:      an array of basin ids per grid cell that are associated with the basin
        :param basin_areas:    an array of basin areas per grid cell that are associated with the basin
        :param precip:         precipitation in mm/month
        :param pet:            PET in mm/month
        :param obs:            Observed runoff in mm/month
        :param tmin:           minimum temperature in degrees C
        :param runoff_spinup:  the number of months from time 0 in the input data to use as spin up
        :param set_calibrate:  0 to calibrate to observed runoff, 1 to calibrate to observed streamflow
        :param obs_unit:       the unit of the input data
        :param out_dir:        calibrated parameters output directory
        :param router_func:    objective function for calibrating routing

        """

        # load calibration data
        self.calib_data = DataCalibration(config_obj=config_obj, cal_observed=cal_observed, purpose_file=purpose_file,
                                          capacity_file=capacity_file, hp_release_file=hp_release_file,
                                          water_consumption_file=water_consumption_file,
                                          instream_natural_flow_file=instream_natural_flow_file,
                                          initial_channel_storage_natural_file=initial_channel_storage_natural_file,
                                          sm_file=sm_file, mtif_natural_file=mtif_natural_file,
                                          maxtif_natural_file=maxtif_natural_file,
                                          total_demand_cumecs_file=total_demand_cumecs_file,
                                          grdc_xanthos_coord_index_file=grdc_xanthos_coord_index_file)

        if config_obj is None:
            self.start_year = start_year
            self.end_year = end_year
            self.nmonths = (self.end_year - self.start_year + 1) * 12

            # load reference data
            self.reference_data = DataReference(nmonths=self.nmonths)

        else:
            self.start_year = self.config_obj.StartYear
            self.end_year = self.config_obj.EndYear
            self.nmonths = self.config_obj.nmonths
            self.reference_data = DataReference(config=config_obj)

        self.nmonths = (self.end_year - self.start_year + 1) * 12
        self.basin_num = basin_num
        self.basin_ids = basin_ids
        self.basin_areas = basin_areas
        self.precip = precip
        self.pet = pet
        self.obs = obs
        self.tmin = tmin
        self.runoff_spinup = runoff_spinup
        self.router_func = router_func
        self.set_calibrate = set_calibrate
        self.obs_unit = obs_unit
        self.out_dir = out_dir
        self.gen = gen
        self.nChains = nchains

        # Minimum temperature is optional; if not provided, the snow components
        # of the model is effectively removed, so remove the model parameter
        # for snow (M)
        self.nosnow = self.tmin is None

        # index for gauge station locations
        self.grdcData_info = np.copy(self.calib_data.grdc_coord_index)
        self.grdc_xanthosID = self.grdcData_info[np.where(self.grdcData_info[:, 0] == self.basin_num), 1][0][0] - 1

        # routing inputs: wdirr, irrmean, tifl, ppose, cpa,dscells
        self.wdirr = np.copy(self.calib_data.total_demand_cumecs)
        self.irrmean = np.mean(self.wdirr, axis=1)  # mean demand
        self.ppose = np.copy(self.calib_data.purpose)
        self.cpa = np.copy(self.calib_data.capacity) * 10**6  # m3
        self.HP_Release = np.copy(self.calib_data.hp_release)
        self.WConsumption = np.copy(self.calib_data.water_consumption)
        self.chs_ini = np.copy(self.calib_data.ini_channel_storage)
        self.Initial_instream_flow = np.copy(self.calib_data.instream_flow_natural)
        self.SM = np.squeeze(np.copy(self.calib_data.sm))
        self.mtifl_natural = np.copy(self.calib_data.mtif_natural)
        self.maxtifl_natural = np.copy(self.calib_data.maxtif_natural)

        # routing data
        self.yr_imth_dys = helper.set_month_arrays(self.nmonths, self.start_year, self.end_year)
        self.map_index = umth.sub2ind([Calibrate.NGRIDROW, Calibrate.NGRIDCOL],
                                      self.reference_data.coords[:, 4].astype(int) - 1,
                                      self.reference_data.coords[:, 3].astype(int) - 1)

        self.routing_data = DataMrtm()

        # set number of parameter combinations
        self.ModelPerformance = os.path.join(self.out_dir, 'Basins_Result', f"basin_calibration_{self.basin_num}")

        # set the bounds; if the values are exactly 0 or 1 the model returns nan
        self.bounds = [(Calibrate.LB, Calibrate.UB), (Calibrate.LB, 8 - Calibrate.LB), (Calibrate.LB, Calibrate.UB),
                       (Calibrate.LB, Calibrate.UB), (Calibrate.LB, Calibrate.UB), (Calibrate.LBA, Calibrate.UB2),
                       (Calibrate.LB1, Calibrate.UB)]

        if self.nosnow:
            self.bounds.pop()  # remove calibration parameter M

        # set up parameters
        self.params = np.zeros((1, len(self.bounds)))

        # get grid indices for the current basin
        self.basin_idx = np.where(self.basin_ids == self.basin_num)
        self.bsn_areas = self.basin_areas[self.basin_idx]

        # transpose data for use in the ABCD model
        self.bsn_PET = self.pet[self.basin_idx]
        self.bsn_P = self.precip[self.basin_idx]

        # if no tmin provided, just ensure it is larger than the rain threshold
        if self.nosnow:
            self.bsn_TMIN = None
        else:
            self.bsn_TMIN = self.tmin[self.basin_idx]

        # initial soil moisture
        self.bsn_SM = self.SM[self.basin_idx]

        # Unit conversion for runoff case
        if self.obs_unit == "km3_per_mth":
            self.conversion = self.bsn_areas * 1e-6
        elif self.obs_unit == "mm_per_mth":
            self.conversion = 1.0

        # select basin and remove extra years from observed runoff data
        self.bsn_obs = self.obs[np.where(self.obs[:, 0] == self.basin_num)][: self.nmonths, 1]
        self.bsn_Robs_calib = self.bsn_obs[0:120]
        self.bsn_Robs_valid = self.bsn_obs[120:240]

        # residence time in hr
        Lst = self.rt_data.flow_dist
        Vst = self.rt_data.str_velocity
        Vst[Vst < 0.01] = 0.01
        grid_size = np.sqrt(self.GridConstants.area) * 1000
        nn = np.where(Lst < grid_size)[0]
        Lst[nn] = grid_size[nn]

        # residence time in hr
        self.Tr = np.divide(Lst, Vst) / 3600

        # best parameter from preceeding simulations
        self.best_params = None

    def bestParams_combination(self):
        self.best_params = None
        if os.path.isfile(self.ModelPerformance + ".csv"):
            results = get_calib_data(self.ModelPerformance + ".csv")
            ObjFuns = 1
            try:
                best = results[np.nanargmax([r[0] for r in results])]
                self.kge_cal = best[0]

                x = [best[p] for p in range(ObjFuns, ObjFuns + len(self.bounds))]
                self.best_params = np.array(x)
            except ValueError:
                pass
            return self.best_params.astype(np.float)

    def parameters(self):
        '''Returns ABCD Params'''
        if os.path.isfile(self.ModelPerformance + ".csv"):
            if not os.stat(self.ModelPerformance + ".csv").st_size == 0:
                self.best_params = self.bestParams_combination()

            else:
                self.best_params = None

        names = [0, 1, 2, 3, 4, 5, 6]  # parameters order
        Tr_basinmin = min(self.Tr[self.basin_idx])
        alpha = 0.95 * Tr_basinmin / 3
        if alpha < self.bounds[5][1]:
            self.bounds[5] = (1, alpha)

        if self.best_params is None:
            params = [spotpy.parameter.Uniform(
                self.bounds[p][0], self.bounds[p][1])
                for p in names]
        else:
            params = [spotpy.parameter.Uniform(
                self.bounds[p][0], self.bounds[p][1], optguess=self.best_params[p])
                for p in names]

        return spotpy.parameter.generate(params)

    # ABCD model and mrtm routing model : this function provides simulated streamflow
    def simulation(self, pars):
        if self.set_calibrate == 0:
            ncell_in_basin = len(self.basin_idx[0])
            self.target_basin_ids = np.zeros(ncell_in_basin)
            he = ABCD(
                pars,
                self.bsn_SM,
                self.bsn_PET,
                self.bsn_P,
                self.bsn_TMIN,
                self.target_basin_ids,
                self.n_months,
                self.runoff_spinup,
                method="dist",
            )
            he.emulate()
            self.qsim_with_validation = np.nansum(he.rsim[0:240, :] * self.conversion, 1)
            return np.nansum(he.rsim[0:120, :] * self.conversion, 1)
        else:
            ncell_in_basin = len(self.basin_ids)
            self.target_basin_ids = np.zeros(ncell_in_basin)
            he = ABCD(pars,
                      self.SM,
                      self.pet,
                      self.precip,
                      self.tmin,
                      self.target_basin_ids,
                      self.n_months,
                      self.runoff_spinup,
                      method="dist",
                      )
            he.emulate()

            self.runoff = he.rsim.T
            runoff_out = self.runoff[self.basin_idx, :]

            # load routing data

            self.flow_dist = self.rt_data.flow_dist
            grid_size = np.sqrt(self.GridConstants.area) * 1000
            nn_grids = np.where(self.flow_dist < grid_size)[0]
            self.flow_dist[nn_grids] = grid_size[nn_grids]

            self.flow_dir = self.rt_data.flow_dir
            self.str_v = (self.rt_data.str_velocity)
            self.str_v[self.str_v < 0.01] = 0.01
            self.str_velocity = self.str_v * pars[5]
            self.chs_prev = self.rt_data.chs_prev

            self.dsid = routing_mod.downstream(self.rt_data.coords, self.flow_dir, self.s)
            self.upid = routing_mod.upstream(self.rt_data.coords, self.dsid, self.s)
            self.um, self.up = routing_mod.upstream_genmatrix(self.upid)

            # routing initializations
            self.chs_prev = np.squeeze(self.chs_ini)  # initial channel storage in m3
            self.instream_flow = np.squeeze(self.Initial_instream_flow)  # initial channel flow in m3 /s
            self.Sini = 0.5 * np.squeeze(self.cpa)  # initial reservoir storage in m3
            self.res_prev = np.squeeze(self.cpa)  # reservoir storage from t-1 in m3
            # self.instream_flow = self.rt_data.instream_flow
            # self.chs_prev = self.rt_data.chs_prev

            # Preallocation for variables to be returned from routing module
            self.Avg_ChFlow = np.zeros_like(self.precip)  # average channel flow m3/s
            self.ChStorage = np.zeros_like(self.precip)  # channel storage
            self.ResStorage = np.zeros_like(self.precip)  # reservoir storage
            self.Qout_res_avg = np.zeros_like(self.precip)  # reservoir storage
            self.Qin_res_avg = np.zeros_like(self.precip)  # reservoir storage

            # Reservoir flag
            self.res_flag = 1  # 1 if with reservoir, 0 without reservoir
            # routing time step
            self.routing_timestep = 3 * 3600  # seconds
            #############################################################################
            # spin up run
            # for nm in range(self.s.routing_spinup):
            mm_month = 120
            for nm in range(mm_month):
                sr = routing_mod.streamrouting(self.flow_dist,
                                               self.chs_prev,
                                               self.instream_flow,
                                               self.str_velocity,
                                               self.runoff[:, nm],
                                               self.basin_areas,
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
                                               pars[6],
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

                # update data
                self.chs_prev = np.copy(self.ChStorage[:, nm])
                self.res_prev = np.copy(self.ResStorage[:, nm])

                # update the reservoir storage at beginning of year
                if np.mod(nm, 12) == 11:
                    self.Sini = self.ResStorage[:, nm]

            # simulation run
            mm_month = 240
            for nm in range(mm_month):
                # for nm in range(self.s.nmonths ):
                sr = routing_mod.streamrouting(self.flow_dist,
                                               self.chs_prev,
                                               self.instream_flow,
                                               self.str_velocity,
                                               self.runoff[:, nm],
                                               self.basin_areas,
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
                                               pars[6],
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

        self.qsim_with_validation = self.Avg_ChFlow[self.grdc_xanthosID, 0:240]
        return self.Avg_ChFlow[self.grdc_xanthosID, 0:120]

    # Objective function to be minimized (if sceua is used) and maximized (all others)
    @staticmethod
    def objectivefunction(simulation, evaluation):
        """Calculates Model Performance"""
        kge = spotpy.objectivefunctions.kge(evaluation, simulation)
        return kge

    # observed streamflow data
    def evaluation(self):
        return self.bsn_Robs_calib

    def save(self, objectivefunctions, parameter, simulations):
        line = str(objectivefunctions) + ',' + str(parameter).strip('[]') + ',' + str(simulations).strip(
            '[]') + '\n'
        self.database.write(line)

    # calibration set up
    def calibrate_basin(self):
        """
        This function is to calibrate the distributed ABCD + water management model against the GRDC to
        obtain optimized parameters of ABCD(a, b, c, d, m) and Water management (beta and c)

        """
        # self.best_params = self.bestParams_combination()
        sampler = spotpy.algorithms.demcz(
            self, dbname=self.ModelPerformance, dbformat="csv", dbappend=True,
            save_sim=False)  # , parallel='mpi')#demcz  sceua

        sampler.sample(self.gen)  # ,nChains=self.nChains)
        print('result with optimal params')
        self.best_params = self.bestParams_combination()
        # self.simulation(self.best_params)
        qsimulated, qobserved, kge_cal, kge_val = self.calibration_run(self.best_params)

    def calibration_run(self, x):
        qobs_cal = self.bsn_obs[0:120]
        qobs_val = self.bsn_obs[121:240]
        qsim_cal = self.simulation(x)
        qsimulated = self.qsim_with_validation
        kge_cal = spotpy.objectivefunctions.kge(qobs_cal, qsimulated[0:120])
        kge_val = spotpy.objectivefunctions.kge(qobs_val, qsimulated[121:240])

        print("kge_cal: {}, kge_val: {}".format(kge_cal, kge_val))

        return qsimulated, self.bsn_obs, kge_cal, kge_val


def process_basin(basin_num, settings, data, pet, router_function=None):
    """
    Process single basin.
    """

    cal = Calibrate(basin_num=settings.basin_num,
                    set_calibrate=settings.set_calibrate,
                    obs_unit=settings.obs_unit,
                    basin_ids=data.basin_ids,
                    basin_areas=data.area,
                    precip=data.precip,
                    pet=pet,
                    obs=data.cal_obs,
                    tmin=data.tmin,
                    n_months=settings.nmonths,
                    runoff_spinup=settings.runoff_spinup,
                    router_func=router_function,
                    out_dir=settings.calib_out_dir
                    )
    cal.calibrate_basin()


def calibrate_all(settings, data, pet):
    """Run calibration for ABCD model for a basins."""
    print("\tCalibrating Basin:  {}".format(settings.basin_num))
    process_basin(settings.basin_num, settings, data, pet, router_function=None)


def get_calib_data(performance_file):
    """Read in HDF data"""
    import pandas as pd
    perf_calib11 = (pd.read_csv(performance_file, skiprows=1)).to_numpy()
    perf_calib = perf_calib11[0:perf_calib11.shape[0], :]
    perf_calib[:, 0] = perf_calib[:, 0]
    return perf_calib
