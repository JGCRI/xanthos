"""
CalibrateManaged the ABCD model using water management rules.

@authors: HongYi Li : hli57@uh.edu,
         University of Houston
         Guta Abeshu : gwabeshu@uh.edu
         University of Houston

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2018, Battelle Memorial Institute
"""

import os
import spotpy
import numpy as np
import pandas as pd
import matplotlib.pylab as plt

from xanthos.runoff.abcd_managed import AbcdManaged
from xanthos.data_reader.data_mrtm_managed import DataMrtmManaged
from xanthos.data_reader.data_calibration_managed import DataCalibrationManaged
from xanthos.data_reader.data_reference import DataReference
from xanthos.data_reader.data_abcd import DataAbcd

import xanthos.routing.mrtm_managed as routing_mod
import xanthos.utils.general as helper
import xanthos.utils.math as umth


class CalibrateManaged:
    """CalibrateManaged the ABCD runoff module.

    :params repetitions:                    maximum number of function evaluations allowed during optimization
    :type repetitions:                      int

    """

    NCELL = 67420
    NGRIDROW = 360
    NGRIDCOL = 720

    LB = 1e-4
    UB = 1 - LB
    LB1 = 0.5
    LB2 = 1
    LBA = 1.0
    UB2 = 10.0

    def __init__(self,
                 basin_num,
                 basin_ids,
                 basin_areas,
                 precip,
                 pet,
                 obs,
                 tmin,
                 runoff_spinup,
                 set_calibrate,
                 obs_unit,
                 out_dir,
                 nmonths=None,
                 router_func=None,
                 config_obj=None,
                 cal_observed=None,
                 purpose_file=None,
                 capacity_file=None,
                 hp_release_file=None,
                 water_consumption_file=None,
                 instream_flow_natural_file=None,
                 initial_channel_storage_natural_file=None,
                 sm_file=None,
                 mtif_natural_file=None,
                 maxtif_natural_file=None,
                 total_demand_cumecs_file=None,
                 grdc_coord_index_file=None,
                 start_year=None,
                 end_year=None,
                 repetitions=100,
                 nchains=10,
                 flow_distance_file=None,
                 flow_direction_file=None,
                 stream_velocity_file=None,
                 historical_mode="True",
                 hist_channel_storage_file=None,
                 hist_channel_storage_varname=None,
                 routing_spinup=None):
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

        if config_obj is None:
            self.start_year = start_year
            self.end_year = end_year
            
            if nmonths is None:
                self.nmonths = (self.end_year - self.start_year + 1) * 12
            else:
                self.nmonths = nmonths

            # load reference data
            self.reference_data = DataReference(nmonths=self.nmonths)

            self.flow_distance_file = flow_distance_file
            self.flow_direction_file = flow_direction_file
            self.stream_velocity_file = stream_velocity_file
            self.historical_mode = historical_mode
            self.hist_channel_storage_file = hist_channel_storage_file
            self.hist_channel_storage_varname = hist_channel_storage_varname
            self.routing_spinup = routing_spinup
            self.repetitions = repetitions

        else:
            self.start_year = config_obj.StartYear
            self.end_year = config_obj.EndYear
            self.nmonths = config_obj.nmonths
            self.reference_data = DataReference(config=config_obj)
            self.routing_spinup = config_obj.routing_spinup
            self.repetitions = config_obj.repetitions

            self.flow_distance_file = config_obj.flow_distance
            self.flow_direction_file = config_obj.flow_direction
            self.stream_velocity_file = config_obj.strm_veloc
            self.historical_mode = config_obj.HistFlag

            self.hist_channel_storage_file = config_obj.ChStorageFile
            self.hist_channel_storage_varname = config_obj.ChStorageVarName

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
        self.nChains = nchains

        # load calibration data
        self.calib_data = DataCalibrationManaged(config_obj=config_obj,
                                                 cal_observed=cal_observed,
                                                 purpose_file=purpose_file,
                                                 capacity_file=capacity_file,
                                                 hp_release_file=hp_release_file,
                                                 water_consumption_file=water_consumption_file,
                                                 instream_flow_natural_file=instream_flow_natural_file,
                                                 initial_channel_storage_natural_file=initial_channel_storage_natural_file,
                                                 sm_file=sm_file,
                                                 mtif_natural_file=mtif_natural_file,
                                                 maxtif_natural_file=maxtif_natural_file,
                                                 total_demand_cumecs_file=total_demand_cumecs_file,
                                                 grdc_coord_index_file=grdc_coord_index_file,
                                                 start_year=self.start_year,
                                                 end_year=self.end_year)

        # Minimum temperature is optional; if not provided, the snow components
        # of the model is effectively removed, so remove the model parameter
        # for snow (M)
        self.nosnow = self.tmin is None

        # index for gauge station locations
        self.grdcData_info = np.copy(self.calib_data.grdc_coord_index_file)
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
        self.map_index = umth.sub2ind([CalibrateManaged.NGRIDROW, CalibrateManaged.NGRIDCOL],
                                      self.reference_data.coords[:, 4].astype(int) - 1,
                                      self.reference_data.coords[:, 3].astype(int) - 1)

        self.routing_data = DataMrtmManaged(start_year=self.start_year,
                                            end_year=self.end_year,
                                            flow_distance_file=self.flow_distance_file,
                                            flow_direction_file=self.flow_direction_file,
                                            stream_velocity_file=self.stream_velocity_file,
                                            historical_mode=self.historical_mode,
                                            hist_channel_storage_file=self.hist_channel_storage_file,
                                            hist_channel_storage_varname=self.hist_channel_storage_varname)

        # set number of parameter combinations
        self.ModelPerformance = os.path.join(self.out_dir, f"basin_calibration_{self.basin_num}")

        # set the bounds; if the values are exactly 0 or 1 the model returns nan
        self.bounds = [(CalibrateManaged.LB, CalibrateManaged.UB),
                       (CalibrateManaged.LB, 8 - CalibrateManaged.LB),
                       (CalibrateManaged.LB, CalibrateManaged.UB),
                       (CalibrateManaged.LB, CalibrateManaged.UB),
                       (CalibrateManaged.LB, CalibrateManaged.UB),
                       (CalibrateManaged.LBA, CalibrateManaged.UB2),
                       (CalibrateManaged.LB1, CalibrateManaged.UB)]

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
        Lst = self.routing_data.flow_dist
        Vst = self.routing_data.str_velocity
        Vst[Vst < 0.01] = 0.01
        grid_size = np.sqrt(self.reference_data.area) * 1000
        nn = np.where(Lst < grid_size)[0]
        Lst[nn] = grid_size[nn]
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
            params = [spotpy.parameter.Uniform(self.bounds[p][0], self.bounds[p][1], optguess=self.best_params[p]) for p in names]

        return spotpy.parameter.generate(params)

    def simulation(self, pars):
        """ABCD model and mrtm routing model : this function provides simulated streamflow"""

        if self.set_calibrate == 0:

            ncell_in_basin = len(self.basin_idx[0])

            self.target_basin_ids = np.zeros(ncell_in_basin)

            he = AbcdManaged(pars=pars,
                             soil_water_initial=self.bsn_SM,
                             pet=self.bsn_PET,
                             precip=self.bsn_P,
                             tmin=self.bsn_TMIN,
                             basin_ids=self.target_basin_ids,
                             process_steps=self.nmonths,
                             spinup_steps=self.runoff_spinup,
                             method="dist")

            he.emulate()

            self.qsim_with_validation = np.nansum(he.rsim[0:240, :] * self.conversion, 1)

            return np.nansum(he.rsim[0:120, :] * self.conversion, 1)

        else:

            ncell_in_basin = len(self.basin_ids)

            self.target_basin_ids = np.zeros(ncell_in_basin)

            he = AbcdManaged(pars=pars,
                             soil_water_initial=self.SM,
                             pet=self.pet,
                             precip=self.precip,
                             tmin=self.tmin,
                             basin_ids=self.target_basin_ids,
                             process_steps=self.nmonths,
                             spinup_steps=self.runoff_spinup,
                             method="dist")

            he.emulate()

            self.runoff = he.rsim.T

            runoff_out = self.runoff[self.basin_idx, :]

            # load routing data
            self.flow_dist = self.routing_data.flow_dist
            grid_size = np.sqrt(self.reference_data.area) * 1000
            nn_grids = np.where(self.flow_dist < grid_size)[0]
            self.flow_dist[nn_grids] = grid_size[nn_grids]

            self.flow_dir = self.routing_data.flow_dir
            self.str_v = (self.routing_data.str_velocity)
            self.str_v[self.str_v < 0.01] = 0.01
            self.str_velocity = self.str_v * pars[5]
            self.chs_prev = self.routing_data.chs_prev

            self.dsid = routing_mod.downstream(self.reference_data.coords, self.flow_dir, CalibrateManaged.NGRIDROW,
                                               CalibrateManaged.NGRIDCOL)
            self.upid = routing_mod.upstream(self.reference_data.coords, self.dsid, CalibrateManaged.NGRIDROW,
                                               CalibrateManaged.NGRIDCOL)
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

            # spin up run
            for nm in range(self.routing_spinup):

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


    @staticmethod
    def objectivefunction(simulation, evaluation, method='sceua'):
        """Calculates Model Performance.
        Objective function to be minimized (if sceua is used) and maximized (all others)

        """
        # sceua requires minimization which will result in a negative KGE
        if method == 'sceua':
            multiplier = -1
        else:
            multiplier = 1

        return spotpy.objectivefunctions.kge(evaluation, simulation) * multiplier

    def evaluation(self):
        """observed streamflow data"""

        return self.bsn_Robs_calib

    def save(self, objectivefunctions, parameter, simulations):
        line = str(objectivefunctions) + ',' + str(parameter).strip('[]') + ',' + str(simulations).strip('[]') + '\n'
        self.database.write(line)

    # calibration set up
    def calibrate_basin(self):
        """This function is to calibrate the distributed ABCD + water management model against the GRDC to
        obtain optimized parameters of ABCD(a, b, c, d, m) and Water management (beta and c)

        """
        sampler = spotpy.algorithms.sceua(self,
                                          dbname=self.ModelPerformance,
                                          dbformat="csv",
                                          dbappend=False,
                                          save_sim=False)

        sampler.sample(self.repetitions, ngs=50, kstop=100, peps=1e-4, pcento=1e-4)

        print('Result using optimal parameters:')
        self.best_params = self.bestParams_combination()

        print(self.best_params)

        qsimulated, qobserved, kge_cal, kge_val = self.calibration_run(self.best_params)

    def calibration_run(self, x):

        qobs_cal = self.bsn_obs[0:120]
        qobs_val = self.bsn_obs[121:240]
        qsim_cal = self.simulation(x)
        qsimulated = self.qsim_with_validation

        # KGE of the calibration period
        kge_cal = spotpy.objectivefunctions.kge(qobs_cal, qsimulated[0:120])

        # KGE of the validation period
        kge_val = spotpy.objectivefunctions.kge(qobs_val, qsimulated[121:240])

        print("Calibration KGE: {}, Validation KGE: {}".format(kge_cal, kge_val))

        return qsimulated, self.bsn_obs, kge_cal, kge_val


def process_basin(config_obj, calibration_data, pet, router_function=None):
    """Process single basin."""

    # load ABCD runoff module data
    data_abcd = DataAbcd(config=config_obj)

    cal = CalibrateManaged(config_obj=config_obj,
                           basin_num=config_obj.basin_num,
                           set_calibrate=config_obj.set_calibrate,
                           obs_unit=config_obj.obs_unit,
                           basin_ids=calibration_data.basin_ids,
                           basin_areas=calibration_data.area,
                           precip=data_abcd.precip,
                           pet=pet,
                           obs=calibration_data.cal_obs,
                           tmin=data_abcd.tmin,
                           nmonths=config_obj.nmonths,
                           runoff_spinup=config_obj.runoff_spinup,
                           router_func=router_function,
                           out_dir=config_obj.calib_out_dir,
                           start_year=config_obj.StartYear,
                           end_year=config_obj.EndYear
                           )
    cal.calibrate_basin()


def calibrate_all(settings, calibration_data, pet, router_fn=None):
    """Run calibration for ABCD model for a basins."""

    print("\tCalibrating Basin:  {}".format(settings.basin_num))

    process_basin(settings, calibration_data, pet, router_function=router_fn)


def get_calib_data(performance_file):
    """Read in HDF data"""

    perf_calib11 = (pd.read_csv(performance_file, skiprows=1)).to_numpy()
    perf_calib = perf_calib11[0:perf_calib11.shape[0], :]
    perf_calib[:, 0] = perf_calib[:, 0]
    
    return perf_calib


def plot_kge(calibration_result_file, output_file_name, dpi=300, figsize=(9, 5)):
    """Plot the KGE result of a calibrated basin"""

    # load results
    results = spotpy.analyser.load_csv_results(calibration_result_file)

    # create plot
    fig = plt.figure(1, figsize=figsize)
    plt.plot(results['like1'] * -1)
    plt.ylabel('KGE')
    plt.xlabel('Repetition')
    fig.savefig(output_file_name, dpi=dpi)

    return plt
