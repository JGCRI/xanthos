import numpy as np
import time
from scipy.optimize import differential_evolution
from joblib import Parallel, delayed
from xanthos.runoff.abcd import ABCD


class Calibrate:

    def __init__(self, basin_num, basin_ids, basin_areas, precip, pet, obs, tmin, n_months, runoff_spinup, n_basins,
                 n_jobs, set_calibrate, obs_unit, out_dir, router_func=None):

        self.basin_num = basin_num
        self.basin_ids = basin_ids
        self.basin_areas = basin_areas
        self.precip = precip
        self.pet = pet
        self.obs = obs
        self.tmin = tmin
        self.n_months = n_months
        self.runoff_spinup = runoff_spinup
        self.router_func = router_func
        self.n_basins = n_basins
        self.n_jobs = n_jobs
        self.set_calibrate = set_calibrate
        self.obs_unit = obs_unit
        self.out_dir = out_dir

        # set the bounds; if the values are exactly 0 or 1 the model returns nan
        LB = 1e-4
        UB = 1 - LB
        self.bounds = [(LB, UB), (LB, 8 - LB), (LB, UB), (LB, UB), (LB, UB)]

        # create arrays to hold outputs
        self.all_pars = np.zeros((1, 5))
        self.kge_vals = np.zeros(1)

        # get grid indices for the current basin
        self.basin_idx = np.where(self.basin_ids == self.basin_num)
        self.bsn_areas = self.basin_areas[self.basin_idx]

        # transpose data for use in the ABCD model
        self.bsn_PET = self.pet[self.basin_idx]  # .T
        self.bsn_P = self.precip[self.basin_idx]  # .T
        self.bsn_TMIN = self.tmin[self.basin_idx]  # .T

        # select basin and remove extra years from observed runoff data
        self.bsn_Robs = self.obs[np.where(self.obs[:, 0] == basin_num)][:self.n_months, 1]

    def calibrate_basin(self, popsize=15):
        """
        This function is to calibrate the distributed ABCD model against the target global hydrological model (GHM) to
        obtain optimized parameters (a, b, c, d, m).

        :param basin_num:               basin number as an integer
        :param basin_ids:               an array of basin ids per grid cell that are associated with the basin
        :param basin_areas:             an array of basin areas per grid cell that are associated with the basin
        :param P:                       precipitation in mm/month
        :param PET:                     PET in mm/month
        :param Robs:                    Observed runoff in mm/month
        :param TMIN:                    minimum temperature in degrees C
        :param n_months:                the number of months in processing period
        :param spinup_steps:            the number of months from time 0 in the input data to use as spin up
        :param popsize:                 default 15; the number of random samples you are taking in each generation.  Higher
                                        should be better but takes longer
        :param polish:                  default False; if True, uses the LBGFT at end to see if result can be improved
        :return:
        """
        # record how long each basin takes to solve
        st = time.time()
        pars = differential_evolution(objective_kge,
                    bounds=self.bounds,
                    args=(basin_runoff, self.set_calibrate, self.bsn_PET, self.bsn_P, self.bsn_TMIN, self.n_months,
                          self.runoff_spinup, self.obs_unit, self.bsn_areas, self.bsn_Robs, self.basin_idx, self.precip.shape,
                          self.router_func),
                    popsize=popsize,
                    polish=False)

        # extract calibrated parameters for A,B,C,D,M; KGE values
        pars, ed, nfev = pars.x, pars.fun, pars.nfev

        self.all_pars[0, :] = pars
        self.kge_vals[0] = 1 - ed

        print("Finished calibration for basin {0} which contains {1} grid cells".format(self.basin_num, self.basin_idx[0].shape[0]))
        print("Population size:  {}".format(popsize))
        print("pars (a,b,c,d,m): {}".format(pars))
        print("KGE:  {}".format(1 - ed))
        print("Number of function evaluations:  {}".format(nfev))
        print("Calibration time (seconds):  {}".format(time.time() - st))

        np.save('{}/kge_result_basin_{}.npy'.format(self.out_dir, self.basin_num), self.kge_vals)
        np.save('{}/abcdm_parameters_basin_{}.npy'.format(self.out_dir, self.basin_num), self.all_pars)


def basin_runoff(pars, set_calibrate, pet, precip, tmin, n_months, runoff_spinup, obs_unit, bsn_areas, basin_idx, arr_shp, routing_func=None):

    # if calibrating against observed runoff
    if set_calibrate == 0:
        he = ABCD(pars, pet, precip, tmin, n_months, runoff_spinup, method='dist')
        he.emulate()

        if obs_unit == 'km3_per_mth':

            # convert from mm to km3
            return np.nansum(he.rsim * bsn_areas * 1e-6, 1)

        else:
            return he.rsim

    # if calibrating against observed streamflow
    else:
        he = ABCD(pars, pet, precip, tmin, n_months, runoff_spinup, method='dist')
        he.emulate()

        # add rsim data from the basin back into a global array of zeros to be used by the router
        rsim = np.zeros(shape=arr_shp)
        np.put(rsim, basin_idx, he.rsim)

        return routing_func(rsim)


def objective_kge(pars, model_func, set_calibrate, pet, precip, tmin, n_months, runoff_spinup, obs_unit, bsn_areas, bsn_Robs, basin_idx, arr_shp, routing_func=None):
    """
    Kling-Gupta efficiency between simulated and observed.

    :param pars:                array of A, B, C, D, M parameters for the ABCD model
    :param P:                   precipitation
    :param PET:                 PET
    :param Robs:                Observed runoff
    :param TMIN:                minimum temerature
    :param grid_sizes:          grid cell area per grid
    :param n_months:            number of processing months
    :param spinup_steps:        number of months from time 0 of the input data to use as spin up
    :return efficiency:         KGE
    """
    modelled = model_func(pars, set_calibrate, pet, precip, tmin, n_months, runoff_spinup, obs_unit, bsn_areas, basin_idx, arr_shp, routing_func)
    observed = bsn_Robs

    # Calculate KGE
    sd_modelled = np.std(modelled)
    sd_observed = np.std(observed)
    m_modelled = np.mean(modelled)
    m_observed = np.mean(observed)

    # alpha
    relvar = sd_modelled / sd_observed

    # beta
    bias = m_modelled / m_observed

    # r
    corrcoef = np.corrcoef(observed, modelled)[1, 0]

    ed = (((corrcoef - 1)**2) + ((relvar - 1)**2) + ((bias - 1)**2))**0.5

    return ed

# def objective_kge(pars, P, PET, Robs, TMIN, grid_sizes, n_months, spinup_steps):
#     """
#     Kling-Gupta efficiency between simulated and observed.
#
#     :param pars:                array of A, B, C, D, M parameters for the ABCD model
#     :param P:                   precipitation
#     :param PET:                 PET
#     :param Robs:                Observed runoff
#     :param TMIN:                minimum temerature
#     :param grid_sizes:          grid cell area per grid
#     :param n_months:            number of processing months
#     :param spinup_steps:        number of months from time 0 of the input data to use as spin up
#     :return efficiency:         KGE
#     """
#     # instantiate and run the model
#     he = ABCD(pars, PET, P, TMIN, n_months, spinup_steps, method='dist')
#     he.emulate()
#
#     Rsim = np.nansum(he.rsim * grid_sizes * 1e-6, 1) # convert from mm to km3
#
#     # Use simulations of the last 20 years for model performance evaluation
#     modelled = Rsim
#     observed = Robs
#
#     # Calculate KGE
#     sd_modelled = np.std(modelled)
#     sd_observed = np.std(observed)
#     m_modelled = np.mean(modelled)
#     m_observed = np.mean(observed)
#
#     # alpha
#     relvar = sd_modelled / sd_observed
#
#     # beta
#     bias = m_modelled / m_observed
#
#     # r
#     corrcoef = np.corrcoef(observed, modelled)[1,0]
#
#     ed = (((corrcoef - 1)**2) + ((relvar - 1)**2) + ((bias - 1)**2))**0.5
#
#     return ed


# def calibrate_basin(basin_num, basin_ids, basin_areas, P, PET, Robs, TMIN, n_months, spinup_factor, popsize=15, out_file=None):
#     """
#     This function is to calibrate the distributed ABCD model against the target global hydrological model (GHM) to
#     obtain optimized parameters (a, b, c, d, m).
#
#     :param basin_num:               basin number as an integer
#     :param basin_ids:               an array of basin ids per grid cell that are associated with the basin
#     :param basin_areas:             an array of basin areas per grid cell that are associated with the basin
#     :param P:                       precipitation in mm/month
#     :param PET:                     PET in mm/month
#     :param Robs:                    Observed runoff in mm/month
#     :param TMIN:                    minimum temperature in degrees C
#     :param n_months:                the number of months in processing period
#     :param spinup_steps:            the number of months from time 0 in the input data to use as spin up
#     :param popsize:                 default 15; the number of random samples you are taking in each generation.  Higher
#                                     should be better but takes longer
#     :param polish:                  default False; if True, uses the LBGFT at end to see if result can be improved
#     :param out_file:
#     :return:
#     """
#     # set the bounds; if the values are exactly 0 or 1 the model returns nan
#     LB = 1e-4
#     UB = 1 - LB
#     bounds = [(LB, UB), (LB, 8-LB), (LB, UB), (LB, UB), (LB, UB)]
#
#     # create arrays to hold outputs
#     all_pars = np.zeros((1, 5))
#     kge_vals = np.zeros(1)
#
#     # get grid indices for the current basin
#     basin_idx = np.where(basin_ids == basin_num)
#     bsn_areas = basin_areas[basin_idx]
#
#     # transpose data for use in the ABCD model
#     bsn_PET = PET[basin_idx]#.T
#     bsn_P = P[basin_idx]#.T
#     bsn_TMIN = TMIN[basin_idx]#.T
#
#     # select basin and remove extra years from observed runoff data
#     bsn_Robs = Robs[np.where(Robs[:, 0] == basin_num)]
#     bsn_Robs = bsn_Robs[:n_months, 1]
#
#     # record how long each basin takes to solve
#     st = time.time()
#     pars = differential_evolution(objective_kge,
#                 bounds=bounds,
#                 args=(bsn_P, bsn_PET, bsn_Robs, bsn_TMIN, bsn_areas, n_months, spinup_factor),
#                 popsize=popsize,
#                 polish=False)
#
#     pars, ed, nfev = pars.x, pars.fun, pars.nfev
#
#     all_pars[0, :] = pars
#     kge_vals[0] = 1 - ed
#
#     print("Finished calibration for basin {0} which contains {1} grid cells".format(basin_num, basin_idx[0].shape[0]))
#     print("Population size:  {}".format(popsize))
#     print("pars (a,b,c,d,m): {}".format(pars))
#     print("KGE:  {}".format(1 - ed))
#     print("Number of function evaluations:  {}".format(nfev))
#
#     # out_dir = '/pic/projects/GCAM/gcam_hydrology/process/outputs'
#     out_dir = '/users/d3y010/Desktop'
#
#     np.save('{}/kge_diffevo_basin_{}.npy'.format(out_dir, basin_num), kge_vals)
#
#     np.save('{}/outpars_diffevo_basin_{}.npy'.format(out_dir, basin_num), all_pars)


def calibrate(basin_ids, basin_areas, precip, pet, obs, tmin, n_months, spinup_steps, n_basins, n_jobs):
    """
    Run calibration for ABCD model parameters.

    :param basin_num:
    :param basin_ids:
    :param basin_areas:
    :param precip:
    :param pet:
    :param obs:
    :param tmin:
    :param n_months:
    :param spinup_steps:
    :param n_basins:
    :return:
    """

    Parallel(n_jobs=n_jobs)(delayed(calibrate_basin)(basin_num, basin_ids, basin_areas, precip, pet, obs, tmin,
                                                     n_months, spinup_steps) for basin_num in range(1, n_basins + 1, 1))