import numpy as np
import time
from scipy.optimize import differential_evolution

from xanthos.runoff.abcd import ABCD


class Calibrate:

    def __init__(self, basin_num, basin_ids, basin_areas, precip, pet, obs, tmin, n_months, runoff_spinup,
                 set_calibrate, obs_unit, out_dir, router_func=None):

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
        self.set_calibrate = set_calibrate
        self.obs_unit = obs_unit
        self.out_dir = out_dir

        # Minimum temperature is optional; if not provided, the snow components
        # of the model is effectively removed, so remove the model parameter
        # for snow (M)
        self.nosnow = self.tmin is None

        # set the bounds; if the values are exactly 0 or 1 the model returns nan
        LB = 1e-4
        UB = 1 - LB
        self.bounds = [(LB, UB), (LB, 8 - LB), (LB, UB), (LB, UB), (LB, UB)]

        if self.nosnow:
            self.bounds.pop()  # remove calibration parameter M

        # create arrays to hold outputs
        self.all_pars = np.zeros((1, len(self.bounds)))
        self.kge_vals = np.zeros(1)

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

        print("\t\tFinished calibration for basin {0} which contains {1} grid cells.".format(self.basin_num, self.basin_idx[0].shape[0]))
        print("\t\tPopulation size:  {}".format(popsize))
        print("\t\tParameter values (a,b,c,d,m):  {}".format(pars))
        print("\t\tKGE:  {}".format(1 - ed))
        print("\t\tNumber of function evaluations:  {}".format(nfev))
        print("\t\tCalibration time (seconds):  {}".format(time.time() - st))

        np.save('{}/kge_result_basin_{}.npy'.format(self.out_dir, self.basin_num), self.kge_vals)
        np.save('{}/abcdm_parameters_basin_{}.npy'.format(self.out_dir, self.basin_num), self.all_pars)


def basin_runoff(pars, set_calibrate, pet, precip, tmin, n_months, runoff_spinup, obs_unit, bsn_areas, basin_idx, arr_shp, routing_func=None):
    """
    Calc

    :return:
    """

    # if calibrating against observed runoff
    if set_calibrate == 0:
        he = ABCD(pars, pet, precip, tmin, n_months, runoff_spinup, method='dist')
        he.emulate()

        if obs_unit == 'km3_per_mth':

            # convert from mm to km3
            return np.nansum(he.rsim * bsn_areas * 1e-6, 1)

        elif obs_unit == 'mm_per_mth':
            return np.nansum(he.rsim, 1)

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


def process_basin(basin_num, settings, data, pet, router_function=None):
    """
    Process single basin.
    """

    cal = Calibrate(basin_num=basin_num,
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
                    out_dir=settings.calib_out_dir)

    cal.calibrate_basin()


def calibrate_all(settings, data, pet, router_function):
    """
    Run calibration for ABCD model for all basins in parallel.
    """
    for basin_num in range(1, settings.n_basins + 1, 1):

        print("\tCalibrating Basin:  {}".format(basin_num))

        process_basin(basin_num, settings, data, pet, router_function)
