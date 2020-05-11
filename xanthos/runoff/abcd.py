"""
ABCD runoff model.

@author   Chris R. Vernon, Caleb J. Braun
@email:   chris.vernon@pnnl.gov
@Project: Xanthos 2.0

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2018, Battelle Memorial Institute
"""

import logging
import numpy as np
from joblib import Parallel, delayed


class ABCD:
    """
    A hydrology emulator.

    | Reference:
    |
    | Liu, Y., Hejazi, M.A., Li, H., Zhang, X., (2017), A Hydrological Emulator for Global
    |    Applications, Geoscientific Model Development Discussions, DOI: 10.5194/gmd-2017-113
    |
    | Martinez, G. F., & Gupta, H. V. (2010). Toward improved identification
    |    of hydrological models: A diagnostic evaluation of the 'abcd' monthly
    |    water balance model for the conterminous United States. Water Resources Research, 46(8).

    @:param prm     Object containing calibrated data
    @:param hist    Object containing Watch data
                    hist.
    @:return a      Float from 0-1
    @:return b      Float from 100-1000 mm
    @:return c      Float from 0-1
    @:return d      Float from 0-1
    @:return m      Float from 0-1
    """

    def __init__(self, pars, pet, precip, tmin, basin_ids, process_steps, spinup_steps, method='dist'):

        # are we running with the snow component?
        self.nosnow = tmin is None

        # extract the a, b, c, d, and optionally m parameters
        self.a = pars[:, 0]
        self.b = pars[:, 1] * 1000
        self.c = pars[:, 2]
        self.d = pars[:, 3]
        self.m = pars[:, 4] if not self.nosnow else 0

        # cache parameter values used in abcd_dist calculations
        self.a_times2 = self.a * 2
        self.b_over_a = self.b / self.a
        self.d_plus_1 = self.d + 1

        self.basin_ids = basin_ids

        # set processing method and steps
        self.method = method
        self.steps = process_steps
        self.spinup_steps = spinup_steps

        # PET and precipitation as input, transpose to row based
        self.pet = pet.T[0:self.steps, :]
        self.precip = precip.T[0:self.steps, :]

        # initialize arrays for spinup values
        self.pet0 = self.pet[0:self.spinup_steps, :].copy()
        self.precip0 = self.precip[0:self.spinup_steps, :].copy()

        # temperature minimum as input, if provided
        if self.nosnow:
            self.tmin = None
            self.tmin0 = None
        else:
            self.tmin = tmin.T[0:self.steps, :]
            self.tmin0 = self.tmin[0:self.spinup_steps, :].copy()

        # values [initial runoff, soil moisture storage, groundwater storage
        self.inv = np.array([20, 100, 500])
        self.soil_water_storage0 = self.inv[1]
        self.groundwater_storage0 = self.inv[2]

        # populate param arrays
        self.snm = None
        self.actual_et = None
        self.groundwater_storage = None
        self.soil_water_storage = None
        self.water_avail = None
        self.et_op = None
        self.snowpack = None
        self.rsim = None
        self.rain = None
        self.snow = None

        self.SN0 = 0
        self.TRAIN = 2.5
        self.TSNOW = 0.6

    def set_actual_et(self, p, frac=0.6):
        """
        Set the initial value of actual evapotranspiration (ET) at 60% of precipitation.

        :return: array
        """
        arr = np.zeros_like(p)

        arr[0, :] = p[0, :] * frac

        self.actual_et = arr

    def set_snowpack(self, p):
        """
        Set the initial value of accumulated snow water at 10% of precipitation.

        :return: array
        """
        arr = np.zeros_like(p)

        arr[0, :] = p[0, :] / 10

        self.snowpack = arr

    def set_rsim(self, p, v):
        """
        Set the initial streamflow.

        :param p:   precipitation
        :param v:   initial runoff

        :return: array
        """
        arr = np.zeros_like(p)

        arr[0, :] = v

        self.rsim = arr

    def set_rain_and_snow(self, p, tmin):
        """Assign rain and snow arrays."""
        # we only need rain array if running without snow
        if self.nosnow:
            self.rain = p
            return

        # construct snow and rain arrays like precip shape
        self.rain = np.zeros_like(p)
        self.snow = np.zeros_like(p)

        # get the indices of each value meeting criteria
        allrain = tmin > self.TRAIN
        rainorsnow = (tmin <= self.TRAIN) & (tmin >= self.TSNOW)
        allsnow = tmin < self.TSNOW

        # populate the snow array
        self.snow[rainorsnow] = p[rainorsnow] * (self.TRAIN - tmin[rainorsnow]) / (self.TRAIN - self.TSNOW)

        if allrain.any():
            self.rain[allrain] = p[allrain]
            self.snow[allrain] = 0

        if rainorsnow.any():
            self.rain[rainorsnow] = p[rainorsnow] - self.snow[rainorsnow]

        if allsnow.any():
            self.rain[allsnow] = 0
            self.snow[allsnow] = p[allsnow]

    def abcd_dist(self, i, pet, tmin):
        """
        Run the ABCD model calculations.

        @:param i       Current month
        @:param pet     Potential Evapotranspiration
        @:param tmin    Monthly minimum temperature (if running with snow)
        """
        if not self.nosnow:
            if i == 0:
                self.snowpack[i, :] = self.SN0 + self.snow[i, :]
            else:
                self.snowpack[i, :] = self.snowpack[i - 1, :] + self.snow[i, :]

            # select only snow, intermediate, or only rain for each case
            allrain = tmin[i, :] > self.TRAIN
            rainorsnow = (tmin[i, :] <= self.TRAIN) & (tmin[i, :] >= self.TSNOW)
            allsnow = tmin[i, :] < self.TSNOW

            # estimate snowmelt (SNM)
            self.snm[i, allrain] = self.snowpack[i, allrain] * self.m[allrain]
            self.snm[i, rainorsnow] = (self.snowpack[i, rainorsnow] * self.m[rainorsnow]) * \
                ((self.TRAIN - tmin[i, rainorsnow]) / (self.TRAIN - self.TSNOW))
            self.snm[i, allsnow] = 0

            # accumulated snow water equivalent
            self.snowpack[i, :] -= self.snm[i, :]

        # get available water
        if i == 0:
            self.water_avail[i, :] = self.rain[i, :] + self.soil_water_storage0
        else:
            self.water_avail[i, :] = self.rain[i, :] + self.soil_water_storage[i - 1, :] + self.snm[i, :]

        # ET opportunity
        rpt = (self.water_avail[i, :] + self.b)
        rpt_over_pt2 = rpt / self.a_times2
        self.et_op[i, :] = rpt_over_pt2 - np.sqrt(np.square(rpt_over_pt2) - (self.water_avail[i, :] * self.b_over_a))

        # soil water storage
        self.soil_water_storage[i, :] = self.et_op[i, :] * np.exp(-pet[i, :].real / self.b)

        # get the difference between available water and ET opportunity
        awet = (self.water_avail[i, :] - self.et_op[i, :])
        c_x_awet = self.c * awet

        # groundwater storage
        if i == 0:
            self.groundwater_storage[i, :] = (self.groundwater_storage0 + c_x_awet) / self.d_plus_1
        else:
            self.groundwater_storage[i, :] = (self.groundwater_storage[i - 1, :] + c_x_awet) / self.d_plus_1

        # populate arrays
        self.actual_et[i, :] = self.et_op[i, :] - self.soil_water_storage[i, :]
        self.actual_et[i, :] = np.maximum(0, self.actual_et[i, :])
        self.actual_et[i, :] = np.minimum(pet[i, :].real, self.actual_et[i, :])
        self.soil_water_storage[i, :] = self.et_op[i, :] - self.actual_et[i, :]
        self.rsim[i, :] = (awet - c_x_awet) + self.d * self.groundwater_storage[i, :]

    def init_arrays(self, p, tmin):
        """Initialize arrays based on spin-up or simulation run status."""
        # construct simulated runoff, actual evapotranspiration, and snowpack arrays
        self.set_rsim(p, self.inv[0])
        self.set_actual_et(p)
        self.set_snowpack(p)

        self.snm = np.zeros_like(p)
        self.groundwater_storage = np.zeros_like(p)
        self.soil_water_storage = np.zeros_like(p)
        self.water_avail = np.zeros_like(p)
        self.et_op = np.zeros_like(p)

        # partition snow and rain
        self.set_rain_and_snow(p, tmin)

    def set_vals(self):
        """
        Set and reset initial values.

        Reset initial values for runoff [0], soil moisture [1], and groundwater
        storage [2] based on spin-up as the average of the last three Decembers.
        """
        try:
            # get Decembers from end of array (-1=last december, -13=two years ago, -25=three years ago)
            dec_idx = [-1, -13, -25]

            rsim_rollover = self.rsim[dec_idx, :]
            sm_rollover = self.soil_water_storage[dec_idx, :]
            gs_rollover = self.groundwater_storage[dec_idx, :]

        except IndexError:
            logging.exception(
                'Spin-up steps must produce at least 10 years spin-up. Your '
                'spin-up only consist of {} months. Please reconfigure and try '
                'again.'.format(self.spinup_steps))
            raise

        # initial values arrays (1d, where length is the sum of number of cells for all basins running)
        ro = np.empty(self.basin_ids.shape)
        sm = np.empty(self.basin_ids.shape)
        gs = np.empty(self.basin_ids.shape)

        # initial values are set by taking the mean value by basin, so here we have to go basin-by-basin
        for i in np.unique(self.basin_ids):
            b_idx = (i == self.basin_ids)
            ro[b_idx] = np.mean(np.nanmean(rsim_rollover[:, b_idx], axis=1))
            sm[b_idx] = np.mean(np.nanmean(sm_rollover[:, b_idx], axis=1))
            gs[b_idx] = np.mean(np.nanmean(gs_rollover[:, b_idx], axis=1))

        self.inv = np.array([ro, sm, gs])
        self.soil_water_storage0 = self.inv[1]
        self.groundwater_storage0 = self.inv[2]

    def spinup(self):
        """Run spin-up using initial values."""
        # initialize arrays for spinup
        self.init_arrays(self.precip0, self.tmin0)

        # run spin-up with initial settings and calibrated a, b, c, d, m params
        for i in range(0, self.spinup_steps, 1):
            self.abcd_dist(i, self.pet0, self.tmin0)

        # reset initial runoff, soil moisture, and groundwater storage values to the average of the last 3 Decembers
        self.set_vals()

    def simulate(self):
        """Run simulation using spin-up values."""
        # initialize arrays for simulation from spin-up
        self.init_arrays(self.precip, self.tmin)

        # process with first pass parameters
        for i in range(0, self.steps, 1):
            self.abcd_dist(i, self.pet, self.tmin)

    def emulate(self):
        """Run hydrologic emulator."""
        # run spin-up
        self.spinup()

        # run simulation
        self.simulate()


def _run_basins(basin_nums, pars_abcdm, basin_ids, pet, precip, tmin, n_months, spinup_steps, method='dist'):
    """
    Run the ABCD model for each basin.

    :param basin_nums:      The numbers of the target basins (1d NumPy array)
    :param basin_ids:       Basin ID Map: 67420 x 1, 235 Basins
    :param n_months:        The number of months to process
    :param spinup_steps:    How many times to tile the historic months by
    :param method:          Either 'dist' for distributed, or 'lump' for lumped processing
    :return                 A NumPy array
    """
    # get the grid cell indices of the target basins (integer 1d array)
    basin_indices = np.where(np.isin(basin_ids, basin_nums))

    # pass basin ids to model for setting basin-specific initial values
    target_basins_ids = basin_ids[basin_indices]

    # import ABCD parameters for the target basins (constant for all months)
    pars_by_cell = pars_abcdm[basin_ids - 1]
    pars = pars_by_cell[basin_indices]

    # extract data for the selected basins only (over all months)
    _pet = pet[basin_indices]
    _precip = precip[basin_indices]

    # tmin is optional, but the snow component of the model will not be used without it
    if tmin is None:
        _tmin = tmin
    else:
        _tmin = tmin[basin_indices]

    # instantiate the model
    he = ABCD(pars, _pet, _precip, _tmin, target_basins_ids, n_months, spinup_steps, method=method)

    # run it
    he.emulate()

    # stack outputs
    vals = np.hstack([he.pet.T, he.actual_et.T, he.rsim.T, he.soil_water_storage.T])

    return vals


def abcd_parallel(n_basins, pars, basin_ids, pet, precip, tmin, n_months, spinup_steps, jobs=-1):
    """
    Run ABCD model on basins in parallel.

    This model can run any number of basins at a time.  Splitting them into
    chunks and running them in parallel greatly speeds things up.

    :param pars:                Array of abcdm parameters by grid cell
    :param n_basins:            How many basins to run
    :return:                    A list of NumPy arrays
    """
    # rough optimization for dividing basins across threads
    if jobs < 1:
        n_chunks = 8
    else:
        n_chunks = jobs * 2

    # Split the range of basins into separate chunks for each thread to process
    min_basin = min(basin_ids)
    basin_ranges = np.array_split(np.arange(min_basin, min_basin + n_basins), n_chunks)

    logging.info("\t\tProcessing spin-up and simulation for basins {}...{}".format(min_basin, n_basins))

    rslts = Parallel(n_jobs=jobs, backend="threading")(delayed(_run_basins)
                                                       (i, pars, basin_ids, pet, precip,
                                                        tmin, n_months, spinup_steps) for i in basin_ranges)

    out = np.empty((len(basin_ids), n_months * 4))  # 4 output variables for every month

    # combine parallel results back in the correct order
    for i, br in enumerate(basin_ranges):
        basin_indices = np.where(np.isin(basin_ids, br))
        out[basin_indices, :] = rslts[i]

    return out


def abcd_execute(n_basins, basin_ids, pet, precip, tmin, calib_file, n_months, spinup_steps, jobs):
    """
    Run the ABCD model.

    :param n_basins:          How many basin to run
    :param basin_ids:         Basin ID Map: 67420 x 1, 235 Basins
    :param pet:               Potential Evapotranspiration for each cell
    :param precip:            Precipitation for each cell
    :param tmin:              Monthly average minimum temperature (optional)
    :param calib_file:        Path to .npy file containing calibrated abcdm parameters
    :param n_months:          The number of months to process
    :param spinup_steps:      How many times to tile the historic months by
    :param jobs:              The number of jobs to use when running basins parallel
                              (-2, all but one core; -1, all cores; 8, 8 cores)

    :return         A NumPy arrays for coordinates (long, lat) and simulated runoff,
                    PET, AET and soil moisture with the shape (grid cells, value per month).
    """
    # read parameters from calibration
    prm = np.load(calib_file)

    # run all basins at once in parallel
    rslts = abcd_parallel(n_basins=n_basins, pars=prm, basin_ids=basin_ids, pet=pet, precip=precip,
                          tmin=tmin, n_months=n_months, spinup_steps=spinup_steps, jobs=jobs)

    # build array to pass to router
    _pet, _aet, _q, _sav = np.split(rslts, 4, axis=1)

    return _pet, _aet, _q, _sav
