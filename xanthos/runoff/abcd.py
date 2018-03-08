"""
ABCD runoff model.

@author   Chris R. Vernon
@email:   chris.vernon@pnnl.gov
@Project: Xanthos 2.0

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute
"""

import numpy as np
import os
import sys
from joblib import Parallel, delayed


class ABCD:
    """
    | Hydrology emulator
    |
    | Reference:
    |
    | Liu, Y., Hejazi, M.A., Li, H., Zhang, X., (2017), A Hydrological Emulator for Global
    |    Applications, Geoscientific Model Development Discussions, DOI: 10.5194/gmd-2017-113
    |
    | Martinez, G. F., & Gupta, H. V. (2010). Toward improved identification
    |    of hydrological models: A diagnostic evaluation of the 'abcs' monthly
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

    def __init__(self, pars, pet, precip, tmin, process_steps, spinup_factor, method='dist'):

        # load data
        # self.prm = prm
        # self.hist = hist

        # set processing method and steps
        self.method = method
        self.steps = process_steps
        self.spinup_factor = spinup_factor
        self.spinup_steps = process_steps * spinup_factor

        # assign attributes
        self.a = pars[0]
        self.b = pars[1] * 1000
        self.c = pars[2]
        self.d = pars[3]
        self.m = pars[4]

        # values coordinate to [a, b, c, d, m]
        # self.lb = np.array([0.001, 0.1, 0, 0, 0])
        # self.ub = np.array([1, 4, 1, 1, 1])
        # self.nn = len(self.lb)

        # values [initial runoff, soil moisture storage, groundwater storage
        self.inv = np.array([20, 100, 500])
        self.s0 = self.inv[1]
        self.g0 = self.inv[2]

        if method == 'lump':

            # get mean array of each item
            self.pet = np.nanmean(pet[:, 0:self.steps], axis=0)
            self.p = np.nanmean(precip[:, 0:self.steps], axis=0)
            self.tmin = np.nanmean(tmin[:, 0:self.steps], axis=0)

            # self.robs = np.nanmean(hist.robs[:, 0:self.steps], axis=0)
            # self.maxstor = np.nanmean(hist.maxstor, axis=0)
            # self.bfi = np.nanmean(hist.bfi[:, 0:self.steps], axis=0)

        else:

            # PET as input
            self.pet = pet[:, 0:self.steps]

            # precipitation as input
            self.p = precip[:, 0:self.steps]

            # temperature minimum as input
            self.tmin = tmin[:, 0:self.steps]

            # self.robs = hist.robs[:, 0:self.steps]
            # self.maxstor = hist.maxstor
            # self.bfi = hist.bfi

        # set start points
        self.pet0 = np.tile(self.pet, self.spinup_factor)
        self.p0 = np.tile(self.p, self.spinup_factor)
        self.tmin0 = np.tile(self.tmin, self.spinup_factor)

        # self.robs0 = np.tile(self.robs, self.spinup_factor)
        # self.maxstor0 = np.tile(self.maxstor, self.spinup_factor)

        # populate param arrays
        self.snm = None
        self.ea = None
        self.re = None
        self.dr = None
        self.base = None
        self.g = None
        self.s = None
        self.w = None
        self.y = None
        self.xs = None
        self.rsim = None
        self.rain = None
        self.snow = None
        self.sn0 = 0
        self.train = 2.5
        self.tsnow = 0.6

    def set_ea(self, p, frac=0.6):
        """
        Set the initial value of actual evapotranspiration (ET) at 60% of precipitation.

        :return: array
        """
        arr = np.zeros_like(p)

        if self.method == 'lump':
            arr[0] = p[0] * frac

        else:
            arr[:, 0] = p[:, 0] * frac

        return arr

    def set_xs(self, p):
        """
        Set the initial value of accumulated snow water at 10% of precipitation.

        :return: array
        """
        arr = np.zeros_like(p)

        if self.method == 'lump':
            arr[0] = p[0] / 10

        else:
            arr[:, 0] = p[:, 0] / 10

        return arr

    def set_rsim(self, p, v):
        """
        Set the initial streamflow.

        :return: array
        """
        arr = np.zeros_like(p)

        if self.method == 'lump':
            arr[0] = v

        else:
            arr[:, 0] = v

        return arr

    def get_rs(self, p, tmin):
        """
        Assign rain and snow arrays.
        """

        # construct snow and rain arrays like precip shape
        self.rain = np.zeros_like(p)
        self.snow = np.zeros_like(p)

        # get the index of each value meeting criteria
        allrain = np.nonzero(tmin > self.train)
        rainorsnow = np.nonzero((tmin <= self.train) & (tmin >= self.tsnow))
        allsnow = np.nonzero(tmin < self.tsnow)

        # populate the snow array
        self.snow[rainorsnow] = p[rainorsnow] * ((self.train - tmin[rainorsnow])) / (self.train - self.tsnow)

        if allrain:
            self.rain[allrain] = p[allrain]
            self.snow[allrain] = 0

        if rainorsnow:
            self.rain[rainorsnow] = p[rainorsnow] - self.snow[rainorsnow]

        if allsnow:
            self.rain[allsnow] = 0
            self.snow[allsnow] = p[allsnow]

    def abcd_lump(self, i, pet, tmin):

        if i == 0:
            self.xs[i] = self.sn0 + self.snow[i]
        else:
            self.xs[i] = self.xs[i-1] + self.snow[i]

        # select only snow, intermediate, or only rain for each case
        allrain = np.nonzero(tmin[i] > self.train)
        rainorsnow = np.nonzero((tmin[i] <= self.train) & (tmin[i] >= self.tsnow))
        allsnow = np.nonzero(tmin[i] < self.tsnow)

        # estimate snowmelt (SNM)
        self.snm[allrain] = self.xs[allrain] * self.m
        self.snm[rainorsnow] = (self.xs[rainorsnow] * self.m) * ((self.train - tmin[rainorsnow]) / (self.train - self.tsnow))
        self.snm[allsnow] = 0

        # accumulated snow water equivalent
        self.xs[i] = self.xs[i] - self.snm[i]

        # get available water
        if i == 0:
            self.w[i] = self.rain[i] + self.s0
        else:
            self.w[i] = self.rain[i] + self.s[i-1] + self.snm[i]

        # ET opportunity
        self.y[i] = (self.w[i] + self.b) / (2 * self.a) - np.power(np.power((self.w[i] + self.b) / (2 * self.a), 2) - (self.w[i] * self.b / self.a), 0.5)

        # soil water storage
        self.s[i] = self.y[i] * np.exp(-pet[i] / self.b)

        # get the difference between available water and ET opportunity
        awet = (self.w[i] - self.y[i])

        # groundwater storage
        if i == 0:
            self.g[i] = (self.g0 + self.c * awet) / (1 + self.d)
        else:
            self.g[i] = (self.g[i-1] + self.c * awet) / (1 + self.d)

        # populate arrays
        self.ea[i] = self.y[i] - self.s[i]
        self.ea[i] = np.maximum(0, self.ea[i])
        self.ea[i] = np.minimum(pet[i], self.ea[i])
        self.s[i] = self.y[i] - self.ea[i]
        self.re[i] = self.c * awet
        self.dr[i] = (1 - self.c) * awet
        self.rsim[i] = (1 - self.c) * awet + self.d * self.g[i]
        self.base[i] = self.d * self.g[i]

    def abcd_dist(self, i, pet, tmin):

        if i == 0:
            self.xs[:, i] = self.sn0 + self.snow[:, i]
        else:
            self.xs[:, i] = self.xs[:, i-1] + self.snow[:, i]

        # select only snow, intermediate, or only rain for each case
        allrain = np.nonzero(tmin[:, i] > self.train)
        rainorsnow = np.nonzero((tmin[:, i] <= self.train) & (tmin[:, i] >= self.tsnow))
        allsnow = np.nonzero(tmin[:, i] < self.tsnow)

        # estimate snowmelt (SNM)
        self.snm[allrain, i] = self.xs[allrain, i] * self.m
        self.snm[rainorsnow, i] = (self.xs[rainorsnow, i] * self.m) * ((self.train - tmin[rainorsnow, i]) / (self.train - self.tsnow))
        self.snm[allsnow, i] = 0

        # accumulated snow water equivalent
        self.xs[:, i] = self.xs[:, i] - self.snm[:, i]

        # get available water
        if i == 0:
            self.w[:, i] = self.rain[:, i] + self.s0
        else:
            self.w[:, i] = self.rain[:, i] + self.s[:, i-1] + self.snm[:, i]

        # ET opportunity
        self.y[:, i] = (self.w[:, i] + self.b) / (2 * self.a) - np.power(np.power((self.w[:, i] + self.b) / (2 * self.a), 2) - (self.w[:, i] * self.b / self.a), 0.5)

        # soil water storage
        self.s[:, i] = self.y[:, i] * np.exp(-pet[:, i].real / self.b)

        # get the difference between available water and ET opportunity
        awet = (self.w[:, i] - self.y[:, i])

        # groundwater storage
        if i == 0:
            self.g[:, i] = (self.g0 + self.c * awet) / (1 + self.d)
        else:
            self.g[:, i] = (self.g[:, i-1] + self.c * awet) / (1 + self.d)

        # populate arrays
        self.ea[:, i] = self.y[:, i] - self.s[:, i]
        self.ea[:, i] = np.maximum(0, self.ea[:, i])
        self.ea[:, i] = np.minimum(pet[:, i].real, self.ea[:, i])
        self.s[:, i] = self.y[:, i] - self.ea[:, i]
        self.re[:, i] = self.c * awet
        self.dr[:, i] = (1 - self.c) * awet
        self.rsim[:, i] = (1 - self.c) * awet + self.d * self.g[:, i]
        self.base[:, i] = self.d * self.g[:, i]

    def init_arrays(self, p, tmin):
        """
        Initialize arrays based on spin-up or simulation run status.
        """
        # construct simulated runoff array
        self.rsim = self.set_rsim(p, self.inv[0])

        self.snm = np.zeros_like(p)
        self.ea = self.set_ea(p)
        self.re = np.zeros_like(p)
        self.dr = np.zeros_like(p)
        self.base = np.zeros_like(p)
        self.g = np.zeros_like(p)
        self.s = np.zeros_like(p)
        self.w = np.zeros_like(p)
        self.y = np.zeros_like(p)
        self.xs = self.set_xs(p)

        # partition snow and rain
        self.get_rs(p, tmin)

    def set_vals(self):
        """
        Reset initial values for runoff [0], soil moisture [1], and groundwater storage [2] based on spin-up as the
        average of the last three Decembers.
        """
        try:
            # get Decembers from end of array (-1=last december, -13=two years ago, -25=three years ago)
            dec_idx = [-1, -13, -25]

        except IndexError:
            print('Spin-up factor must produce at least 10 years spin-up.')
            print('Your factor ({0}) when multiplied by the total years ({1}) yeids a spin-up of {2} years.'.format(self.spinup_factor, self.steps/12, self.spinup_steps/12))
            print('Please reconfigure and try again.')
            sys.exit()

        if self.method == 'lump':
            rsim_rollover = self.rsim[dec_idx]
            sm_rollover = self.s[dec_idx]
            gs_rollover = self.g[dec_idx]
        else:
            rsim_rollover = self.rsim[:, dec_idx]
            sm_rollover = self.s[:, dec_idx]
            gs_rollover = self.g[:, dec_idx]

        ro = np.mean(np.nanmean(rsim_rollover, axis=0))
        sm = np.mean(np.nanmean(sm_rollover, axis=0))
        gs = np.mean(np.nanmean(gs_rollover, axis=0))

        self.inv = np.array([ro, sm, gs])
        self.s0 = self.inv[1]
        self.g0 = self.inv[2]

    def spinup(self):
        """
        Run spin-up using initial values.
        """
        # initialize arrays for spinup
        self.init_arrays(self.p0, self.tmin0)

        # run spin-up with initial settings and calibrated a, b, c, d, m params
        for i in range(0, self.spinup_steps, 1):

            if self.method == 'lump':
                self.abcd_lump(i, self.pet0, self.tmin0)

            else:
                self.abcd_dist(i, self.pet0, self.tmin0)

        # reset initial runoff, soil moisture, and groundwater storage values to the average of the last three Decembers
        self.set_vals()

    def simulate(self):
        """
        Run simulation using spin-up values.
        """
        # initialize arrays for simulation from  spin-up
        self.init_arrays(self.p, self.tmin)

        # process with first pass parameters
        for i in range(0, self.steps, 1):

            if self.method == 'lump':
                self.abcd_lump(i, self.pet, self.tmin)

            else:
                self.abcd_dist(i, self.pet, self.tmin)

    def emulate(self):
        """
        Run hydrologic emulator.
        """

        # run spin-up
        self.spinup()

        # run simulation
        self.simulate()


def _run_basin(basin_num, pars_abcdm, basin_ids, pet, precip, tmin, n_months, spinup_factor=2, method='dist'):
    """
    Run the ABCD model for each basin.

    :param basin_num:       The number of the target basin
    :param n_months:        The number of months to process
    :param spinup_factor    How many times to tile the historic months by
    :param method:          Either 'dist' for distributed, or 'lump' for lumped processing
    :param verbose:         True if all ABCD parameters are to be exported, False if only coords and rsim (Default)
    :return                 A NumPy array
    """

    # import ABCD parameters for the target basin
    pars = pars_abcdm[basin_num-1]

    # get index where equals basin
    basin_idx = np.where(basin_ids == basin_num)

    # extract data where equals basin
    _pet = pet[basin_idx]
    _precip = precip[basin_idx]
    _tmin = tmin[basin_idx]

    # instantiate the model
    he = ABCD(pars, _pet, _precip, _tmin, n_months, spinup_factor, method=method)

    # run it
    he.emulate()

    # stack outputs
    vals = np.hstack([he.pet, he.ea, he.rsim, he.s])

    return vals

def abcd_parallel(num_basins, pars, basin_ids, pet, precip, tmin, n_months, spinup_factor=1, jobs=-1):
    """
    This model is made to run a basin at a time.  Running them in parallel
    greatly speeds things up. The user can choose to output just the coordinates
    and simulated runoff, or the whole suite of outputs from the ABCD model.

    :param num_basins:          How many basin to run
    :return:                    A list of NumPy arrays
    """
    rslts = Parallel(n_jobs=jobs)(delayed(_run_basin)(i, pars, basin_ids, pet, precip, tmin, n_months, spinup_factor) for i in range(1, num_basins + 1, 1))
    return(rslts)

def abcd_outputs(rslts, n_months, basin_ids, ncells):
    """
    Load each saved array of outputs from the ABCD runs and stack them.

    :param rslts:   Python list containing ABCD outputs in order.
    :return         A NumPy array for coordinates (long, lat) and simulated runoff
                    with the shape (grid cells, value per month).
    """
    _aet = np.zeros(shape=(ncells, n_months))
    _q = np.zeros(shape=(ncells, n_months))
    _sav = np.zeros(shape=(ncells, n_months))
    _pet = np.zeros(shape=(ncells, n_months))

    # for each file, load it and then stack it
    for idx, arr in enumerate(rslts):
        # get basin idx...
        basin_idx = np.where(basin_ids == idx + 1)

        start = 0
        end = start + n_months

        # slice out pet
        _pet[basin_idx] = arr[:, start:end]

        start = end
        end = start + n_months

        # slice out aet
        _aet[basin_idx] = arr[:, start:end]

        start = end
        end = start + n_months

        # slice out q
        _q[basin_idx] = arr[:, start:end]

        start = end
        end = start + n_months

        # slice out soil moisture storage
        _sav[basin_idx] = arr[:, start:end]

    return _pet, _aet, _q, _sav


def abcd_execute(n_basins, basin_ids, pet, precip, tmin, calib_file, n_months, spinup_factor, jobs):
    """
    Run the ABCD model.
    :return:
    """
    # read parameters from calibration
    prm = np.load(calib_file)

    # run all basins at once in parallel
    all_bsns = abcd_parallel(num_basins=n_basins, pars=prm, basin_ids=basin_ids,
               pet=pet, precip=precip, tmin=tmin, n_months=n_months,
               spinup_factor=spinup_factor, jobs=jobs)

    # build array to pass to router
    _pet, _aet, _q, _sav = abcd_outputs(all_bsns, n_months=n_months, basin_ids=basin_ids, ncells=67420)

    return _pet, _aet, _q, _sav
