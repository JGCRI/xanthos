"""
Calculate Monthly PET using the Penman Monteith Method.

@Project: Xanthos V2.0

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute
"""

import logging
import calendar
import numpy as np


class SetData:

    def __init__(self, target_yr, data, ncells, nmonths, nlcs, st_mth_idx, ed_mth_idx, land_cover_years):

        self.c = data

        # subset data by year range
        self.c.tair = self.c.tair_load[:, st_mth_idx:ed_mth_idx]
        self.c.TMIN = self.c.TMIN_load[:, st_mth_idx:ed_mth_idx]
        self.c.rhs = self.c.rhs_load[:, st_mth_idx:ed_mth_idx]
        self.c.wind = self.c.wind_load[:, st_mth_idx:ed_mth_idx]
        self.c.rsds = self.c.rsds_load[:, st_mth_idx:ed_mth_idx]
        self.c.rlds = self.c.rlds_load[:, st_mth_idx:ed_mth_idx]

        # set previous air temp value leaving the first value at 0
        self.c.tairprev = self.c.tairprev_load[:, st_mth_idx:ed_mth_idx]

        # sort land cover years from smallest to largest
        lc_yrs = sorted(land_cover_years)

        # if year exceeds the most current land cover year, use the most current land cover year
        if target_yr >= lc_yrs[-1]:
            lc_yr = lc_yrs.index(lc_yrs[-1])

        # else use the land cover data from the year of the nearest decade
        else:
            lc_yr = lc_yrs.index([x for x in lc_yrs if x - target_yr >= -4][0])

        self.c.lct = np.swapaxes(np.tile(self.c.lct_load[:, :, lc_yr], (nmonths, 1, 1)), 0, -1)
        totpct = np.sum(self.c.lct, axis=0)
        self.c.totpct = np.where(totpct == 0, 0.01, totpct)

        self.ncells = ncells
        self.nmonths = nmonths
        self.nlcs = nlcs

        self.lcr = np.arange(self.nlcs)
        self.yrs = int(self.nmonths / 12)

        # account for leap year
        if calendar.isleap(target_yr):
            self.dys_in_mth = np.array([31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
        else:
            self.dys_in_mth = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])

        self.dz = np.tile(self.dys_in_mth, self.yrs)

        # build index of months for a year
        self.dm = np.arange(12)

        self.mth_idx = np.tile(np.arange(12), self.yrs)
        self.mth_le_5 = np.where(self.mth_idx <= 5)
        self.mth_gr_5 = np.where(self.mth_idx > 5)

        alpx = np.tile(self.c.alpha, self.yrs)
        alpd = np.tile(alpx, (self.ncells, 1, 1))
        self.alpha = np.swapaxes(alpd, 0, -2)

        # PENMAN-MONTEITH CONSTANTS
        self.lambda1 = 2.46e6   # latent heat of vaporization (J kg-1)
        self.cp = 1006          # specific heat capacity of air (J kg-1 k-1)
        self.sigma = 4.9e-3     # Stephan-Boltzmann constant (J m-2 K-4 d-1)
        self.sigma2 = 5.67e-8   # Stephan-Boltzmann constant (J m-2 K-4 s-1)
        self.gamma = 0.67       # Psychrometric constant (mbar k-1),mbar degree-1 is same meaning with mbar k-1
        self.k = -0.5           # extinction coefficient

        self.esx = np.tile(6.10588 * np.exp(17.32491 * self.c.tair / (self.c.tair + 238.102)), (self.nlcs, 1, 1))

        # vapor pressure
        self.vapx = np.multiply(self.esx, self.c.rhs / 100)

        # mbar degree-1 is same meaning with mbar k-1
        self.sx = (238.1 * 17.325 * self.esx / np.power((self.c.tair + 238.1), 2))

        self.rsdsx = np.tile(self.c.rsds, (nlcs, 1, 1))

        # net solar radiation(J m-2 mon-1)
        self.rsnx = (1 - self.alpha) * self.rsdsx * 86400 * self.dz

        self.tairx = np.tile(self.c.tair, (nlcs, 1, 1))

        # see http:#en.wikipedia.org/wiki/Wind_profile_power_law
        self.wind2x = np.tile(self.c.wind * np.power(2 / 10, 0.11), (nlcs, 1, 1))


def calc_mtmin(v, nlcs, nmonths, ncells):
    tmin = np.tile(v.c.TMIN, (nlcs, 1, 1))
    tminopen = np.swapaxes(np.tile(v.c.Tminopen, (nmonths, ncells, 1)), 0, -1)
    tminclose = np.swapaxes(np.tile(v.c.Tminclose, (nmonths, ncells, 1)), 0, -1)

    mtmin = np.zeros_like(tmin)
    mtmin[tmin >= tminopen] = 1.0
    mtmin[tmin <= tminclose] = 0.1
    xi = (tmin < tminopen) & (tmin > tminclose)
    xii = (tmin - tminclose) / (tminopen - tminclose)
    mtmin = np.where(xi, xii, mtmin)

    return mtmin


def calc_vpd(v, nmonths, ncells, vap):
    vpdopen = np.swapaxes(np.tile(v.c.VPDopen, (nmonths, ncells, 1)), 0, -1)
    vpdclose = np.swapaxes(np.tile(v.c.VPDclose, (nmonths, ncells, 1)), 0, -1)

    vpd = v.esx - vap
    mvpd = vpd.copy()
    mvpd = np.where(vpd <= vpdopen, 1.0, mvpd)
    mvpd = np.where(vpd >= vpdclose, 0.1, mvpd)
    vi = (vpd > vpdopen) & (vpd < vpdclose)
    vii = (vpdclose - vpd) / (vpdclose - vpdopen)
    mvpd = np.where(vi, vii, mvpd)

    return mvpd, vpd, vpdopen, vpdclose


def calc_rtotc(v, nmonths, ncells, vpd, vpdopen, vpdclose):
    rtotc = np.zeros_like(vpd)
    rblmax = np.tile(v.c.RBLmax, (nmonths, ncells, 1))
    rblmax = np.swapaxes(rblmax, 0, -1)
    rblmin = np.tile(v.c.RBLmin, (nmonths, ncells, 1))
    rblmin = np.swapaxes(rblmin, 0, -1)

    rtotc = np.where(vpd <= vpdopen, rblmax, rtotc)
    rtotc = np.where(vpd >= vpdclose, rblmin, rtotc)
    ri = (vpd > vpdopen) & (vpd < vpdclose)
    rii = rblmax - (rblmax - rblmin) * (vpdclose - vpd) / (vpdclose - vpdopen)
    rtotc = np.where(ri, rii, rtotc)

    return rtotc


def calc_a(v, nlcs, nmonths, ncells, dz):
    """
    Calculate available energy.

    Net outgoiong long-wave radiation(J m-2 mon-1), remove f,
    emiss here is land surface emissivity.
    """
    dm = np.tile(v.dm, (nlcs, ncells, int(nmonths / 12)))  # TODO: remove
    emiss = np.swapaxes(np.tile(v.c.emiss, (nmonths, ncells, 1)), 0, -1)
    rlds = np.tile(v.c.rlds, (nlcs, 1, 1))
    rnl = v.sigma * np.power(v.tairx + 273, 4.0) * emiss * dz - rlds * 86400 * dz
    rn = ((1 - v.alpha) * v.rsdsx) * 86400 * dz - rnl
    a = rn / (86400 * dz)

    return a


def calc_fwet(rh):
    fwet = np.where(rh < 70, 0, rh)
    fwet = np.where(rh >= 70, np.power(rh / 100, 8), fwet)
    fwet = np.where(rh >= 80, np.power(rh / 100, 10), fwet)
    fwet = np.where(rh >= 90, np.power(rh / 100, 12), fwet)
    fwet = np.where(rh >= 95, np.power(rh / 100, 16), fwet)

    return fwet


def calc_lai(v, ncells, nmonths):
    yrs = int(nmonths / 12)

    lai = np.swapaxes(np.tile(v.c.lai, (ncells, 1, yrs)), 0, -2)
    laimin = np.swapaxes(np.tile(v.c.laimin, (ncells, 1, yrs)), 0, -2)
    laimax = np.swapaxes(np.tile(v.c.laimax, (ncells, 1, yrs)), 0, -2)

    return lai, laimin, laimax


def calc_p(v, nlcs):
    elevx = np.tile(v.c.elev, (nlcs, 1, 1))
    p = 101325 * np.power((1 - 0.0065 * elevx / 288.15), 5.2558)

    return p


def calc_cc(lai, fwet, gs1, rcx, gcu):
    cc = np.where((gs1 + 1 / rcx + gcu) < 0.0001, 10000,
                  np.where(fwet == 1, 0.00001, np.where(lai < 0.0001, 0.00001, 0)))
    ccx = np.where(cc == 0, 1 / rcx * (gs1 + gcu) * lai * (1 - fwet) / (gs1 + 1 / rcx + gcu), cc)

    return ccx


def calc_rcx(v, ncells, nmonths):

    return np.swapaxes(np.swapaxes(np.tile(v.c.rc, (ncells, nmonths, 1)), 0, -1), 1, -1)


def calc_rh(v, nlcs):
    rh = np.tile(v.c.rhs, (nlcs, 1, 1))
    rh[rh > 99.9999] = 99.9

    return rh


def calc_g(v, nlcs):
    g = 1.6198 * (v.tairx - np.tile(v.c.tairprev, (nlcs, 1, 1)))
    g[:, :, 0] = 0

    return g


def calc_clx(v, ncells, nmonths):
    return np.swapaxes(np.swapaxes(np.tile(v.c.cL, (ncells, nmonths, 1)), 0, -1), 1, -1)


def et_veg(v, nlcs, nmonths, ncells):
    dz = np.tile(v.dz, (nlcs, ncells, 1))

    # for et_veg
    p = calc_p(v, nlcs)

    rcorr = p / (101300 * np.power((273.15 + v.tairx) / 293.15, 1.75))
    gcu = 0.00001 * rcorr

    # calculate mTmin
    mtmin = calc_mtmin(v, nlcs, nmonths, ncells)

    vap = np.multiply(v.esx, v.c.rhs / 100)

    # calculate mVPD
    mvpd, vpd, vpdopen, vpdclose = calc_vpd(v, nmonths, ncells, vap)

    clx = calc_clx(v, ncells, nmonths)

    gs1 = clx * mtmin * mvpd * rcorr

    rh = calc_rh(v, nlcs)

    # calculate rtotc
    rtotc = calc_rtotc(v, nmonths, ncells, vpd, vpdopen, vpdclose)

    # w m-2
    g = calc_g(v, nlcs)

    # available energy(w m-2)
    a = calc_a(v, nlcs, nmonths, ncells, dz)

    lai, laimin, laimax = calc_lai(v, ncells, nmonths)

    fc_denom = np.exp(-0.5 * laimin) - np.exp(-0.5 * laimax)
    fc_denom = np.where(fc_denom == 0.0, 1, fc_denom)

    fc = (np.exp(-0.5 * laimin) - np.exp(-0.5 * lai)) / fc_denom
    fc = np.where(fc > 1, 1, fc)

    ac = fc * a

    # consider the G
    asoil = (1 - fc) * a - g

    rtot = rtotc * rcorr
    rtot = np.where(rtot > 80, 80, rtot)

    rho = p / ((v.tairx + 273.15) * 287.058)

    rr = rho * v.cp / (4.0 * v.sigma2 * np.power((v.tairx + 273.15), 3))

    rcx = calc_rcx(v, ncells, nmonths)

    ra = rcx * rr / (rcx + rr)
    ra = np.where(ra > rtot, rtot, ra)

    fwet = calc_fwet(rh)

    cc = calc_cc(lai, fwet, gs1, rcx, gcu)

    rs = np.zeros_like(lai)
    rs = np.where(cc == 0, 100000, 1 / cc)

    rslimit = np.tile(v.c.rslimit, (ncells, nmonths, 1))
    rslimit = np.swapaxes(rslimit, 0, -1)
    rslimit = np.swapaxes(rslimit, 1, -1)

    rs = np.where(rs > rslimit, rslimit, rs)

    rrc = rho * v.cp / (4 * v.sigma2 * np.power((273.15 + v.tairx), 3))

    # account for division by 0 warning
    lai_fwet = np.where(lai * fwet == 0, 1, lai * fwet)
    rhc = np.where(lai > 0.00001, rcx / lai_fwet, rslimit)
    rvc = np.where(lai > 0.00001, rcx / lai_fwet, rslimit)

    rhc = np.where(rhc > rslimit, rslimit, rhc)
    rvc = np.where(rvc > rslimit, rslimit, rvc)

    rhrc = rhc * rrc / (rhc + rrc)
    rhrc = np.where(rhrc > rtot, rtot, rhrc)

    apres = dz * 86400 * (v.sx * ac + rho * v.cp * (v.esx - vap) * fc / rhrc) * fwet / (
        (v.sx + p * 0.01 * v.cp * rvc / (v.lambda1 * 0.622 * rhrc)) * v.lambda1)

    ewet_c = np.zeros_like(lai)
    ewet_c = np.where(rh >= 70, apres, ewet_c)

    rasoil = rtot * rr / (rtot + rr)

    ewet_soil = 86400 * dz * (v.sx * asoil + rho * v.cp * (1 - fc) * (v.esx - vap) / rasoil) * fwet / (
        (v.sx + v.gamma * rtot / rasoil) * v.lambda1)
    esoilpot = 86400 * dz * (v.sx * asoil + rho * v.cp * (1 - fc) * (v.esx - vap) / rasoil) * (1 - fwet) / (
        (v.sx + v.gamma * rtot / rasoil) * v.lambda1)

    beta = np.tile(v.c.beta, (ncells, nmonths, 1))
    beta = np.swapaxes(beta, 0, -1)
    beta = np.swapaxes(beta, 1, -1)

    esoil = ewet_soil + esoilpot * np.power((rh / 100), (v.esx - vap) / beta)

    # ET from the plants (mm mon-1)
    trans = dz * 86400 * (v.sx * ac + rho * v.cp * (v.esx - vap) * fc / ra) * (1 - fwet) / (
        (v.sx + v.gamma * (1 + rs / ra)) * v.lambda1)
    trans = np.where(fc == 0, 0, trans)

    eet = trans + ewet_c + esoil

    eet = np.where(eet < 0.0, 0.0, eet)

    return eet


def et_water(v, nmonths, emiss):
    # net outgoiong long-wave radiation(J m-2 mon-1), emiss here is land surface emissivity
    rnlx = np.tile(v.sigma * np.power(v.c.tair + 273, 4.0) * emiss * v.dz - v.c.rlds * 86400 * v.dz, (v.nlcs, 1, 1))

    # net all-wave radiation(J m-2 mon-1)
    rnx = ((1 - v.alpha) * v.rsdsx) * 86400 * v.dz - rnlx
    rnx[rnx < 0] = 0.0

    qtx = np.zeros(shape=(v.nlcs, v.ncells, nmonths))
    qtx[:, :, v.mth_le_5] = 0.5 * v.rsnx[:, :, v.mth_le_5] - 0.8 * rnlx[:, :, v.mth_le_5]
    qtx[:, :, v.mth_gr_5] = 0.5 * v.rsnx[:, :, v.mth_gr_5] - 1.3 * rnlx[:, :, v.mth_gr_5]

    # (w m-2)
    ax = (rnx - qtx) / (86400 * v.dz)
    ax[ax < 0] = 0

    rn2x = (rnx / (86400 * v.dz))

    ewetx = rn2x * v.dz * 0.6 / 2845
    ewety = v.dz * 86400 * (v.sx * ax + v.gamma * 6.43 * (0.5 + 0.54 * v.wind2x) * (v.esx - v.vapx)) / (
        (v.sx + v.gamma) * v.lambda1)
    ewet = np.where(v.tairx < -1, ewetx, ewety)
    ewet[ewet < 0.0] = 0.0

    return ewet[0, :, :]


def et_snow(v, emiss):
    # net outgoiong long-wave radiation(J m-2 mon-1), emiss here is land surface emissivity
    rnlx = np.tile(v.sigma * np.power(v.c.tair + 273, 4.0) * emiss * v.dz - v.c.rlds * 86400 * v.dz, (v.nlcs, 1, 1))

    # net all-wave radiation(J m-2 mon-1)
    rnx = ((1 - v.alpha) * v.rsdsx) * 86400 * v.dz - rnlx
    rnx[rnx < 0] = 0.0

    rn2x = (rnx / (86400 * v.dz))

    ewet = rn2x * v.dz * 0.6 / 2845
    ewet[ewet < 0.0] = 0.0

    return ewet[6, :, :]


def set_years(y, end_yr, process_nyrs):
    if y < end_yr:
        st_yr = y
    else:
        st_yr = end_yr

    if (st_yr + process_nyrs) > end_yr:
        ed_yr = end_yr
    else:
        ed_yr = y + process_nyrs

    return st_yr, ed_yr


def run_pmpet(data, ncells, nlcs, start_yr, end_yr, water_idx, snow_idx, land_cover_years, process_nyrs=1):
    """
    Run Penman-Monteith PET.

    :param data:                Object containing temperature, humidity, wind, radiation, and land cover data:
                                    - surface air temperature in degrees Celsius, 2D [ncells, nmonths]
                                    - minimum surface air temperature in degrees Celsus, 2D [ncells, nmonths]
                                    - relative humidity (rhs) in percent
                                    - wind speed in meters/second
                                    - surface downwelling shortwave radiation (rsds) in W m-2
                                    - surface longwave radiation (rlds) in W m-2
                                    - land cover data [n_cells, n_land_classes, n_yrs]
    :param ncells:              Integer number of grid cells in the data
    :param nlcs:                Number of land classes in the input data
    :param start_yr:            Start year of the project in YYYY format
    :param end_yr:              End year of the project in YYYY format
    :param water_idx:           index of water in the land cover data
    :param snow_idx:            index of snow in the land cover data
    :param land_cover_years     list of integer years in YYYY format that are contained in the land cover data
    :param process_nyrs:                    int. How many years to process at once.  Default is 1 due to local RAM c
                                            RAM constraints.  This can be increased for faster run times on systems with
                                            larger RAM available.

    :return:                    PET array [ncells, nmonths]

    """

    # calculate total project months
    tot_yrs = end_yr - start_yr
    if tot_yrs == 0:
        tot_nmonths = 12
    else:
        tot_nmonths = tot_yrs * 12 + 12

    # set output array
    out_arr = np.zeros(shape=(ncells, tot_nmonths))

    st_mth_idx = 0
    ed_mth_idx = 0
    for idx, y in enumerate(range(start_yr, end_yr + 1, process_nyrs)):

        # get start year and end year for time slice
        st_yr, ed_yr = set_years(y, end_yr, process_nyrs)

        logging.info("\t\tProcessing Years:  {} to {}".format(st_yr, ed_yr))

        # set the number of months for the time slice
        yx = ed_yr - st_yr
        if yx == 0:
            nmonths = 12
        else:
            nmonths = yx * 12

        # set end month index
        ed_mth_idx += nmonths

        # set data based on time slice
        v = SetData(y, data, ncells, nmonths, nlcs, st_mth_idx, ed_mth_idx, land_cover_years)

        # array to hold outputs for the time slice
        arr = np.zeros(shape=(nlcs, ncells, nmonths))

        # calculate vegetation
        veg = et_veg(v, nlcs, nmonths, ncells)
        arr[:, :, :] = veg

        # calculate water
        wat = et_water(v, nmonths, emiss=0.98)
        arr[water_idx, :, :] = wat

        # calculate snow
        snow = et_snow(v, emiss=0.85)
        arr[snow_idx, :, :] = snow

        # apply the land cover data
        arr *= v.c.lct

        # evaluate by total percent of land cover
        fin = np.sum(arr, axis=0) / v.c.totpct

        out_arr[:, st_mth_idx:ed_mth_idx] = fin

        # advance start month index for whole time period
        st_mth_idx += nmonths

    return out_arr
