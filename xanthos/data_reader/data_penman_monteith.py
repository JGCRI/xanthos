import numpy as np
from xanthos.data_reader.data_utils import DataUtils


class DataPenmanMonteith(DataUtils):

    def __init__(self, config_obj=None, params_file=None, albedo_file=None, lai_file=None, laimax_file=None,
                 laimin_file=None, tas_file=None, tasmin_file=None, rhs_file=None, wind_file=None, rsds_file=None,
                 rlds_file=None, lulc_file=None, elev_file=None, start_yr=None, through_yr=None):

        # enbale passing of args for config
        if config_obj is None:
            self.pm_params = params_file
            self.pm_alpha = albedo_file
            self.pm_lai = lai_file
            self.pm_laimax = laimax_file
            self.pm_laimin = laimin_file
            self.pm_tas = tas_file
            self.pm_tmin = tasmin_file
            self.pm_rhs = rhs_file
            self.pm_wind = wind_file
            self.pm_rsds = rsds_file
            self.pm_rlds = rlds_file
            self.pm_lct = lulc_file
            self.pm_elev = elev_file
            nmonths = (through_yr - start_yr + 1) * 12

        else:
            self.pm_params = config_obj.pm_params
            self.pm_alpha = config_obj.pm_alpha
            self.pm_lai = config_obj.pm_lai
            self.pm_laimax = config_obj.pm_laimax
            self.pm_laimin = config_obj.pm_laimin
            self.pm_tas = config_obj.pm_tas
            self.pm_tmin = config_obj.pm_tmin
            self.pm_rhs = config_obj.pm_rhs
            self.pm_wind = config_obj.pm_wind
            self.pm_rsds = config_obj.pm_rsds
            self.pm_rlds = config_obj.pm_rlds
            self.pm_lct = config_obj.pm_lct
            self.pm_elev = config_obj.pm_elev
            nmonths = config_obj.nmonths

        super().__init__(nmonths=nmonths)

        # values from literature
        et_params = np.genfromtxt(self.pm_params, delimiter=',')

        self.cL = et_params[:, 0]
        self.beta = et_params[:, 1]
        self.rslimit = et_params[:, 2]

        # correlation coefficient for calculating emissivity (range from 0.34 to 0.44), not the ABCD parameter a
        self.ae = et_params[:, 3]

        # correlation coefficient for calculating emissivity (range from -0.14 to -0.25), not the ABCD parameter b
        self.be = et_params[:, 4]
        self.Tminopen = et_params[:, 5]
        self.Tminclose = et_params[:, 6]
        self.VPDclose = et_params[:, 7]
        self.VPDopen = et_params[:, 8]
        self.RBLmin = et_params[:, 9]
        self.RBLmax = et_params[:, 10]
        self.rc = et_params[:, 11]
        self.emiss = et_params[:, 12]

        # 2-d, rows:11 land cover types, cols:12 months
        self.alpha = np.genfromtxt(self.pm_alpha, delimiter=',')
        self.lai = np.genfromtxt(self.pm_lai, delimiter=',')
        self.laimax = np.genfromtxt(self.pm_laimax, delimiter=',')
        self.laimin = np.genfromtxt(self.pm_laimin, delimiter=',')

        # convert missing values to 0
        self.tair_load = self.load_to_array(self.pm_tas, 'pm_tas', nan_to_num=True)
        self.TMIN_load = self.load_to_array(self.pm_tmin, 'pm_tmin', nan_to_num=True)
        self.rhs_load = self.load_to_array(self.pm_rhs, 'pm_rhs', nan_to_num=True)
        self.wind_load = self.load_to_array(self.pm_wind, 'pm_wind', nan_to_num=True)
        self.rsds_load = self.load_to_array(self.pm_rsds, 'pm_rsds', nan_to_num=True)
        self.rlds_load = self.load_to_array(self.pm_rlds, 'pm_rlds', nan_to_num=True)

        # set previous air temp value leaving the first value at 0
        self.tairprev_load = np.zeros_like(self.tair_load)
        self.tairprev_load[1:, :] = self.tair_load[:-1, :]

        # use land cover for each target year (nan to 0);  1-d:67420 cells, 2-d: land cover type, 3-d: years
        self.lct_load = np.nan_to_num(np.load(self.pm_lct))

        # static data
        self.elev = np.nan_to_num(np.load(self.pm_elev))
