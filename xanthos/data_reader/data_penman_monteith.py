import numpy as np
from xanthos.data_reader.data_utils import DataUtils


class DataPenmanMonteith(DataUtils):

    def __init__(self, config):

        # TODO: remove config imports from data classes
        self.config = config

        # TODO create parameters that can be used instead of a config file
        self._pm_params = None


        # values from literature
        et_params = np.genfromtxt(self.config.pm_params, delimiter=',')

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
        self.alpha = np.genfromtxt(self.config.pm_alpha, delimiter=',')
        self.lai = np.genfromtxt(self.config.pm_lai, delimiter=',')
        self.laimax = np.genfromtxt(self.config.pm_laimax, delimiter=',')
        self.laimin = np.genfromtxt(self.config.pm_laimin, delimiter=',')

        # convert missing values to 0
        self.tair_load = self.load_to_array(self.config.pm_tas, 'pm_tas', nan_to_num=True)
        self.TMIN_load = self.load_to_array(self.config.pm_tmin, 'pm_tmin', nan_to_num=True)
        self.rhs_load = self.load_to_array(self.config.pm_rhs, 'pm_rhs', nan_to_num=True)
        self.wind_load = self.load_to_array(self.config.pm_wind, 'pm_wind', nan_to_num=True)
        self.rsds_load = self.load_to_array(self.config.pm_rsds, 'pm_rsds', nan_to_num=True)
        self.rlds_load = self.load_to_array(self.config.pm_rlds, 'pm_rlds', nan_to_num=True)

        # set previous air temp value leaving the first value at 0
        self.tairprev_load = np.zeros_like(self.tair_load)
        self.tairprev_load[1:, :] = self.tair_load[:-1, :]

        # use land cover for each target year (nan to 0);  1-d:67420 cells, 2-d: land cover type, 3-d: years
        self.lct_load = np.nan_to_num(np.load(self.config.pm_lct))

        # static data
        self.elev = np.nan_to_num(np.load(self.config.pm_elev))
