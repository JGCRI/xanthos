import scipy.io as sio

from xanthos.data_reader.data_utils import DataUtils


class DataCalibration(DataUtils):
    """Load data for calibration that uses streamflow and accounts for water management."""

    def __init__(self, config_obj=None, cal_observed=None, purpose_file=None, capacity_file=None, hp_release_file=None,
                         water_consumption_file=None, instream_natural_flow_file=None,
                         initial_channel_storage_natural_file=None, sm_file=None, mtif_natural_file=None,
                         maxtif_natural_file=None, total_demand_cumecs_file=None,
                         grdc_xanthos_coord_index_file=None, start_year=None, end_year=None):

        if config_obj is None:

            self.start_year = start_year
            self.end_year = end_year
            self.nmonths = (self.end_year - self.start_year + 1) * 12

            super().__init__(nmonths=self.nmonths)

            # use basin-level flow as target for calibration; select only columns for basin number and runoff
            self.cal_obs = self.load_data(cal_observed, 0)[:, [0, 3]]

            # load dam and other input data
            self.purpose = sio.loadmat(purpose_file)['purpose3']
            self.capacity = sio.loadmat(capacity_file)['capacity']
            self.hp_release = sio.loadmat(hp_release_file)['HEP_Release']
            self.water_consumption = sio.loadmat(water_consumption_file)['WConsumption']
            self.instream_flow_natural = sio.loadmat(instream_natural_flow_file)['Initial_instream_flow_Natural']
            self.ini_channel_storage = sio.loadmat(initial_channel_storage_natural_file)['Initial_Channel_Storage_Natural']
            self.sm = sio.loadmat(sm_file)['SM']
            self.mtif_natural = sio.loadmat(mtif_natural_file)['mtif_natural']
            self.maxtif_natural = sio.loadmat(maxtif_natural_file)['maxtif_natural']
            self.total_demand_cumecs = sio.loadmat(total_demand_cumecs_file)['Total_Demand_cumecs']
            self.grdc_coord_index = sio.loadmat(grdc_xanthos_coord_index_file)['GRDC_xanthosCoordIndx']

        else:

            self.start_year = config_obj.StartYear
            self.end_year = config_obj.EndYear
            self.nmonths = config_obj.nmonths

            super().__init__(nmonths=self.nmonths)

            # use basin-level flow as target for calibration; select only columns for basin number and runoff
            self.cal_obs = self.load_data(self.config_obj.cal_observed, 0)[:, [0, 3]]

            # load dam and other input data
            self.purpose = sio.loadmat(config_obj.purpose_file)['purpose3']
            self.capacity = sio.loadmat(config_obj.capacity_file)['capacity']
            self.hp_release = sio.loadmat(config_obj.hp_release_file)['HEP_Release']
            self.water_consumption = sio.loadmat(config_obj.water_consumption_file)['WConsumption']
            self.instream_flow_natural = sio.loadmat(config_obj.instream_natural_flow_file)['Initial_instream_flow_Natural']
            self.ini_channel_storage = sio.loadmat(config_obj.initial_channel_storage_natural_file)['Initial_Channel_Storage_Natural']
            self.sm = sio.loadmat(config_obj.sm_file)['SM']
            self.mtif_natural = sio.loadmat(config_obj.mtif_natural_file)['mtif_natural']
            self.maxtif_natural = sio.loadmat(config_obj.maxtif_natural_file)['maxtif_natural']
            self.total_demand_cumecs = sio.loadmat(config_obj.total_demand_cumecs_file)['Total_Demand_cumecs']
            self.grdc_coord_index = sio.loadmat(config_obj.grdc_xanthos_coord_index_file)['GRDC_xanthosCoordIndx']
