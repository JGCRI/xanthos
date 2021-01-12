import scipy.io as sio


class DataCalibration:
    """Load data for calibration"""

    def __init__(self, config_obj=None, cal_observed=None, purpose_file=None, capacity_file=None, hp_release_file=None,
                         water_consumption_file=None, instream_natural_flow_file=None,
                         initial_channel_storage_natural_file=None, sm_file=None, mtif_natural_file=None,
                         maxtif_natural_file=None, total_demand_cumecs_file=None,
                         grdc_xanthos_coord_index_file=None):

        if config_obj is None:

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



