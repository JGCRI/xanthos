import numpy as np
from xanthos.data_reader.data_reference import DataReference


class DataCalibrationManaged(DataReference):
    """Load data for calibration that uses streamflow and accounts for water management."""

    def __init__(self,
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
                 end_year=None):

        if config_obj is None:

            self.start_year = start_year
            self.end_year = end_year
            self.nmonths = (self.end_year - self.start_year + 1) * 12

            super().__init__(nmonths=self.nmonths)

            # use basin-level flow as target for calibration; select only columns for basin number and runoff
            self.cal_obs = self.load_data(cal_observed, 0)[:, [0, 3]]

            # load dam and other input data
            self.purpose = np.load(purpose_file)
            self.capacity = np.load(capacity_file)
            self.hp_release = np.load(hp_release_file)
            self.water_consumption = np.load(water_consumption_file)
            self.instream_flow_natural = np.load(instream_flow_natural_file)
            self.ini_channel_storage = np.load(initial_channel_storage_natural_file)
            self.sm = np.load(sm_file)
            self.mtif_natural = np.load(mtif_natural_file)
            self.maxtif_natural = np.load(maxtif_natural_file)
            self.total_demand_cumecs = np.load(total_demand_cumecs_file)
            self.grdc_coord_index_file = np.load(grdc_coord_index_file)

        else:

            self.config_obj = config_obj
            self.start_year = config_obj.StartYear
            self.end_year = config_obj.EndYear
            self.nmonths = config_obj.nmonths

            super().__init__(nmonths=self.nmonths)

            # use basin-level flow as target for calibration; select only columns for basin number and runoff
            self.cal_obs = self.load_data(self.config_obj.cal_observed, 0)[:, [0, 3]]

            # load dam and other input data
            self.purpose = np.load(self.config_obj.purpose_file)
            self.capacity = np.load(self.config_obj.capacity_file)
            self.hp_release = np.load(self.config_obj.hp_release_file)
            self.water_consumption = np.load(self.config_obj.water_consumption_file)
            self.instream_flow_natural = np.load(self.config_obj.instream_flow_natural_file)
            self.ini_channel_storage = np.load(self.config_obj.initial_channel_storage_natural_file)
            self.sm = np.load(self.config_obj.sm_file)
            self.mtif_natural = np.load(self.config_obj.mtif_natural_file)
            self.maxtif_natural = np.load(self.config_obj.maxtif_natural_file)
            self.total_demand_cumecs = np.load(self.config_obj.total_demand_cumecs_file)
            self.grdc_coord_index_file = np.load(self.config_obj.grdc_coord_index_file)
