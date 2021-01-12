from xanthos.data_reader.data_utils import DataUtils


class DataCalibration(DataUtils):
    """Load data for calibration"""

    def __init__(self, config):

        super().__init__(config_obj=config)

        # use basin-level flow as target for calibration; select only columns for basin number and runoff
        self.cal_obs = self.load_data(self.config.cal_observed, 0)[:, [0, 3]]
