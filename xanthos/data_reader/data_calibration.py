from xanthos.data_reader.data_utils import DataUtils


class DataCalibration(DataUtils):
    """Load data for calibration"""

    def __init__(self, config_obj=None, cal_observed=None, start_year=None, end_year=None):

        if config_obj is None:
            self.start_year = start_year
            self.end_year = end_year
            self.nmonths = (self.end_year - self.start_year + 1) * 12

            super().__init__(nmonths=self.nmonths)

            # use basin-level flow as target for calibration; select only columns for basin number and runoff
            self.cal_obs = self.load_data(cal_observed, 0)[:, [0, 3]]

        else:
            self.start_year = config_obj.StartYear
            self.end_year = config_obj.EndYear
            self.nmonths = config_obj.nmonths

            super().__init__(nmonths=self.nmonths)

            # use basin-level flow as target for calibration; select only columns for basin number and runoff
            self.cal_obs = self.load_data(config_obj.cal_observed, 0)[:, [0, 3]]
