from xanthos.data_reader.data_utils import DataUtils


class DataHargreaves(DataUtils):

    def __init__(self, config):

        super().__init__(config_obj=config)

        # average monthly temperature in degree C
        self.temp = self.load_to_array(self.config.TemperatureFile, var_name=self.config.TempVarName)

        # monthly average of daily temperature range in degree c
        self.dtr = self.load_to_array(self.config.DailyTemperatureRangeFile, var_name=self.config.DTRVarName, neg_to_zero=True)
