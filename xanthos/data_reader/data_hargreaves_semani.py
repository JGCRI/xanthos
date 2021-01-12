from xanthos.data_reader.data_utils import DataUtils


class DataHargreavesSemani(DataUtils):

    def __init__(self, config):

        super().__init__(config_obj=config)

        self.hs_tas = self.load_to_array(self.config.hs_tas)
        self.hs_tmin = self.load_to_array(self.config.hs_tmin)
        self.hs_tmax = self.load_to_array(self.config.hs_tmax)
