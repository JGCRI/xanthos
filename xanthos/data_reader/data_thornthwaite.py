from xanthos.data_reader.data_utils import DataUtils


class DataThornthwaite(DataUtils):

    def __init__(self, config):

        super().__init__(config_obj=config)

        self.tair = self.load_to_array(self.config.trn_tas, 'trn_tas', nan_to_num=True, warn_nan=True)
