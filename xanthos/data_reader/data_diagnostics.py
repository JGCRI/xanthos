from xanthos.data_reader.data_utils import DataUtils


class DataDiagnostics(DataUtils):

    def __init__(self, config):

        super().__init__(config_obj=config)

        self.vic = self.load_data(self.config.VICDataFile, 0, "q")
        self.unh = self.load_data(self.config.UNHDataFile, 0, "q")
        self.wbmd = self.load_data(self.config.WBMDataFile, 0, 'q')
        self.wbmc = self.load_data(self.config.WBMCDataFile, 0, "q")
