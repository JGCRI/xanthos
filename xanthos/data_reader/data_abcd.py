import logging

from xanthos.data_reader.data_utils import DataUtils


class DataAbcd(DataUtils):

    def __init__(self, config):

        super().__init__(config_obj=config)

        # monthly precipitation mm/mth
        self.precip = self.load_to_array(self.config.PrecipitationFile, var_name=self.config.PrecipVarName, warn_nan=True)

        # monthly average minimum daily temperature degree C (optional)
        if self.config.TempMinFile is None:
            logging.info('TempMinFile variable not found for the ABCD runoff module; Snowmelt will not be accounted for.')
            self.tmin = None

        else:
            self.tmin = self.load_to_array(self.config.TempMinFile, var_name=self.config.TempMinVarName, nan_to_num=True, warn_nan=True)
