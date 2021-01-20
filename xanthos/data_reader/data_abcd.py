import logging

from xanthos.data_reader.data_utils import DataUtils


class DataAbcd(DataUtils):

    def __init__(self, config=None, precip_file=None, precip_variable_name=None, tasmin_file=None,
                 tasmin_variable_name=None, start_yr=None, through_yr=None):

        # enable loading from arguments
        if config is None:
            self.precip_file = precip_file
            self.precip_variable_name = precip_variable_name
            self.tasmin_file = tasmin_file
            self.tasmin_variable_name = tasmin_variable_name
            nmonths = (through_yr - start_yr + 1) * 12

        else:
            self.precip_file = config.PrecipitationFile
            self.precip_variable_name = config.PrecipVarName
            self.tasmin_file = config.TempMinFile
            self.tasmin_variable_name = config.TempMinVarName
            nmonths = config.nmonths

        super().__init__(nmonths=nmonths)

        # monthly precipitation mm/mth
        self.precip = self.load_to_array(self.precip_file, var_name=self.precip_variable_name, warn_nan=True)

        # monthly average minimum daily temperature degree C (optional)
        if self.tasmin_file is None:
            logging.info('TempMinFile variable not found for the ABCD runoff module; Snowmelt will not be accounted for.')
            self.tmin = None

        else:
            self.tmin = self.load_to_array(self.tasmin_file, var_name=self.tasmin_variable_name, nan_to_num=True, warn_nan=True)
