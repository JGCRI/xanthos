import logging

from xanthos.data_reader.data_utils import DataUtils


class DataAbcd(DataUtils):

    def __init__(self, config_obj=None, precip_file=None, precip_variable_name=None, tasmin_file=None,
                 tasmin_variable_name=None):

        # enable loading from arguments
        if config_obj is None:
            self.precip_file = precip_file
            self.precip_variable_name = precip_variable_name
            self.tasmin_file = tasmin_file
            self.tasmin_variable_name = tasmin_variable_name

        else:
            self.precip_file = config_obj.precip_file
            self.precip_variable_name = config_obj.precip_variable_name
            self.tasmin_file = config_obj.tasmin_file
            self.tasmin_variable_name = config_obj.tasmin_variable_name

        # monthly precipitation mm/mth
        self.precip = self.load_to_array(self.precip_file, var_name=self.precip_variable_name, warn_nan=True)

        # monthly average minimum daily temperature degree C (optional)
        if self.tasmin_file is None:
            logging.info('TempMinFile variable not found for the ABCD runoff module; Snowmelt will not be accounted for.')
            self.tmin = None

        else:
            self.tmin = self.load_to_array(self.tasmin_file, var_name=self.tasmin_variable_name, nan_to_num=True, warn_nan=True)
