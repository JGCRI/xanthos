"""
Model interface for Xanthos.

@author   Chris R. Vernon
@email:   chris.vernon@pnnl.gov
@Project: Xanthos 2.0

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute
"""

import argparse
import os
import sys
import logging
from xanthos.data_reader.ini_reader import ConfigReader
import xanthos.configurations as mods


class Xanthos:
    """
    An extensible global hydrologic model.

    Xanthos is an open-source hydrologic model, designed to quantify
    and analyze global water availability. Xanthos simulates historical
    and future global water availability on a monthly time step at a
    spatial resolution of 0.5 geographic degrees.
    """

    def __init__(self, ini):
        """
        Initialize Xanthos.

        :param ini:     path to configuration file
        """
        self.ini = ini
        self.config = None

    @staticmethod
    def make_dir(pth):
        """Create dir if not exists."""
        if not os.path.exists(pth):
            os.makedirs(pth)

    def init_log(self):
        """
        Initialize project-wide logger.

        The logger outputs to both stdout and a file.
        """
        log_format = logging.Formatter('%(levelname)s: %(message)s')
        log_level = logging.DEBUG
        log_file = os.path.join(self.config.OutputFolder, 'logfile.log')

        logger = logging.getLogger()
        logger.setLevel(log_level)

        # logger console handler
        c_handler = logging.StreamHandler(sys.stdout)
        c_handler.setLevel(log_level)
        c_handler.setFormatter(log_format)
        logger.addHandler(c_handler)

        # logger file handler
        f_handler = logging.FileHandler(log_file)
        c_handler.setLevel(log_level)
        c_handler.setFormatter(log_format)
        logger.addHandler(f_handler)

    def stage(self, mem_args):
        """Set up run."""
        self.config = ConfigReader(self.ini)

        self.config.update(mem_args)

        # create output directory
        self.make_dir(self.config.OutputFolder)

        self.init_log()

    def execute(self, args={}):
        """
        Instantiate and write log file.

        @:param args:   Dictionary of config parameters
        """
        # stage data
        self.stage(args)

        self.config.log_info()

        # run selected model configuration
        if self.config.mod_cfg == 'hargreaves_gwam_mrtm':
            results = mods.hargreaves_gwam_mrtm(self.config)

        elif self.config.mod_cfg == 'hargreaves_abcd_mrtm':
            results = mods.hargreaves_abcd_mrtm(self.config)

        elif self.config.mod_cfg == 'pm_gwam_mrtm':
            results = mods.pm_gwam_mrtm(self.config)

        elif self.config.mod_cfg == 'pm_abcd_mrtm':
            results = mods.pm_abcd_mrtm(self.config)

        elif self.config.mod_cfg == 'hs_abcd_mrtm':
            results = mods.hs_abcd_mrtm(self.config)

        elif self.config.mod_cfg == 'thornthwaite_abcd_mrtm':
            results = mods.thornthwaite_abcd_mrtm(self.config)

        elif self.config.mod_cfg == 'none_gwam_mrtm':
            results = mods.none_gwam_mrtm(self.config)

        elif self.config.mod_cfg == 'none_abcd_mrtm':
            results = mods.none_abcd_mrtm(self.config)

        elif self.config.mod_cfg == 'none_none_mrtm':
            results = mods.none_none_mrtm(self.config)

        logging.info("End of {0}".format(self.config.ProjectName))

        return results


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-config_file', type=str, help='Full path with file name to INI configuration file.')
    args = parser.parse_args()

    xth = Xanthos(args.config_file)
    xth.execute()
    del xth
