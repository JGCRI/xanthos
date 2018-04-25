"""
Model interface for Xanthos

@author   Chris R. Vernon
@email:   chris.vernon@pnnl.gov
@Project: Xanthos 2.0

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute
"""

import argparse
import os
import sys
from data_reader.ini_reader import ConfigReader
from utils.logging import Logger
import configurations as mods


class Xanthos:

    def __init__(self, ini):
        self.ini = ini
        self.config = None
        self.log_file = None

    @staticmethod
    def make_dir(pth):
        """
        Create dir if not exists
        """
        if not os.path.exists(pth):
            os.makedirs(pth)

    def stage(self):
        """
        Set up run.
        """
        self.config = ConfigReader(self.ini)

        # create output directory
        self.make_dir(self.config.OutputFolder)

        # instantiate and write log file
        sys.stdout = Logger()
        self.log_file = os.path.join(self.config.OutputFolder, 'logfile.log')

    def execute(self):
        """
        Instantiate and write log file.
        """

        # stage data
        self.stage()

        with open(self.log_file, 'w') as sys.stdout.log:
            self.config.log_info()

            # run selected model configuration
            eval('mods.{0}(self.config)'.format(self.config.mod_cfg))

            print("End of {0}".format(self.config.ProjectName))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-config_file', type=str, help='Full path with file name to INI configuration file.')
    args = parser.parse_args()

    xth = Xanthos(args.config_file)
    xth.execute()
    del xth
