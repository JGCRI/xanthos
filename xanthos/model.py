"""
Model interface for Xanthos

@author   Chris R. Vernon
@email:   chris.vernon@pnnl.gov
@Project: Xanthos 2.0

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute
"""

import os
import sys
from xanthos.data_reader.ini_reader import ConfigReader
from xanthos.utils.logging import Logger
import xanthos.configurations as mods


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

    # run from terminal
    args = sys.argv[1:]

    if len(args) > 1:
        print('USAGE:  One argument should be passed. Full path file name with extension for config file.')
        print('Exiting...')
        raise IOError()

    ini = args[0]

    if os.path.isfile is False:
        print('ERROR:  Config file not found.')
        print('You entered:  {0}'.format(ini))
        print('Please enter a full path file name with extension to config file and retry.')
        raise IOError()

    xth = Xanthos(ini=ini)
    xth.execute()
    del xth
