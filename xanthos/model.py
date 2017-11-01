"""
@Date: 10/01/2016
@author: Xinya Li (xinya.li@pnl.gov)
@Project: Xanthos V1.0

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute
"""

import os
import sys
from xanthos.data_reader.IniReader import ConfigReader
from xanthos.Utils.Logging import Logger
from gcam_hydro import Hydro


class Xanthos:

    def __init__(self, ini):

        self.ini = ini
        self.s = None
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
        self.s = ConfigReader(self.ini)

        # create output directory
        self.make_dir(self.s.OutputFolder)

        # instantiate and write log file
        sys.stdout = Logger()
        self.log_file = os.path.join(self.s.OutputFolder, 'logfile.log')

    def execute(self):
        """
        Instantiate and write log file.
        """

        # stage data
        self.stage()

        with open(self.log_file, 'w') as sys.stdout.log:

            # self.s.log_info()

            # instantiate Hydro class
            h = Hydro(self.s)

            # run model
            h.process()

            print("End of {0}".format(self.s.ProjectName))
