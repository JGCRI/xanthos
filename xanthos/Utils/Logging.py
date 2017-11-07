"""
Created on Jun 9, 2016
@author: lixi729

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute


This is the class to define logging process:

print command will write to screen and a logging file
"""

import sys


class Logger:

    def __init__(self):

        self.terminal = sys.stdout
        self.log = None

    def write(self, message):

        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        """
        Flush method is needed for python 3 compatibility.
        This handles the flush command by doing nothing
        """
        pass