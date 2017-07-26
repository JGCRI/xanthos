"""
@Date: 10/01/2016
@author: Xinya Li (xinya.li@pnl.gov)
@Project: Xanthos V1.0

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute
"""

import os, sys
# include the folder of source code
# sys.path.insert(0, './Source Code')
import DataReader.IniReader as IniReader
from Utils.Logging import Logger
from gcam_hydro import Hydro as GCAM_Hydro

# Read simulator settings from configuration file
settingFile = 'config_hist.ini' # Use config_hist.ini first for historical case; Then config_futu.ini for future case
settings = IniReader.getSimulatorSettings(settingFile)
# Check if outputFolder exist or not
if not os.path.exists(settings.OutputFolder):
    os.makedirs(settings.OutputFolder)

# Setup the log file, save it into the output folder
sys.stdout = Logger()
sys.stdout.log = open(settings.OutputFolder + "logfile.log", "w")

IniReader.PrintInfo(settings)

# Call the main function
GCAM_Hydro(settings)

print "End of", settings.ProjectName

sys.stdout.log.close()