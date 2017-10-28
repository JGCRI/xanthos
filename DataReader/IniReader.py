'''
Read in settings from configuration file *.ini
Created on Oct 4, 2016

@author: lixi729
@email: xinya.li@pnl.gov
@Project: Xanthos V1.0


License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute

'''

import os
import sys
from configobj import ConfigObj
from ConfigSettings import ConfigSettings


#def getSimulatorSettings(iniFile):

class ConfigReader:

    def __init__(self, ini):

        c = ConfigObj(ini)

        p = c['Project']
        m = c['Climate']
        g = c['GriddedMap']
        d = c['Diagnostics']
        t = c['TimeSeriesPlot']
        a = c['AccessibleWater']

        # project dirs
        self.root = p['RootDir']
        self.ProjectName = p['ProjectName']
        self.InputFolder = os.path.join(self.root, p['InputFolder'])
        self.OutputFolder = os.path.join(self.root, '{}/{}'.format(p['OutputFolder'], self.ProjectName))
        self.ClimateFolder = os.path.join(self.InputFolder, p['ClimateDir'])

        # project settings
        self.HistFlag = p['HistFlag']
        self.StartYear = p['StartYear']
        self.EndYear = p['EndYear']
        self.nmonths = int((self.EndYear - self.StartYear + 1) * 12)
        self.SpinUp = int(p['SpinUp'])
        self.OutputFormat = int(p['OutputFormat'])
        self.OutputUnit = int(p['OutputUnit'])
        self.OutputInYear = int(p['OutputInYear'])
        self.AggregateRunoffBasin = int(p['AggregateRunoffBasin'])
        self.AggregateRunoffCountry = int(p['AggregateRunoffCountry'])
        self.AggregateRunoffGCAMRegion = int(p['AggregateRunoffGCAMRegion'])
        self.PerformDiagnostics = int(p['PerformDiagnostics'])
        self.CreateTimeSeriesPlot = int(p['CreateTimeSeriesPlot'])
        self.CalculateAccessibleWater = int(p['CalculateAccessibleWater'])

        # climate
        self.OutputNameStr = m['ClimateScenario']
        self.PrecipitationFile = os.path.join(self.ClimateFolder, m['PrecipitationFile'])
        self.PrecipVarName = m['PrecipVarName']
        self.TemperatureFile = os.path.join(self.ClimateFolder, m['TemperatureFile'])
        self.TempVarName = m['TempVarName']
        self.DailyTemperatureRangeFile = os.path.join(self.ClimateFolder, m['DailyTemperatureRangeFile'])
        self.DTRVarName = m['DTRVarName']

        if s.HistFlag == 'False':
            try:
                s.ChStorageFile = config['Climate']['ChStorageFile']
                s.ChStorageVarName = config['Climate']['ChStorageVarName']
                s.SavFile = config['Climate']['SavFile']
                s.SavVarName = config['Climate']['SavVarName']
            except:
                print("Error! ChStorageFile and ChStorageVarName are not defined for Future Mode.")
                sys.exit()

        s.Area = s.InputFolder + config['GriddedMap']['Area']
        s.Coord = s.InputFolder + config['GriddedMap']['Coord']
        s.FlowDis = s.InputFolder + config['GriddedMap']['FlowDis']
        s.FlowDir = s.InputFolder + config['GriddedMap']['FlowDir']
        s.BasinIDs = s.InputFolder + config['GriddedMap']['BasinIDs']
        s.BasinNames = s.InputFolder + config['GriddedMap']['BasinNames']
        s.GCAMRegionIDs = s.InputFolder + config['GriddedMap']['GCAMRegionIDs']
        s.GCAMRegionNames = s.InputFolder + config['GriddedMap']['GCAMRegionNames']
        s.CountryIDs = s.InputFolder + config['GriddedMap']['CountryIDs']
        s.CountryNames = s.InputFolder + config['GriddedMap']['CountryNames']
        s.MaxSoilMois = s.InputFolder + config['GriddedMap']['MaxSoilMois']
        s.LakesMSM = s.InputFolder + config['GriddedMap']['LakesMSM']
        s.AdditWaterMSM = s.InputFolder + config['GriddedMap']['AdditWaterMSM']

        if s.PerformDiagnostics:
            s.VICDataFile = config['Diagnostics']['VICDataFile']
            s.UNHDataFile = config['Diagnostics']['UNHDataFile']
            s.WBMDataFile = config['Diagnostics']['WBMDataFile']
            s.WBMCDataFile = config['Diagnostics']['WBMCDataFile']
            s.DiagnosticScale = int(config['Diagnostics']['Scale'])

        if s.CreateTimeSeriesPlot:
            s.TimeSeriesScale = int(config['TimeSeriesPlot']['Scale'])
            try:
                l = int(config['TimeSeriesPlot']['MapID'])
            except:
                l = map(int, config['TimeSeriesPlot']['MapID'])  # list
                s.TimeSeriesMapID = l

        if s.CalculateAccessibleWater:
            s.ResCapacityFile = config['AccessibleWater']['ResCapacityFile']
            s.BfiFile = config['AccessibleWater']['BfiFile']
            s.HistEndYear = int(config['AccessibleWater']['HistEndYear'])
            s.GCAM_StartYear = int(config['AccessibleWater']['GCAM_StartYear'])
            s.GCAM_EndYear = int(config['AccessibleWater']['GCAM_EndYear'])
            s.GCAM_YearStep = int(config['AccessibleWater']['GCAM_YearStep'])
            s.MovingMeanWindow = int(config['AccessibleWater']['MovingMeanWindow'])
            s.Env_FlowPercent = float(config['AccessibleWater']['Env_FlowPercent'])

            if s.StartYear > s.GCAM_StartYear or s.EndYear < s.GCAM_EndYear:
                print("Error! For accessible water, range of GCAM years are outside the range of years in climate data.")
                sys.exit()

        def PrintInfo(self, settings):

            print 'ProjectName :', settings.ProjectName
            print 'InputFolder :', settings.InputFolder
            print 'OutputFolder:', settings.OutputFolder
            print 'Precipitation File          :', settings.PrecipitationFile
            print 'Temperature File            :', settings.TemperatureFile
            print 'Daily Temperature Range File:', settings.DailyTemperatureRangeFile
            print 'StartYear - End Year        : ' + str(settings.StartYear) + " - " + str(settings.EndYear)
            print 'Number of Months            :', settings.nmonths

            if settings.HistFlag.lower() in ['true', 't', 'yes', 'y', '1']:
                print 'Historic Mode!'
            else:
                print 'Future Mode!'
                print 'Historic Soil Moisture File:', settings.SavFile
                print 'Historic Channel Storage File:', settings.ChStorageFile
            if settings.SpinUp > 0:
                print 'Spin-up     :Initialize the model using the first ' + str(settings.SpinUp) + ' years'
            if settings.PerformDiagnostics:
                print 'Diagnostics will be performed using the data file: ' + settings.VICDataFile

