'''
Read in settings from configuration file *.ini
Created on Oct 4, 2016

@author: lixi729
@email: xinya.li@pnl.gov
@Project: Xanthos V1.0


License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute

'''
import sys
from configobj import ConfigObj  # install configobj package in Python
from ConfigSettings import ConfigSettings


def getSimulatorSettings(iniFile):
    config = ConfigObj(iniFile)
    settings = ConfigSettings()

    settings.ProjectName = config['Project']['ProjectName']
    settings.InputFolder = AddSlashToDir(config['Project']['InputFolder'])
    settings.OutputFolder = AddSlashToDir(config['Project']['OutputFolder'])
    settings.OutputFormat = int(config['Project']['OutputFormat'])
    settings.OutputUnit = int(config['Project']['OutputUnit'])
    settings.OutputInYear = int(config['Project']['OutputInYear'])
    settings.AggregateRunoffBasin = int(config['Project']['AggregateRunoffBasin'])
    settings.AggregateRunoffCountry = int(config['Project']['AggregateRunoffCountry'])
    settings.AggregateRunoffGCAMRegion = int(config['Project']['AggregateRunoffGCAMRegion'])
    settings.PerformDiagnostics = int(config['Project']['PerformDiagnostics'])
    settings.CreateTimeSeriesPlot = int(config['Project']['CreateTimeSeriesPlot'])
    settings.CalculateAccessibleWater = int(config['Project']['CalculateAccessibleWater'])

    settings.OutputNameStr, config['Climate']['HistFlag'], settings.StartYear, settings.EndYear = CheckClimateDataNames(
        config)
    settings.nmonths = int((settings.EndYear - settings.StartYear + 1) * 12)

    ClimateFolder = AddSlashToDir(config['Climate']['Folder'])
    settings.PrecipitationFile = ClimateFolder + config['Climate']['PrecipitationFile']
    settings.PrecipVarName = config['Climate']['PrecipVarName']
    settings.TemperatureFile = ClimateFolder + config['Climate']['TemperatureFile']
    settings.TempVarName = config['Climate']['TempVarName']
    settings.DailyTemperatureRangeFile = ClimateFolder + config['Climate']['DailyTemperatureRangeFile']
    settings.DTRVarName = config['Climate']['DTRVarName']
    settings.HistFlag = config['Climate']['HistFlag']
    settings.SpinUp = int(config['Climate']['SpinUp'])
    if settings.HistFlag == 'False':
        try:
            settings.ChStorageFile = config['Climate']['ChStorageFile']
            settings.ChStorageVarName = config['Climate']['ChStorageVarName']
            settings.SavFile = config['Climate']['SavFile']
            settings.SavVarName = config['Climate']['SavVarName']
        except:
            print("Error! ChStorageFile and ChStorageVarName are not defined for Future Mode.")
            sys.exit()

    settings.Area = settings.InputFolder + config['GriddedMap']['Area']
    settings.Coord = settings.InputFolder + config['GriddedMap']['Coord']
    settings.FlowDis = settings.InputFolder + config['GriddedMap']['FlowDis']
    settings.FlowDir = settings.InputFolder + config['GriddedMap']['FlowDir']
    settings.BasinIDs = settings.InputFolder + config['GriddedMap']['BasinIDs']
    settings.BasinNames = settings.InputFolder + config['GriddedMap']['BasinNames']
    settings.GCAMRegionIDs = settings.InputFolder + config['GriddedMap']['GCAMRegionIDs']
    settings.GCAMRegionNames = settings.InputFolder + config['GriddedMap']['GCAMRegionNames']
    settings.CountryIDs = settings.InputFolder + config['GriddedMap']['CountryIDs']
    settings.CountryNames = settings.InputFolder + config['GriddedMap']['CountryNames']
    settings.MaxSoilMois = settings.InputFolder + config['GriddedMap']['MaxSoilMois']
    settings.LakesMSM = settings.InputFolder + config['GriddedMap']['LakesMSM']
    settings.AdditWaterMSM = settings.InputFolder + config['GriddedMap']['AdditWaterMSM']

    if settings.PerformDiagnostics:
        settings.VICDataFile = config['Diagnostics']['VICDataFile']
        settings.UNHDataFile = config['Diagnostics']['UNHDataFile']
        settings.WBMDataFile = config['Diagnostics']['WBMDataFile']
        settings.WBMCDataFile = config['Diagnostics']['WBMCDataFile']
        settings.DiagnosticScale = int(config['Diagnostics']['Scale'])

    if settings.CreateTimeSeriesPlot:
        settings.TimeSeriesScale = int(config['TimeSeriesPlot']['Scale'])
        try:
            l = int(config['TimeSeriesPlot']['MapID'])
        except:
            l = map(int, config['TimeSeriesPlot']['MapID'])  # list
        settings.TimeSeriesMapID = l

    if settings.CalculateAccessibleWater:
        settings.ResCapacityFile = config['AccessibleWater']['ResCapacityFile']
        settings.BfiFile = config['AccessibleWater']['BfiFile']
        settings.HistEndYear = int(config['AccessibleWater']['HistEndYear'])
        settings.GCAM_StartYear = int(config['AccessibleWater']['GCAM_StartYear'])
        settings.GCAM_EndYear = int(config['AccessibleWater']['GCAM_EndYear'])
        settings.GCAM_YearStep = int(config['AccessibleWater']['GCAM_YearStep'])
        settings.MovingMeanWindow = int(config['AccessibleWater']['MovingMeanWindow'])
        settings.Env_FlowPercent = float(config['AccessibleWater']['Env_FlowPercent'])

        if settings.StartYear > settings.GCAM_StartYear or settings.EndYear < settings.GCAM_EndYear:
            print("Error! For accessible water, range of GCAM years are outside the range of years in climate data.")
            sys.exit()

    return settings


def PrintInfo(settings):
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


def CheckClimateDataNames(config):
    a = config['Climate']['PrecipitationFile'].split(".")[0]
    b = config['Climate']['TemperatureFile'].split(".")[0]
    c = config['Climate']['DailyTemperatureRangeFile'].split(".")[0]
    flag = config['Climate']['HistFlag']

    a = "_".join(a.split("_")[1:])
    b = "_".join(b.split("_")[1:])
    c = "_".join(c.split("_")[1:])

    if a == b == c:
        #             if (not "histor" in a) and (flag.lower() in ['true', 't', 'yes', 'y', '1']) :
        #                 print "Warning! Climate data should be in Future Mode. HistFlag = False in ini file."
        #                 flag = 'False'
        #             elif ("histor" in a) and (flag.lower() in ['false', 'f', 'no', 'n', '0']) :
        #                 print "Warning! Climate data should be in Historic Mode. HistFlag = True in ini file."
        #                 flag = 'True'
        if (flag.lower() in ['true', 't', 'yes', 'y', '1']):
            flag = 'True'
        elif (flag.lower() in ['false', 'f', 'no', 'n', '0']):
            flag = 'False'

        startyear = int(a.split("_")[-2][:4])
        endyear = int(a.split("_")[-1][:4])

        return a, flag, startyear, endyear

    else:
        print "Error! Precipitation, Temperature and Daily Temperature Range files are not from the same case ."
        sys.exit()


def AddSlashToDir(string):
    string = string.rstrip('/') + '/'

    return string