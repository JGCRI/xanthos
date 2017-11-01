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


class ConfigReader:

    def __init__(self, ini):

        c = ConfigObj(ini)

        p = c['Project']
        m = c['Climate']
        r = c['Reference']
        g = c['Routing']
        d = c['Diagnostics']
        t = c['TimeSeriesPlot']
        a = c['AccessibleWater']
        ha = c['HydropowerActual']
        hp = c['HydropowerPotential']

        # project dirs
        self.root = p['RootDir']
        self.ProjectName = p['ProjectName']
        self.InputFolder = os.path.join(self.root, p['InputFolder'])
        self.OutputFolder = os.path.join(self.root, '{}/{}'.format(p['OutputFolder'], self.ProjectName))
        self.ClimateFolder = os.path.join(self.InputFolder, p['ClimateDir'])
        self.Reference = os.path.join(self.InputFolder, p['RefDir'])
        self.RoutingDir = os.path.join(self.InputFolder, p['RoutingDir'])
        self.DiagDir = os.path.join(self.InputFolder, p['DiagDir'])
        self.AccWatDir = os.path.join(self.InputFolder, p['AccWatDir'])
        self.HydActDir = os.path.join(self.InputFolder, p['HydActDir'])

        self.ncell = 67420
        self.ngridrow = 360
        self.ngridcol = 720
        self.OutputUnitStr = None

        # project settings
        self.HistFlag = p['HistFlag']
        self.StartYear = int(p['StartYear'])
        self.EndYear = int(p['EndYear'])
        self.nmonths = (self.EndYear - self.StartYear + 1) * 12
        self.SpinUp = int(p['SpinUp'])
        self.pet = p['pet']
        self.runoff = p['runoff']
        self.routing = p['routing']
        self.OutputFormat = int(p['OutputFormat'])
        self.OutputUnit = int(p['OutputUnit'])
        self.OutputInYear = int(p['OutputInYear'])
        self.AggregateRunoffBasin = int(p['AggregateRunoffBasin'])
        self.AggregateRunoffCountry = int(p['AggregateRunoffCountry'])
        self.AggregateRunoffGCAMRegion = int(p['AggregateRunoffGCAMRegion'])
        self.PerformDiagnostics = int(p['PerformDiagnostics'])
        self.CreateTimeSeriesPlot = int(p['CreateTimeSeriesPlot'])
        self.CalculateAccessibleWater = int(p['CalculateAccessibleWater'])
        self.CalculateHydropowerPotential = int(p['CalculateHydropowerPotential'])
        self.CalculateHydropowerActual = int(p['CalculateHydropowerActual'])

        # climate data
        self.OutputNameStr = m['ClimateScenario']
        self.PrecipitationFile = os.path.join(self.ClimateFolder, m['PrecipitationFile'])
        self.PrecipVarName = m['PrecipVarName']
        self.TemperatureFile = os.path.join(self.ClimateFolder, m['TemperatureFile'])
        self.TempVarName = m['TempVarName']
        self.DailyTemperatureRangeFile = os.path.join(self.ClimateFolder, m['DailyTemperatureRangeFile'])
        self.DTRVarName = m['DTRVarName']

        self.ChStorageFile = None
        self.ChStorageVarName = None
        self.SavFile = None
        self.SavVarName = None

        if self.HistFlag == 'False':
            try:
                self.ChStorageFile = m['ChStorageFile']
                self.ChStorageVarName = m['ChStorageVarName']
                self.SavFile = m['SavFile']
                self.SavVarName = m['SavVarName']

            except KeyError:
                print("Error: ChStorageFile and ChStorageVarName are not defined for Future Mode.")
                sys.exit()

        # reference
        self.Area = os.path.join(self.Reference, r['Area'])
        self.Coord = os.path.join(self.Reference, r['Coord'])
        self.FlowDis = os.path.join(self.Reference, r['FlowDis'])
        self.FlowDir = os.path.join(self.Reference, r['FlowDir'])
        self.BasinIDs = os.path.join(self.Reference, r['BasinIDs'])
        self.BasinNames = os.path.join(self.Reference, r['BasinNames'])
        self.GCAMRegionIDs = os.path.join(self.Reference, r['GCAMRegionIDs'])
        self.GCAMRegionNames = os.path.join(self.Reference, r['GCAMRegionNames'])
        self.CountryIDs = os.path.join(self.Reference, r['CountryIDs'])
        self.CountryNames = os.path.join(self.Reference, r['CountryNames'])
        self.MaxSoilMois = os.path.join(self.Reference, r['MaxSoilMois'])
        self.LakesMSM = os.path.join(self.Reference, r['LakesMSM'])
        self.AdditWaterMSM = os.path.join(self.Reference, r['AdditWaterMSM'])

        self.CellArea = None
        self.ChSlope = None
        self.DrainArea = None
        self.RiversMSM = None

        # routing
        self.ChVeloc = os.path.join(self.RoutingDir, g['ChVeloc'])

        # diagnostics
        if self.PerformDiagnostics:
            self.VICDataFile = os.path.join(self.DiagDir, d['VICDataFile'])
            self.UNHDataFile = os.path.join(self.DiagDir, d['UNHDataFile'])
            self.WBMDataFile = os.path.join(self.DiagDir, d['WBMDataFile'])
            self.WBMCDataFile = os.path.join(self.DiagDir, d['WBMCDataFile'])
            self.DiagnosticScale = int(d['Scale'])

        # plots
        if self.CreateTimeSeriesPlot:
            self.TimeSeriesScale = int(t['Scale'])

            try:
                l = int(t['MapID'])

            except:
                l = map(int, t['MapID'])  # list
                self.TimeSeriesMapID = l

        # accessible water
        if self.CalculateAccessibleWater:
            self.ResCapacityFile = os.path.join(self.AccWatDir, a['ResCapacityFile'])
            self.BfiFile = os.path.join(self.AccWatDir, a['BfiFile'])
            self.HistEndYear = int(a['HistEndYear'])
            self.GCAM_StartYear = self.ck_year(int(a['GCAM_StartYear']))
            self.GCAM_EndYear = int(a['GCAM_EndYear'])
            self.GCAM_YearStep = int(a['GCAM_YearStep'])
            self.MovingMeanWindow = int(a['MovingMeanWindow'])
            self.Env_FlowPercent = float(a['Env_FlowPercent'])

            if self.StartYear > self.GCAM_StartYear or self.EndYear < self.GCAM_EndYear:
                print("Error: Accessible water range of GCAM years are outside the range of years in climate data.")
                sys.exit()

        # hydropower potential
        if self.CalculateHydropowerPotential:
            self.hpot_start_date = hp['hpot_start_date']
            self.q_ex = float(hp['q_ex'])
            self.ef = float(hp['ef'])
            self.GridData = os.path.join(self.HydActDir, ha['GridData'])

        # hydropower actual
        if self.CalculateHydropowerActual:
            self.hact_start_date = ha['hact_start_date']
            self.HydroDamData = os.path.join(self.HydActDir, ha['HydroDamData'])
            self.MissingCap = os.path.join(self.HydActDir, ha['MissingCap'])
            self.rule_curves = os.path.join(self.HydActDir, ha['rule_curves'])
            self.GridData = os.path.join(self.HydActDir, ha['GridData'])
            self.DrainArea = os.path.join(self.HydActDir, ha['DrainArea'])

    def ck_year(self, yr):
        """
        Check to see if the target year is within the bounds of the data.
        """
        if (yr < self.StartYear) or (yr > self.EndYear):
            print("Error:  Accessible water year {0} is outside the range of years in the climate data.".format(yr))
            sys.exit()
        else:
            return yr

    def ck_path(self, pth):
        """
        Check to see if the target path is exists.
        """
        pass

    def log_info(self):

        print('ProjectName : {}'.format(self.ProjectName))
        print('InputFolder : {}'.format(self.InputFolder))
        print('OutputFolder: {}'.format(self.OutputFolder))
        print('Precipitation File          : {}'.format(self.PrecipitationFile))
        print('Temperature File            : {}'.format(self.TemperatureFile))
        print('Daily Temperature Range File: {}'.format(self.DailyTemperatureRangeFile))
        print('StartYear - End Year        : {0}-{1}'.format(self.StartYear, self.EndYear))
        print('Number of Months            : {}'.format(self.nmonths))

        if self.HistFlag.lower() in ['true', 't', 'yes', 'y', '1']:
            print('Running: Historic Mode')

        else:
            print('Running:  Future Mode')
            print('Historic Soil Moisture File: {}'.format(self.SavFile))
            print('Historic Channel Storage File: {}'.format(self.ChStorageFile))

        if self.SpinUp > 0:
            print('Spin-up     :Initialize the model using the first {} years'.format(self.SpinUp))

        if self.PerformDiagnostics:
            print('Diagnostics will be performed using the data file: {}'.format(self.VICDataFile))

