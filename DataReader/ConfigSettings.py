'''
Created on Oct 4, 2016

@author: lixi729
@email: xinya.li@pnl.gov
@Project: Xanthos V1.0


License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute

'''


class ConfigSettings():
    def __init__(self):
        '''Project'''
        self.ProjectName = None
        self.InputFolder = None
        self.OutputFolder = None
        self.ncell = 67420
        self.ngridrow = 360
        self.ngridcol = 720
        self.OutputFormat = 0
        self.OutputUnit = 0
        self.OutputInYear = 0
        self.AggregateRunoffBasin = 0
        self.AggregateRunoffCountry = 0
        self.AggregateRunoffGCAMRegion = 0
        self.PerformDiagnostics = 0
        self.CreateTimeSeriesPlot = 0
        self.CalculateAccessibleWater = 0
        self.OutputUnitStr = None

        '''Climate'''
        self.PrecipitationFile = None
        self.PrecipVarName = None
        self.TemperatureFile = None
        self.TempVarName = None
        self.DailyTemperatureRangeFile = None
        self.DTRVarName = None
        self.HistFlag = None
        self.ChStorageFile = None
        self.ChStorageVarName = None
        self.SavFile = None
        self.SavVarName = None
        self.SpinUp = 0

        '''Map'''
        self.Area = None
        self.Coord = None
        self.CellArea = None
        self.FlowDis = None
        self.FlowDir = None
        self.ChSlope = None
        self.DrainArea = None
        self.BasinIDs = None
        self.BasinNames = None
        self.GCAMRegionIDs = None
        self.CountryIDs = None
        self.CountryNames = None
        self.MaxSoilMois = None
        self.LakesMSM = None
        self.RiversMSM = None
        self.AdditWaterMSM = None

        '''Diagnostics'''
        self.VICDataFile = None
        self.UNHDataFile = None
        self.WBMDataFile = None
        self.WBMCDataFile = None
        self.DiagnosticScale = 0

        '''TimeSeriesPlot'''
        self.TimeSeriesScale = 0
        self.TimeSeriesMapID = 999

        '''AccessibleWater'''
        self.ResCapacityFile = None
        self.BfiFile = None
        self.ResCapacityFile = None
        self.BfiFile = None
        self.HistEndYear = 2005
        self.GCAM_StartYear = 2005
        self.GCAM_EndYear = 2100
        self.GCAM_YearStep = 5
        self.MovingMeanWindow = 9
        self.Env_FlowPercent = 0.1

        '''Defined in code'''
        self.OutputNameStr = None
        self.StartYear = None
        self.EndYear = None
        self.nmonths = None
        self.mapindex = None