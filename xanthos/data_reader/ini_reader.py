"""
Read in settings from configuration file *.ini
Created on Oct 4, 2016

@author: lixi729
@email: xinya.li@pnl.gov
@Project: Xanthos V1.0


License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute

"""

import os
from configobj import ConfigObj


class ConfigReader:

    def __init__(self, ini):

        c = ConfigObj(ini)

        p = c['Project']

        # project dirs
        self.root = p['RootDir']
        self.ProjectName = p['ProjectName']
        self.InputFolder = os.path.join(self.root, p['InputFolder'])
        self.OutDir = self.create_dir(os.path.join(self.root, p['OutputFolder']))
        self.OutputFolder = self.create_dir(os.path.join(self.OutDir, self.ProjectName))

        # see if modules are in config file
        try:
            r = c['Reference']
            self.Reference = os.path.join(self.InputFolder, p['RefDir'])
        except KeyError:
            r = False

        try:
            pt = c['PET']
            self.PET = os.path.join(self.InputFolder, p['pet_dir'])
        except KeyError:
            pt = False

        try:
            ro = c['Runoff']
            self.RunoffDir = os.path.join(self.InputFolder, p['RunoffDir'])
        except KeyError:
            ro = False

        try:
            rt = c['Routing']
            self.RoutingDir = os.path.join(self.InputFolder, p['RoutingDir'])
        except KeyError:
            rt = False

        try:
            clm = c['Climate']
            self.ClimateFolder = os.path.join(self.InputFolder, p['ClimateDir'])
        except KeyError:
            clm = False

        try:
            d = c['Diagnostics']
            self.DiagDir = os.path.join(self.InputFolder, p['DiagDir'])
        except KeyError:
            d = False

        try:
            t = c['TimeSeriesPlot']
        except KeyError:
            t = False

        try:
            a = c['AccessibleWater']
            self.AccWatDir = os.path.join(self.InputFolder, p['AccWatDir'])
        except KeyError:
            a = False

        try:
            ha = c['HydropowerActual']
            self.HydActDir = os.path.join(self.InputFolder, p['HydActDir'])
        except KeyError:
            ha = False

        try:
            hp = c['HydropowerPotential']
        except KeyError:
            hp = False

        # project level settings
        self.ncell = 67420
        self.ngridrow = 360
        self.ngridcol = 720
        self.OutputUnitStr = None
        self.HistFlag = p['HistFlag']
        self.StartYear = int(p['StartYear'])
        self.EndYear = int(p['EndYear'])
        self.nmonths = (self.EndYear - self.StartYear + 1) * 12
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

        # pet
        if pt is not False:
            self.pet_module = pt['pet_module'].lower()

            if self.pet_module == 'hargreaves':
                pass

            elif self.pet_module == 'penman-monteith':
                pass

            elif self.pet_module == 'none':
                pass

            else:
                msg = "ERROR: PET module '{0}' not found. Please check spelling and try again.".format(self.pet_module)
                raise RuntimeError(msg)
        else:
            self.pet_module = 'none'

        print 'MODULE!!!!: ', self.pet_module

        # runoff
        if ro is not False:
            self.runoff_module = ro['runoff_module'].lower()

            if self.runoff_module == 'gwam':

                ro_mod = ro['GWAM']
                self.ro_model_dir = os.path.join(self.RunoffDir, ro_mod['model_dir'])
                self.SpinUp = int(ro_mod['SpinUp']) * 12
                self.MaxSoilMois = os.path.join(self.ro_model_dir, ro_mod['MaxSoilMois'])
                self.LakesMSM = os.path.join(self.ro_model_dir, ro_mod['LakesMSM'])
                self.AdditWaterMSM = os.path.join(self.ro_model_dir, ro_mod['AdditWaterMSM'])

                # channel storage file full path name with extension if running future
                self.ChStorageFile = None
                self.ChStorageVarName = None

                # soil moisture file full path name with extension if running future
                self.SavFile = None
                self.SavVarName = None

                if self.HistFlag == 'False':
                    try:
                        self.ChStorageFile = ro_mod['ChStorageFile']
                        self.ChStorageVarName = ro_mod['ChStorageVarName']
                        self.SavFile = ro_mod['SavFile']
                        self.SavVarName = ro_mod['SavVarName']

                    except KeyError:
                        raise RuntimeError("Error: ChStorageFile and ChStorageVarName are not defined for Future Mode.")

            elif self.runoff_module == 'abcd':

                ro_mod = ro['abcd']
                self.ro_model_dir = os.path.join(self.RunoffDir, ro_mod['model_dir'])
                self.ro_out_dir = self.create_dir(os.path.join(self.OutputFolder, ro_mod['output_dir']))
                self.calib_file = os.path.join(self.ro_model_dir, ro_mod['calib_file'])
                self.SpinUp = int(ro_mod['SpinUp'])
                self.ro_jobs = int(ro_mod['jobs'])

            elif self.runoff_module == 'none':
                pass

            else:
                raise RuntimeError("ERROR: Runoff module '{0}' not found. Please check spelling and try again.".format(self.runoff_module))
        else:
            self.runoff_module = 'none'

        # routing
        if rt is not False:
            self.routing_module = rt['routing_module'].lower()

            if self.routing_module == 'mrtm':
                rt_mod = rt['MRTM']
                self.rt_model_dir = os.path.join(self.RoutingDir, rt_mod['model_dir'])
                self.strm_veloc = os.path.join(self.rt_model_dir, rt_mod['ChVeloc'])
                self.FlowDis = os.path.join(self.rt_model_dir, rt_mod['FlowDis'])
                self.FlowDir = os.path.join(self.rt_model_dir, rt_mod['FlowDir'])

                try:
                    self.routing_spinup = int(rt_mod['routing_spinup'])
                except KeyError:
                    self.routing_spinup = self.nmonths

                try:
                    self.alt_runoff = self.custom_runoff(rt_mod['alt_runoff'])
                except KeyError:
                    self.alt_runoff = None

            elif self.routing_module == 'none':
                pass

            else:
                raise RuntimeError("ERROR: Routing module '{0}' not found. Please check spelling and try again.".format(self.routing_module))

        else:
            self.routing_module = 'none'

        # create model configuration string (pet_runoff_routing)
        self.mod_cfg = '{0}_{1}_{2}'.format(self.pet_module, self.runoff_module, self.routing_module)

        if self.mod_cfg == 'none_none_none':
            raise RuntimeError('No PFT, Runoff, or Routing model selected.')

        # climate
        if clm is not False:
            self.climate = True
            self.OutputNameStr = clm['climate_scenario']
            self.PrecipitationFile = os.path.join(self.ClimateFolder, clm['PrecipitationFile'])
            self.PrecipVarName = clm['PrecipVarName']
            self.TemperatureFile = os.path.join(self.ClimateFolder, clm['TemperatureFile'])
            self.TempVarName = clm['TempVarName']
            self.DailyTemperatureRangeFile = os.path.join(self.ClimateFolder, clm['DailyTemperatureRangeFile'])
            self.DTRVarName = clm['DTRVarName']
        else:
            self.climate = False
            self.OutputNameStr = self.ProjectName.replace(' ', '_').replace(',', '_').replace(';', '_')
            print('WARNING:  No climate data selected for use.')

        # reference
        if r is not False:
            self.Area = os.path.join(self.Reference, r['Area'])
            self.Coord = os.path.join(self.Reference, r['Coord'])
            self.BasinIDs = os.path.join(self.Reference, r['BasinIDs'])
            self.BasinNames = os.path.join(self.Reference, r['BasinNames'])
            self.GCAMRegionIDs = os.path.join(self.Reference, r['GCAMRegionIDs'])
            self.GCAMRegionNames = os.path.join(self.Reference, r['GCAMRegionNames'])
            self.CountryIDs = os.path.join(self.Reference, r['CountryIDs'])
            self.CountryNames = os.path.join(self.Reference, r['CountryNames'])
        else:
            print('WARNING:  No reference data selected for use.')

        self.CellArea = None
        self.ChSlope = None
        self.DrainArea = None
        self.RiversMSM = None

        # diagnostics
        if d is not False:
            if self.PerformDiagnostics:
                self.VICDataFile = os.path.join(self.DiagDir, d['VICDataFile'])
                self.UNHDataFile = os.path.join(self.DiagDir, d['UNHDataFile'])
                self.WBMDataFile = os.path.join(self.DiagDir, d['WBMDataFile'])
                self.WBMCDataFile = os.path.join(self.DiagDir, d['WBMCDataFile'])
                self.DiagnosticScale = int(d['Scale'])

        # plots
        if t is not False:
            if self.CreateTimeSeriesPlot:
                self.TimeSeriesScale = int(t['Scale'])
                self.TimeSeriesMapID = 999

                try:
                    l = int(t['MapID'])
                    self.TimeSeriesMapID = l

                except:
                    # as list
                    l = map(int, t['MapID'])
                    self.TimeSeriesMapID = l

        # accessible water
        if a is not False:
            if self.CalculateAccessibleWater:
                self.ResCapacityFile = os.path.join(self.AccWatDir, a['ResCapacityFile'])
                self.BfiFile = os.path.join(self.AccWatDir, a['BfiFile'])
                self.HistEndYear = int(a['HistEndYear'])
                self.GCAM_StartYear = self.ck_year(int(a['GCAM_StartYear']))
                self.GCAM_EndYear = int(a['GCAM_EndYear'])
                self.GCAM_YearStep = int(a['GCAM_YearStep'])
                self.MovingMeanWindow = int(a['MovingMeanWindow'])
                self.Env_FlowPercent = float(a['Env_FlowPercent'])

                if (self.StartYear > self.GCAM_StartYear) or (self.EndYear < self.GCAM_EndYear):
                    raise RuntimeError("Accessible water range of GCAM years are outside the range of years in climate data.")

        # hydropower potential
        if hp is not False:
            if self.CalculateHydropowerPotential:
                self.hpot_start_date = hp['hpot_start_date']
                self.q_ex = float(hp['q_ex'])
                self.ef = float(hp['ef'])
                self.GridData = os.path.join(self.HydActDir, ha['GridData'])

        # hydropower actual
        if ha is not False:
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
            raise RuntimeError("Accessible water year {0} is outside the range of years in the climate data.".format(yr))
        else:
            return yr

    def custom_runoff(self, f):
        """
        Check for custom runoff file name.  If 'none', return None; else
        return full path to file.

        :param f:
        :return:
        """
        if f == 'none':
            return None
        else:
            return os.path.join(self.rt_model_dir, f)

    def create_dir(self, pth):
        """
        Check to see if the target path is exists.
        """
        if os.path.isdir(pth) is False:
            os.mkdir(pth)
        return pth

    def log_info(self):

        print('ProjectName : {}'.format(self.ProjectName))
        print('InputFolder : {}'.format(self.InputFolder))
        print('OutputFolder: {}'.format(self.OutputFolder))
        print('StartYear - End Year        : {0}-{1}'.format(self.StartYear, self.EndYear))
        print('Number of Months            : {}'.format(self.nmonths))

        if self.climate:
            print('Precipitation File          : {}'.format(self.PrecipitationFile))
            print('Temperature File            : {}'.format(self.TemperatureFile))
            print('Daily Temperature Range File: {}'.format(self.DailyTemperatureRangeFile))

        if self.HistFlag.lower() in ['true', 't', 'yes', 'y', '1']:
            print('Running: Historic Mode')

        else:
            print('Running:  Future Mode')
            try:
                print('Historic Soil Moisture File: {}'.format(self.SavFile))
                print('Historic Channel Storage File: {}'.format(self.ChStorageFile))
            except AttributeError:
                pass

        try:
            print('Diagnostics will be performed using the data file: {}'.format(self.VICDataFile))
        except AttributeError:
            pass
