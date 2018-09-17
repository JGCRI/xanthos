"""
Read in settings from configuration file *.ini

@author: Xinya Li (xinya.li@pnl.gov) and Chris Vernon (chris.vernon@pnnl.gov)

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute
"""

import os
from configobj import ConfigObj


class ValidationException(Exception):
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)


class ConfigReader:

    def __init__(self, ini):

        c = ConfigObj(ini)

        p = c['Project']

        # project dirs
        self.root = p['RootDir']
        self.ProjectName = p['ProjectName']
        self.OutputNameStr = p['ProjectName']
        self.InputFolder = os.path.join(self.root, p['InputFolder'])
        self.OutDir = self.create_dir(os.path.join(self.root, p['OutputFolder']))
        self.OutputFolder = self.create_dir(os.path.join(self.OutDir, self.ProjectName))

        # -------------------------------------------------------------------
        # -------------------------------------------------------------------
        # COMPONENT MODULE CHECK
        # -------------------------------------------------------------------
        # -------------------------------------------------------------------

        # see if modules are in config file
        try:
            r = True
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

        try:
            cal = c['Calibrate']
        except KeyError:
            cal = False

        # -*****************************************************************-
        # ADD NEW COMPONENT MODULE CHECK HERE
        #
        # try:
        #     new_var = c['NewModule']
        # except KeyError:
        #     new_var = False
        #
        # -*****************************************************************-

        # -------------------------------------------------------------------
        # -------------------------------------------------------------------
        # PROJECT LEVEL SETTINGS
        # -------------------------------------------------------------------
        # -------------------------------------------------------------------

        # project level settings
        self.ncell = 67420
        self.ngridrow = 360
        self.ngridcol = 720
        self.n_basins = int(p['n_basins'])
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
        self.calibrate = int(p['Calibrate'])

        # -*****************************************************************-
        # GET MODULE RUN EXPECTATION FROM CONFIG
        #
        # self.RunNewModule = int(p['NewModule'])
        # -*****************************************************************-

        # -------------------------------------------------------------------
        # -------------------------------------------------------------------
        # CONFIGURE PET MODULE
        # -------------------------------------------------------------------
        # -------------------------------------------------------------------

        # PET config
        if pt:
            self.pet_module = pt['pet_module'].lower()

            if self.pet_module == 'hargreaves':
                pet_mod = pt['hargreaves']
                self.pet_dir = os.path.join(self.PET, pt['pet_dir'])

                # climate data
                try:
                    self.TemperatureFile = os.path.join(self.pet_dir, pet_mod['TemperatureFile'])
                except KeyError:
                    print('File path not provided for the TemperatureFile variable in the PET section of the config file.')
                    raise

                try:
                    self.TempVarName = pet_mod['TempVarName']
                except KeyError:
                    self.TempVarName = None

                try:
                    self.DailyTemperatureRangeFile = os.path.join(self.pet_dir, pet_mod['DailyTemperatureRangeFile'])
                except KeyError:
                    print('File path not provided for the DailyTemperatureRangeFile variable in the PET section of the config file.')
                    raise

                try:
                    self.DTRVarName = pet_mod['DTRVarName']
                except KeyError:
                    self.DTRVarName = None

            elif self.pet_module == 'hs':
                pet_mod = pt['hargreaves-samani']
                self.pet_dir = os.path.join(self.PET, pet_mod['pet_dir'])

                # climate data
                self.hs_tas = os.path.join(self.pet_dir, pet_mod['hs_tas'])
                self.hs_tmin = os.path.join(self.pet_dir, pet_mod['hs_tmin'])
                self.hs_tmax = os.path.join(self.pet_dir, pet_mod['hs_tmax'])

            elif self.pet_module == 'pm':
                pet_mod = pt['penman-monteith']
                self.pet_dir = os.path.join(self.PET, pet_mod['pet_dir'])

                # climate data
                self.pm_tas = os.path.join(self.pet_dir, pet_mod['pm_tas'])
                self.pm_tmin = os.path.join(self.pet_dir, pet_mod['pm_tmin'])
                self.pm_rhs = os.path.join(self.pet_dir, pet_mod['pm_rhs'])
                self.pm_rlds = os.path.join(self.pet_dir, pet_mod['pm_rlds'])
                self.pm_rsds = os.path.join(self.pet_dir, pet_mod['pm_rsds'])
                self.pm_wind = os.path.join(self.pet_dir, pet_mod['pm_wind'])

                # land cover data
                self.pm_lct = os.path.join(self.pet_dir, pet_mod['pm_lct'])
                self.pm_nlcs = int(pet_mod['pm_nlcs'])
                self.pm_water_idx = int(pet_mod['pm_water_idx'])
                self.pm_snow_idx = int(pet_mod['pm_snow_idx'])
                self.pm_lc_years = [int(i) for i in pet_mod['pm_lc_years']]

                # built-in data
                self.pm_params = os.path.join(self.pet_dir, 'gcam_ET_para.csv')
                self.pm_alpha = os.path.join(self.pet_dir, 'gcam_albedo.csv')
                self.pm_lai = os.path.join(self.pet_dir, 'gcam_lai.csv')
                self.pm_laimin = os.path.join(self.pet_dir, 'gcam_laimin.csv')
                self.pm_laimax = os.path.join(self.pet_dir, 'gcam_laimax.csv')
                self.pm_elev = os.path.join(self.pet_dir, 'elev.npy')

            elif self.pet_module == 'thornthwaite':
                pet_mod = pt['thornthwaite']
                self.pet_dir = os.path.join(self.PET, pet_mod['pet_dir'])

                # climate data
                self.trn_tas = os.path.join(self.pet_dir, pet_mod['trn_tas'])

            # -*****************************************************************-
            # CONDITIONAL FOR NEW PET MODULE
            #
            # elif self.pet_module == 'new_module':
            #     pet_mod = pt['hargreaves-samani']
            #     self.pet_dir = os.path.join(self.PET, pet_mod['pet_dir'])
            #
            #     # climate data
            #     self.hs_tas = os.path.join(self.pet_dir, pet_mod['hs_tas'])
            #     self.hs_tmin = os.path.join(self.pet_dir, pet_mod['hs_tmin'])
            #     self.hs_tmax = os.path.join(self.pet_dir, pet_mod['hs_tmax'])
            # -*****************************************************************-

            # use your own PET dataset
            elif self.pet_module == 'none':
                try:
                    self.pet_file = pt['pet_file']
                except KeyError:
                    raise "USAGE: Must provide a pet_file variable in the PET config section that contains the full path to an input PET file if not using an existing module."

            else:
                msg = "ERROR: PET module '{0}' not found. Please check spelling and try again.".format(self.pet_module)
                raise ValidationException(msg)
        else:
            self.pet_module = 'none'

            try:
                self.pet_file = pt['pet_file']
            except KeyError:
                raise "USAGE: Must provide a pet_file variable in the PET config section that contains the full path to an input PET file if not using an existing module."

        # -------------------------------------------------------------------
        # -------------------------------------------------------------------
        # CONFIGURE RUNOFF MODULE
        # -------------------------------------------------------------------
        # -------------------------------------------------------------------
        if ro:
            self.runoff_module = ro['runoff_module'].lower()

            if self.runoff_module == 'gwam':

                ro_mod = ro['GWAM']
                self.ro_model_dir = os.path.join(self.RunoffDir, ro_mod['runoff_dir'])
                self.runoff_spinup = int(ro_mod['runoff_spinup'])

                # built in files
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
                        raise ValidationException("Error: ChStorageFile and ChStorageVarName are not defined for Future Mode.")

                try:
                    self.PrecipitationFile = os.path.join(self.ro_model_dir, ro_mod['PrecipitationFile'])
                except KeyError:
                    print('File path not provided for the PrecipitationFile variable in the GCAM runoff section of the config file.')
                    raise

                try:
                    self.PrecipVarName = ro_mod['PrecipVarName']
                except KeyError:
                    self.PrecipVarName = None

                try:
                    self.TemperatureFile = os.path.join(self.ro_model_dir, ro_mod['TemperatureFile'])
                except KeyError:
                    print('File path not provided for the TemperatureFile variable in the GCAM runoff section of the config file.')
                    raise

                try:
                    self.TempVarName = ro_mod['TempVarName']
                except KeyError:
                    self.TempVarName = None

                try:
                    self.DailyTemperatureRangeFile = os.path.join(self.ro_model_dir, ro_mod['DailyTemperatureRangeFile'])
                except KeyError:
                    print('File path not provided for the DailyTemperatureRangeFile variable in the GCAM runoff section of the config file.')
                    raise

                try:
                    self.DTRVarName = ro_mod['DTRVarName']
                except KeyError:
                    self.DTRVarName = None

            elif self.runoff_module == 'abcd':

                ro_mod = ro['abcd']
                self.ro_model_dir = os.path.join(self.RunoffDir, ro_mod['runoff_dir'])
                self.calib_file = os.path.join(self.ro_model_dir, ro_mod['calib_file'])
                self.runoff_spinup = int(ro_mod['runoff_spinup'])
                self.ro_jobs = int(ro_mod['jobs'])

                try:
                    self.PrecipitationFile = ro_mod['PrecipitationFile']
                except KeyError:
                    print('File path not provided for the PrecipitationFile variable in the ABCD runoff section of the config file.')
                    raise


                try:
                    self.PrecipVarName = ro_mod['PrecipVarName']
                except KeyError:
                    self.PrecipVarName = None

                try:
                    self.TempMinFile = ro_mod['TempMinFile']
                except KeyError:
                    print('WARNING: TempMinFile variable not found in the ABCD runoff section of the config file; '
                          'results will not account for snow melt.')
                    self.TempMinFile = None

                try:
                    self.TempMinVarName = ro_mod['TempMinVarName']
                except KeyError:
                    self.TempMinVarName = None

            # -*****************************************************************-
            # CONDITIONAL FOR NEW RUNOFF MODULE
            #
            # elif self.runoff_module == 'new_module':
            #
            #     ro_mod = ro['abcd']
            #     self.ro_model_dir = os.path.join(self.RunoffDir, ro_mod['runoff_dir'])
            #     self.calib_file = os.path.join(self.ro_model_dir, ro_mod['calib_file'])
            #     self.runoff_spinup = int(ro_mod['runoff_spinup'])
            #     self.ro_jobs = int(ro_mod['jobs'])
            #
            #     try:
            #         self.PrecipitationFile = ro_mod['PrecipitationFile']
            #     except KeyError:
            #         print('File path not provided for the PrecipitationFile variable in the ABCD runoff section of the config file.')
            #         raise
            #
            #
            #     try:
            #         self.PrecipVarName = ro_mod['PrecipVarName']
            #     except KeyError:
            #         self.PrecipVarName = None
            #
            #     try:
            #         self.TempMinFile = ro_mod['TempMinFile']
            #     except KeyError:
            #         print('File path not provided for the TempMinFile variable in the ABCD runoff section of the config file.')
            #         raise
            #
            #     try:
            #         self.TempMinVarName = ro_mod['TempMinVarName']
            #     except KeyError:
            #         self.TempMinVarName = None
            # -*****************************************************************-

            elif self.runoff_module == 'none':
                pass

            else:
                raise ValidationException("ERROR: Runoff module '{0}' not found. Please check spelling and try again.".format(self.runoff_module))
        else:
            self.runoff_module = 'none'

        # -------------------------------------------------------------------
        # -------------------------------------------------------------------
        # CONFIGURE RUNOFF MODULE
        # -------------------------------------------------------------------
        # -------------------------------------------------------------------
        if rt:
            self.routing_module = rt['routing_module'].lower()

            if self.routing_module == 'mrtm':
                rt_mod = rt[self.routing_module]
                self.rt_model_dir = os.path.join(self.RoutingDir, rt_mod['routing_dir'])

                # load built-in files [ channel velocity, flow distance, flow direction ]
                self.strm_veloc = os.path.join(self.rt_model_dir, rt_mod['channel_velocity'])
                self.FlowDis = os.path.join(self.rt_model_dir, rt_mod['flow_distance'])
                self.FlowDir = os.path.join(self.rt_model_dir, rt_mod['flow_direction'])

                try:
                    self.routing_spinup = int(rt_mod['routing_spinup'])
                except KeyError:
                    self.routing_spinup = self.nmonths

                try:
                    self.alt_runoff = self.custom_runoff(rt_mod['alt_runoff'])
                except KeyError:
                    self.alt_runoff = None

            # -*****************************************************************-
            # CONDITIONAL FOR NEW ROUTING MODULE
            #
            # elif self.routing_module == 'new_module':
            #     rt_mod = rt[self.routing_module]
            #     self.rt_model_dir = os.path.join(self.RoutingDir, rt_mod['routing_dir'])
            #
            #     # load built-in files [ channel velocity, flow distance, flow direction ]
            #     self.strm_veloc = os.path.join(self.rt_model_dir, rt_mod['channel_velocity'])
            #     self.FlowDis = os.path.join(self.rt_model_dir, rt_mod['flow_distance'])
            #     self.FlowDir = os.path.join(self.rt_model_dir, rt_mod['flow_direction'])
            #
            #     try:
            #         self.routing_spinup = int(rt_mod['routing_spinup'])
            #     except KeyError:
            #         self.routing_spinup = self.nmonths
            #
            #     try:
            #         self.alt_runoff = self.custom_runoff(rt_mod['alt_runoff'])
            #     except KeyError:
            #         self.alt_runoff = None
            # -*****************************************************************-

            elif self.routing_module == 'none':
                pass

            else:
                raise ValidationException("ERROR: Routing module '{0}' not found. Please check spelling and try again.".format(self.routing_module))

        else:
            self.routing_module = 'none'

        # create model configuration string (pet_runoff_routing)
        self.mod_cfg = '{0}_{1}_{2}'.format(self.pet_module, self.runoff_module, self.routing_module)

        if self.mod_cfg == 'none_none_none':
            raise ValidationException('No PFT, Runoff, or Routing model selected.')

        # -------------------------------------------------------------------
        # -------------------------------------------------------------------
        # REFERENCE DATA
        # -------------------------------------------------------------------
        # -------------------------------------------------------------------
        if r:
            self.Area = os.path.join(self.Reference, 'Grid_Areas_ID.csv')
            self.Coord = os.path.join(self.Reference, 'coordinates.csv')
            self.BasinIDs = os.path.join(self.Reference, 'basin.csv')
            self.BasinNames = os.path.join(self.Reference, 'BasinNames235.txt')
            self.GCAMRegionIDs = os.path.join(self.Reference, 'region32_grids.csv')
            self.GCAMRegionNames = os.path.join(self.Reference, 'Rgn32Names.csv')
            self.CountryIDs = os.path.join(self.Reference, 'country.csv')
            self.CountryNames = os.path.join(self.Reference, 'country-names.csv')
        else:
            print('WARNING:  No reference data selected for use.')

        # -------------------------------------------------------------------
        # -------------------------------------------------------------------
        # OPTIONAL POST-PROCESSING MODULES
        # -------------------------------------------------------------------
        # -------------------------------------------------------------------

        # diagnostics
        if d:
            if self.PerformDiagnostics:
                self.VICDataFile = os.path.join(self.DiagDir, d['VICDataFile'])
                self.UNHDataFile = os.path.join(self.DiagDir, d['UNHDataFile'])
                self.WBMDataFile = os.path.join(self.DiagDir, d['WBMDataFile'])
                self.WBMCDataFile = os.path.join(self.DiagDir, d['WBMCDataFile'])
                self.DiagnosticScale = int(d['Scale'])

        # plots
        if t:
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
        if a:
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
                    raise ValidationException("Accessible water range of GCAM years are outside the range of years in climate data.")

        # hydropower potential
        if hp:
            if self.CalculateHydropowerPotential:
                self.hpot_start_date = hp['hpot_start_date']
                self.q_ex = float(hp['q_ex'])
                self.ef = float(hp['ef'])
                self.GridData = os.path.join(self.HydActDir, 'gridData.csv')

        # hydropower actual
        if ha:
            if self.CalculateHydropowerActual:
                self.hact_start_date = ha['hact_start_date']

                # built in data
                self.HydroDamData = os.path.join(self.HydActDir, 'resData_1593.csv')
                self.MissingCap = os.path.join(self.HydActDir, 'simulated_cap_by_country.csv')
                self.rule_curves = os.path.join(self.HydActDir, 'rule_curves_1593.npy')
                self.GridData = os.path.join(self.HydActDir, 'gridData.csv')
                self.DrainArea = os.path.join(self.HydActDir, 'DRT_half_SourceArea_globe_float.txt')

        # calibration mode
        if cal:
            if self.calibrate:
                self.set_calibrate = int(cal['set_calibrate'])
                self.cal_observed = cal['observed']
                self.obs_unit = self.ck_obs_unit(self.set_calibrate, cal['obs_unit'])
                self.calib_out_dir = self.create_dir(cal['calib_out_dir'])

        # -*****************************************************************-
        # CONDITIONAL FOR NEW RUNOFF MODULE
        #
        # if new_module:
        #     if self.calibrate:
        #         self.set_calibrate = int(cal['set_calibrate'])
        #         self.cal_observed = cal['observed']
        #         self.obs_unit = self.ck_obs_unit(self.set_calibrate, cal['obs_unit'])
        #         self.calib_out_dir = self.create_dir(cal['calib_out_dir'])
        # -*****************************************************************-

    @staticmethod
    def ck_obs_unit(set_calib, unit):
        """
        Checks the defined unit of the calibration data input.
        """
        valid_runoff = ('km3_per_mth', 'mm_per_mth')
        valid_streamflow = ('m3_per_sec')

        if set_calib == 0:

            if unit not in valid_runoff:
                raise ValidationException("Calibration data input units '{}' for runoff data not in required units '{}'".format(unit, valid_runoff))

            else:
                return unit

        elif set_calib == 1:

            if unit not in valid_streamflow:
                raise ValidationException("Calibration data input units '{}' for streamflow data not in required units '{}'".format(unit, valid_streamflow))

            else:
                return unit

    def ck_year(self, yr):
        """
        Check to see if the target year is within the bounds of the data.
        """
        if (yr < self.StartYear) or (yr > self.EndYear):
            raise ValidationException("Accessible water year {0} is outside the range of years in the climate data.".format(yr))
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

    @staticmethod
    def create_dir(pth):
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

    def update(self, args):
        """
        Overwrite configuration options

        :@param args:   Dictionary of parameters, where the key is the parameter name
        """
        for k, v in args.items():
            setattr(self, k, v)
