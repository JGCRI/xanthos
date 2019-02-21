"""
Read in settings from configuration file *.ini.

@author: Xinya Li (xinya.li@pnl.gov), Chris R. Vernon (chris.vernon@pnnl.gov)
@Project: Xanthos V2.0

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute
"""

import os
import logging
from configobj import ConfigObj


class ValidationException(Exception):
    """Custom exception for invalid Xanthos inputs."""

    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)


class ConfigReader:
    """Read the Xanthos configuration .ini file."""

    def __init__(self, ini):
        """
        Load values from configuration file.

        :param ini:     path to the config file
        """
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
        try:
            ref = True
            self.Reference = os.path.join(self.InputFolder, p['RefDir'])
        except KeyError:
            ref = False

        try:
            pet_config = c['PET']
            self.PET = os.path.join(self.InputFolder, p['pet_dir'])
        except KeyError:
            pet_config = False

        try:
            runoff_config = c['Runoff']
            self.RunoffDir = os.path.join(self.InputFolder, p['RunoffDir'])
        except KeyError:
            runoff_config = False

        try:
            routing_config = c['Routing']
            self.RoutingDir = os.path.join(self.InputFolder, p['RoutingDir'])
        except KeyError:
            routing_config = False

        try:
            diagnostics_config = c['Diagnostics']
            self.DiagDir = os.path.join(self.InputFolder, p['DiagDir'])
        except KeyError:
            diagnostics_config = False

        try:
            timeseries_config = c['TimeSeriesPlot']
        except KeyError:
            timeseries_config = False

        try:
            drought_config = c['Drought']
        except KeyError:
            drought_config = False

        try:
            acc_water_config = c['AccessibleWater']
            self.AccWatDir = os.path.join(self.InputFolder, p['AccWatDir'])
        except KeyError:
            acc_water_config = False

        try:
            hydro_actual_config = c['HydropowerActual']
            self.HydActDir = os.path.join(self.InputFolder, p['HydActDir'])
        except KeyError:
            hydro_actual_config = False

        try:
            hydro_potential_config = c['HydropowerPotential']
        except KeyError:
            hydro_potential_config = False

        try:
            calibration_config = c['Calibrate']
        except KeyError:
            calibration_config = False

        # -------------------------------------------------------------------
        # -------------------------------------------------------------------
        # PROJECT LEVEL SETTINGS
        # -------------------------------------------------------------------
        # -------------------------------------------------------------------
        self.ncell = 67420
        self.ngridrow = 360
        self.ngridcol = 720
        self.n_basins = int(p['n_basins'])
        self.HistFlag = p['HistFlag']
        self.StartYear = int(p['StartYear'])
        self.EndYear = int(p['EndYear'])
        self.output_vars = p['output_vars']
        self.OutputFormat = int(p['OutputFormat'])
        self.OutputUnit = int(p['OutputUnit'])
        self.OutputInYear = int(p['OutputInYear'])
        self.AggregateRunoffBasin = int(p['AggregateRunoffBasin'])
        self.AggregateRunoffCountry = int(p['AggregateRunoffCountry'])
        self.AggregateRunoffGCAMRegion = int(p['AggregateRunoffGCAMRegion'])
        self.PerformDiagnostics = int(p['PerformDiagnostics'])
        self.CreateTimeSeriesPlot = int(p['CreateTimeSeriesPlot'])
        self.CalculateDroughtStats = int(p['CalculateDroughtStats'])
        self.CalculateAccessibleWater = int(p['CalculateAccessibleWater'])
        self.CalculateHydropowerPotential = int(p['CalculateHydropowerPotential'])
        self.CalculateHydropowerActual = int(p['CalculateHydropowerActual'])
        self.calibrate = int(p['Calibrate'])

        self.nmonths = (self.EndYear - self.StartYear + 1) * 12
        self.OutputUnitStr = '{}per{}'.format(('mm', 'km3')[self.OutputUnit], ('month', 'year')[self.OutputInYear])
        # Better to use this instead with force_list():
        #   http://www.voidspace.org.uk/python/articles/configobj.shtml#validation
        self.output_vars = [self.output_vars] if not isinstance(self.output_vars, list) else self.output_vars

        # -------------------------------------------------------------------
        # -------------------------------------------------------------------
        # REQUIRED MODULES
        # -------------------------------------------------------------------
        # -------------------------------------------------------------------

        # Load pet, runoff, and routing modules
        self.configure_pet(pet_config)
        self.configure_runoff(runoff_config)
        self.configure_routing(routing_config)

        # create model configuration string (pet_runoff_routing)
        self.mod_cfg = '{0}_{1}_{2}'.format(self.pet_module, self.runoff_module, self.routing_module)

        if self.mod_cfg == 'none_none_none':
            raise ValidationException('No PET, Runoff, or Routing model selected.')

        self.configure_reference_data(ref)

        # -------------------------------------------------------------------
        # -------------------------------------------------------------------
        # OPTIONAL POST-PROCESSING MODULES
        # -------------------------------------------------------------------
        # -------------------------------------------------------------------

        # diagnostics
        if diagnostics_config and self.PerformDiagnostics:
            self.configure_diagnostics(diagnostics_config)

        # plots
        if timeseries_config and self.CreateTimeSeriesPlot:
            self.configure_timeseries_plot(timeseries_config)

        # drought statistics
        if drought_config and self.CalculateDroughtStats:
            self.configure_drought_stats(drought_config)

        # accessible water
        if acc_water_config and self.CalculateAccessibleWater:
            self.configure_acc_water(acc_water_config)

        # hydropower potential
        if hydro_potential_config and self.CalculateHydropowerPotential:
            self.configure_hydropower_potential(hydro_potential_config)

        # hydropower actual
        if hydro_actual_config and self.CalculateHydropowerActual:
            self.configure_hydropower_actual(hydro_actual_config)

        # calibration mode
        if calibration_config and self.calibrate:
            self.configure_calibration(calibration_config)

    def configure_pet(self, pet_config):
        """Configure the PET module."""
        if not pet_config:
            self.pet_module = 'none'

            try:
                self.pet_file = pet_config['pet_file']
            except KeyError:
                raise ValidationException(
                    "USAGE: Must provide a pet_file variable in the PET config section that "
                    "contains the full path to an input PET file if not using an existing module."
                )

            return

        # Set up the PET module based on which module was selected
        self.pet_module = pet_config['pet_module'].lower()

        if self.pet_module == 'hargreaves':
            pet_mod = pet_config['hargreaves']
            self.pet_dir = os.path.join(self.PET, pet_mod['pet_dir'])

            # climate data
            try:
                self.TemperatureFile = os.path.join(self.pet_dir, pet_mod['TemperatureFile'])
            except KeyError:
                logging.exception("File path not provided for the TemperatureFile "
                                  "variable in the PET section of the config file.")
                raise

            try:
                self.TempVarName = pet_mod['TempVarName']
            except KeyError:
                self.TempVarName = None

            try:
                self.DailyTemperatureRangeFile = os.path.join(self.pet_dir, pet_mod['DailyTemperatureRangeFile'])
            except KeyError:
                logging.exception("File path not provided for the DailyTemperatureRangeFile "
                                  "variable in the PET section of the config file.")
                raise

            try:
                self.DTRVarName = pet_mod['DTRVarName']
            except KeyError:
                self.DTRVarName = None

        elif self.pet_module == 'hs':
            pet_mod = pet_config['hargreaves-samani']
            self.pet_dir = os.path.join(self.PET, pet_mod['pet_dir'])

            # climate data
            self.hs_tas = os.path.join(self.pet_dir, pet_mod['hs_tas'])
            self.hs_tmin = os.path.join(self.pet_dir, pet_mod['hs_tmin'])
            self.hs_tmax = os.path.join(self.pet_dir, pet_mod['hs_tmax'])

        elif self.pet_module == 'pm':
            pet_mod = pet_config['penman-monteith']
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
            pet_mod = pet_config['thornthwaite']
            self.pet_dir = os.path.join(self.PET, pet_mod['pet_dir'])

            # climate data
            self.trn_tas = os.path.join(self.pet_dir, pet_mod['trn_tas'])

        # use your own PET dataset
        elif self.pet_module == 'none':
            try:
                self.pet_file = pet_config['pet_file']
            except KeyError:
                raise ValidationException(
                    "USAGE: Must provide a pet_file variable in the PET config section that "
                    "contains the full path to an input PET file if not using an existing module."
                )

        else:
            raise ValidationException("ERROR: PET module '{0}' not found. Please check "
                                      "spelling and try again.".format(self.pet_module))

    def configure_runoff(self, runoff_config):
        """Configure runoff module."""
        # Return right away if runoff module not specified
        if not runoff_config:
            self.runoff_module = 'none'
            return

        self.runoff_module = runoff_config['runoff_module'].lower()

        if self.runoff_module == 'gwam':

            ro_mod = runoff_config['gwam']
            self.ro_model_dir = os.path.join(self.RunoffDir, ro_mod['runoff_dir'])
            self.runoff_spinup = int(ro_mod['runoff_spinup'])

            # built in files
            self.max_soil_moisture = os.path.join(self.ro_model_dir, ro_mod['max_soil_moisture'])
            self.lakes_msm = os.path.join(self.ro_model_dir, ro_mod['lakes_msm'])
            self.addit_water_msm = os.path.join(self.ro_model_dir, ro_mod['addit_water_msm'])

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
                    raise ValidationException("Error: ChStorageFile and ChStorageVarName "
                                              "are not defined for Future Mode.")

            try:
                self.PrecipitationFile = os.path.join(self.ro_model_dir, ro_mod['PrecipitationFile'])
            except KeyError:
                logging.exception("File path not provided for the PrecipitationFile variable "
                                  "in the GWAM runoff section of the config file.")
                raise

            try:
                self.PrecipVarName = ro_mod['PrecipVarName']
            except KeyError:
                self.PrecipVarName = None

        elif self.runoff_module == 'abcd':

            ro_mod = runoff_config['abcd']
            self.ro_model_dir = os.path.join(self.RunoffDir, ro_mod['runoff_dir'])
            self.calib_file = os.path.join(self.ro_model_dir, ro_mod['calib_file'])
            self.runoff_spinup = int(ro_mod['runoff_spinup'])
            self.ro_jobs = int(ro_mod['jobs'])

            try:
                self.PrecipitationFile = ro_mod['PrecipitationFile']
            except KeyError:
                logging.exception("File path not provided for the PrecipitationFile "
                                  "variable in the ABCD runoff section of the config file.")
                raise

            try:
                self.PrecipVarName = ro_mod['PrecipVarName']
            except KeyError:
                self.PrecipVarName = None

            try:
                self.TempMinFile = ro_mod['TempMinFile']
            except KeyError:
                self.TempMinFile = None

            try:
                self.TempMinVarName = ro_mod['TempMinVarName']
            except KeyError:
                self.TempMinVarName = None

        elif self.runoff_module == 'none':
            pass

        else:
            raise ValidationException("ERROR: Runoff module '{0}' not found. Please check "
                                      "spelling and try again.".format(self.runoff_module))

    def configure_routing(self, routing_config):
        """Configure routing module."""
        # Return right away if routing module not specified
        if not routing_config:
            self.routing_module = 'none'
            return

        self.routing_module = routing_config['routing_module'].lower()

        if self.routing_module == 'mrtm':
            rt_mod = routing_config[self.routing_module]
            self.rt_model_dir = os.path.join(self.RoutingDir, rt_mod['routing_dir'])

            # load built-in files [ channel velocity, flow distance, flow direction ]
            self.strm_veloc = os.path.join(self.rt_model_dir, rt_mod['channel_velocity'])
            self.flow_distance = os.path.join(self.rt_model_dir, rt_mod['flow_distance'])
            self.flow_direction = os.path.join(self.rt_model_dir, rt_mod['flow_direction'])

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
            raise ValidationException("ERROR: Routing module '{0}' not found. Please check "
                                      "spelling and try again.".format(self.routing_module))

    def configure_reference_data(self, ref):
        """Build reference data file paths."""
        if ref:
            self.Area = os.path.join(self.Reference, 'Grid_Areas_ID.csv')
            self.Coord = os.path.join(self.Reference, 'coordinates.csv')
            self.BasinIDs = os.path.join(self.Reference, 'basin.csv')
            self.BasinNames = os.path.join(self.Reference, 'BasinNames235.txt')
            self.GCAMRegionIDs = os.path.join(self.Reference, 'region32_grids.csv')
            self.GCAMRegionNames = os.path.join(self.Reference, 'Rgn32Names.csv')
            self.CountryIDs = os.path.join(self.Reference, 'country.csv')
            self.CountryNames = os.path.join(self.Reference, 'country-names.csv')
        else:
            logging.warning('No reference data selected for use.')

    def configure_diagnostics(self, diagnostics_config):
        """Configure diagnostics post-processing module."""
        self.VICDataFile = os.path.join(self.DiagDir, diagnostics_config['VICDataFile'])
        self.UNHDataFile = os.path.join(self.DiagDir, diagnostics_config['UNHDataFile'])
        self.WBMDataFile = os.path.join(self.DiagDir, diagnostics_config['WBMDataFile'])
        self.WBMCDataFile = os.path.join(self.DiagDir, diagnostics_config['WBMCDataFile'])
        self.DiagnosticScale = int(diagnostics_config['Scale'])

    def configure_timeseries_plot(self, timeseries_config):
        """Configure time series plots post-processing module."""
        self.TimeSeriesScale = int(timeseries_config['Scale'])
        self.TimeSeriesMapID = 999

        try:
            map_id = int(timeseries_config['MapID'])
            self.TimeSeriesMapID = map_id
        except TypeError:
            # as list
            map_id = list(map(int, timeseries_config['MapID']))
            self.TimeSeriesMapID = map_id

    def configure_drought_stats(self, drought_config):
        """Configure accessible water post-processing module."""
        self.drought_var = drought_config['drought_var']
        self.drought_thresholds = drought_config.get('drought_thresholds')  # optional

        # If no drought file is given, calculate threshold values
        if self.drought_thresholds is None:
            self.threshold_nper = int(drought_config['threshold_nper'])
            self.threshold_start_year = int(drought_config['threshold_start_year'])
            self.threshold_end_year = int(drought_config['threshold_end_year'])
            if (self.StartYear > self.threshold_start_year) or (self.EndYear < self.threshold_end_year):
                raise ValidationException("Drought threshold year range is outside the output year range.")

    def configure_acc_water(self, acc_water_config):
        """Configure accessible water post-processing module."""
        self.ResCapacityFile = os.path.join(self.AccWatDir, acc_water_config['ResCapacityFile'])
        self.BfiFile = os.path.join(self.AccWatDir, acc_water_config['BfiFile'])
        self.HistEndYear = int(acc_water_config['HistEndYear'])
        self.GCAM_StartYear = self.ck_year(int(acc_water_config['GCAM_StartYear']))
        self.GCAM_EndYear = int(acc_water_config['GCAM_EndYear'])
        self.GCAM_YearStep = int(acc_water_config['GCAM_YearStep'])
        self.MovingMeanWindow = int(acc_water_config['MovingMeanWindow'])
        self.Env_FlowPercent = float(acc_water_config['Env_FlowPercent'])

        if (self.StartYear > self.GCAM_StartYear) or (self.EndYear < self.GCAM_EndYear):
            raise ValidationException("Accessible water range of GCAM years are outside "
                                      "the range of years in climate data.")

    def configure_hydropower_potential(self, hydro_potential_config):
        """Configure hydropower potential post-processing module."""
        self.hpot_start_date = hydro_potential_config['hpot_start_date']
        self.q_ex = float(hydro_potential_config['q_ex'])
        self.ef = float(hydro_potential_config['ef'])
        self.GridData = os.path.join(self.HydActDir, 'gridData.csv')

    def configure_hydropower_actual(self, hydro_actual_config):
        """Configure hydropower actual post-processing module."""
        self.hact_start_date = hydro_actual_config['hact_start_date']

        # built in data
        self.HydroDamData = os.path.join(self.HydActDir, 'resData_1593.csv')
        self.MissingCap = os.path.join(self.HydActDir, 'simulated_cap_by_country.csv')
        self.rule_curves = os.path.join(self.HydActDir, 'rule_curves_1593.npy')
        self.GridData = os.path.join(self.HydActDir, 'gridData.csv')
        self.DrainArea = os.path.join(self.HydActDir, 'DRT_half_SourceArea_globe_float.txt')

    def configure_calibration(self, calibration_config):
        """Configure calibration settings."""
        self.set_calibrate = int(calibration_config['set_calibrate'])
        self.cal_observed = calibration_config['observed']
        self.obs_unit = self.ck_obs_unit(self.set_calibrate, calibration_config['obs_unit'])
        self.calib_out_dir = self.create_dir(calibration_config['calib_out_dir'])
        try:
            self.cal_basins = calibration_config['calibration_basins']
            # calibration basins can be specified as either a single
            # number (ex. 3), range (ex. 3-5), or a list (ex. 3,4,5)
            if type(self.cal_basins) is not list:
                self.cal_basins = [self.cal_basins]
        except KeyError:
            self.cal_basins = ['1-{}'.format(self.n_basins)]

    @staticmethod
    def ck_obs_unit(set_calib, unit):
        """Check the defined unit of the calibration data input."""
        valid_runoff = ('km3_per_mth', 'mm_per_mth')
        valid_streamflow = ('m3_per_sec')

        if set_calib == 0:

            if unit not in valid_runoff:
                raise ValidationException("Calibration data input units '{}' for runoff data "
                                          "not in required units '{}'".format(unit, valid_runoff))

            else:
                return unit

        elif set_calib == 1:

            if unit not in valid_streamflow:
                raise ValidationException("Calibration data input units '{}' for streamflow data "
                                          "not in required units '{}'".format(unit, valid_streamflow))

            else:
                return unit

    def ck_year(self, yr):
        """Check to see if the target year is within the bounds of the data."""
        if (yr < self.StartYear) or (yr > self.EndYear):
            raise ValidationException("Accessible water year {0} is outside the range "
                                      "of years in the climate data.".format(yr))
        else:
            return yr

    def custom_runoff(self, f):
        """
        Check for custom runoff file name.

        If 'none', return None; else return full path to file.

        :param f:
        :return:
        """
        if f == 'none':
            return None
        else:
            return os.path.join(self.rt_model_dir, f)

    @staticmethod
    def create_dir(pth):
        """Check to see if the target path is exists and create directory."""
        if os.path.isdir(pth) is False:
            os.mkdir(pth)
        return pth

    def log_info(self):
        """Log project-level details."""
        logging.info('ProjectName : {}'.format(self.ProjectName))
        logging.info('InputFolder : {}'.format(self.InputFolder))
        logging.info('OutputFolder: {}'.format(self.OutputFolder))
        logging.info('StartYear - End Year: {0}-{1}'.format(self.StartYear, self.EndYear))
        logging.info('Number of Months    : {}'.format(self.nmonths))

        if self.HistFlag.lower() in ['true', 't', 'yes', 'y', '1']:
            logging.info('Running: Historic Mode')

        else:
            logging.info('Running: Future Mode')
            try:
                logging.info('Historic Soil Moisture File: {}'.format(self.SavFile))
                logging.info('Historic Channel Storage File: {}'.format(self.ChStorageFile))
            except AttributeError:
                pass

        try:
            logging.info('Diagnostics will be performed using the data file: {}'.format(self.VICDataFile))
        except AttributeError:
            pass

    def update(self, args):
        """
        Overwrite configuration options.

        :@param args:   Dictionary of parameters, where the key is the parameter name
        """
        for k, v in args.items():
            if not hasattr(self, k):
                print('Warning: {} is not a valid parameter'.format(k))
            setattr(self, k, v)
