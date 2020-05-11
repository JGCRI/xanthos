# Available Modules
This documentation details the parameters associated with each module available for use in Xanthos and how to represent them in the configuration file.

## Potential Evapotranspiration (PET)
Config tag:
```ini
[PET]
```

| Variable | Description | Required |
| -------- | ----------- | -------- |
| pet_module | The reference name of the PET model to be used.  Currently supported:  pm (for Panman-Monteith), hargreaves (for Hargreaves), and none (for cases where the user will provide their own PET input file). | True |
| pet_file | Only used if `pet_module` is set to `none`.  This will be a full path file with extension to the input PET data in mm/month.  Array must be 2-D, where (gridcell_idx, month_idx).  File must be saved as a NumPy array file (.npy). | False |


### Penman-Monetith
Config tag:
```ini
[penman-monteith]
```

| Variable | Description | Required |
| ------- | ------- | ------- |
| pet_dir | The directory name of the PET module chosen where its input files are stored | True |
| pm_tas | File name and extension to input mean monthly daily mean surface air temperature file in degree Celsius.  Array must be 2-D, where (gridcell_idx, month_idx).  File must be saved as a NumPy array file (.npy). | True |
| pm_tmin | File name and extension to input mean monthly daily minimum surface air temperature file in degree Celsius.  Array must be 2-D, where (gridcell_idx, month_idx).  File must be saved as a NumPy array file (.npy). | True |
| pm_rhs | File name and extension to input mean monthly relative humidity as a percentage.  Array must be 2-D, where (gridcell_idx, month_idx).   File must be saved as a NumPy array file (.npy). | True |
| pm_rlds | File name and extension to input mean monthly surface downwelling longwave radiation in W m-2.  Array must be 2-D, where (gridcell_idx, month_idx).  File must be saved as a NumPy array file (.npy). | True |
| pm_rsds | File name and extension to input mean monthly surface downwelling shortwave radiation in W m-2.  Array must be 2-D, where (gridcell_idx, month_idx).  File must be saved as a NumPy array file (.npy). | True |
| pm_wind | File name and extension to input mean monthly wind speed in m/s.  Array must be 2-D, where (gridcell_idx, month_idx).  File must be saved as a NumPy array file (.npy). | True |
| pm_lct | File name and extension to input land cover data fraction per 67,420 grid cells.  Array must be 3-D, where (gridcell_idx, landclass_idx, year_idx).  File must be saved as a NumPy array file (.npy). | True |
| pm_nlcs | The number of land cover classes in the input data as an integer. | True |
| pm_water_idx | The number of the index (when starting at 0) of water in the land cover data as an integer. | True |
| pm_snow_idx | The number of the index (when starting at 0) of snow in the land cover data as an integer. | True |
| pm_lc_years | Comma separated years (e.g., 2005, 2010, 2015) represented in the input land cover data. | True |

### Hargreaves-Samani
Config tag:
```ini
[hargreaves-samani]
```

| Variable | Description | Required |
| ------- | ------- | ------- |
| pet_dir | The directory name of the PET module chosen where its input files are stored | True |
| hs_tas | File name and extension to input mean monthly daily mean surface air temperature file in degree Celsius.  Array must be 2-D, where (gridcell_idx, month_idx).  File must be saved as a NumPy array file (.npy). | True |
| hs_tmin | File name and extension to input mean monthly daily minimum surface air temperature file in degree Celsius.  Array must be 2-D, where (gridcell_idx, month_idx).  File must be saved as a NumPy array file (.npy). | True |
| hs_tmax | File name and extension to input mean monthly daily maximum surface air temperature file in degree Celsius.  Array must be 2-D, where (gridcell_idx, month_idx).  File must be saved as a NumPy array file (.npy). | True |

### Hargreaves
Config tag:
```ini
[hargreaves]
```

| Variable | Description | Required |
| ------- | ------- | ------- |
| pet_dir | The directory name of the PET module chosen where its input files are stored | True |
| TemperatureFile | File name and extension to input mean monthly daily mean surface air temperature file in degree Celsius.  Array must be 2-D, where (gridcell_idx, month_idx).  File can be a NumPy array file (.npy) or a NetCDF Classic file (.nc). | True |
| TempVarName | If file is a NetCDF file, this is the variable name of temperature. | False |
| DailyTemperatureRangeFile | File name and extension to input mean monthly daily temperature range file in degree Celsius.  Array must be 2-D, where (gridcell_idx, month_idx).  File can be a NumPy array file (.npy) or a NetCDF Classic file (.nc). | True |
| DTRVarName | If file is a NetCDF file, this is the variable name of DTR. | False |


## Runoff
Config tag:
```ini
[Runoff]
```

| Variable | Description | Required |
| -------- | ----------- | -------- |
| runoff_module | The name of the runoff model to be used.  Currently supported:  abcd (for ABCD), gwam (for GWAM), and none (for cases where the user will provide their own runoff input file). | True |
| runoff_dir | The directory name of the runoff module chosen where its input files are stored | True |
| runoff_file | Only used if `runoff_module` is set to `none`.  This will be a full path file with extension to the input runoff data in mm/month.  Array must be 2-D, where (gridcell_idx, month_idx).  File must be saved as a NumPy array file (.npy). | False |

### ABCD
Config tag:
```ini
[abcd]
```

| Variable | Description | Required |
| -------- | ----------- | -------- |
| runoff_dir | The directory name of the runoff module chosen where its input files are stored | True |
| calib_file | File name with extension of calibrated ABCDM parameters for each basin.  Array must be 2-D, where (basin_idx, param_idx).  File must be saved as a NumPy array file (.npy). | True |
| runoff_spinup | The number of months from the start of the data that will be used to spin-up the data as an integer. | True |
| jobs | The number of jobs to use when running basins parallel (-2, all but one core; -1, all cores; 8, 8 cores). | False |
| TempMinFile | File name and extension to input mean monthly daily minimum surface air temperature file in degree Celsius.  File can be a NumPy array file (.npy) or a NetCDF Classic file (.nc). | True |
| TempMinVarName | If file is a NetCDF file, this is the variable name of temperature. | False |
| PrecipitationFile | File name and extension to input mean monthly preciptation file in mm/month.  Array must be 2-D, where (gridcell_idx, month_idx).  File can be a NumPy array file (.npy) or a NetCDF Classic file (.nc). | True |
| PrecipVarName | If file is a NetCDF file, this is the variable name of precipitation. | False |

### GWAM
Config tag:
```ini
[gwam]
```

| Variable | Description | Required |
| -------- | ----------- | -------- |
| runoff_dir | The directory name of the runoff module chosen where its input files are stored | True |
| PrecipitationFile | File name and extension to input mean monthly precipitation file in mm/month.  Array must be 2-D, where (gridcell_idx, month_idx).  File can be a NumPy array file (.npy) or a NetCDF Classic file (.nc). | True |
| PrecipVarName | If file is a NetCDF file, this is the variable name of precipitation. | False |
| max_soil_moisture | File of maximum soil moisture (gridcells x 1) in mm/month. | True |
| lakes_msm | File of water bodies: assign MSM = 999, 306 x 2, Col 1 is the cell number in 67420. | True |
| addit_water_msm | File of additional water bodies: assign MSM = 999, 421 x 2, Col 1 is the cell number in 67420. | True |

## Routing
Config tag:
```ini
[Routing]
```

| Variable | Description | Required |
| -------- | ----------- | -------- |
| routing_module | The reference name of the routing model to be used.  Currently supported:  mrtm (for MRTM), and none (for cases where the user does not which to run the routing module). | True |

### Modified River Transport Model (MRTM)
Config tag:
```ini
[mrtm]
```

| Variable | Description | Required |
| -------- | ----------- | -------- |
| routing_dir | The directory name of the routing module chosen where its input files are stored | True |
| routing_spinup | The number of months from the start of the data that will be used to spin-up the data as an integer. | True |

# Optional Modules

## Diagnostics
Config tag:
```ini
[Diagnostics]
```

This section was built to compare runoff outputs against that of other models.

| Variable | Description | Required |
| -------- | ----------- | -------- |
| VICDataFile | File name with extension of the input yearly mean runoff in km3/yr from VIC.  Array must be 2-D, where (gridcell_idx, year_idx). File can be a NetCDF Classic file or a CSV file. | True |
| WBMDataFile | File name with extension of the input yearly mean runoff in km3/yr from the WBM model.  Array must be 2-D, where (gridcell_idx, year_idx). File can be a NetCDF Classic file or a CSV file. | True |
| WBMCDataFile | File name with extension of the input yearly mean runoff in km3/yr from the WBMc model.  Array must be 2-D, where (gridcell_idx, year_idx). File can be a NetCDF Classic file or a CSV file. | True |
| UNHDataFile | File name with extension of the input yearly mean runoff in km3/yr from the UNH-GRDC model.  Array must be 2-D, where (gridcell_idx, year_idx). File can be a NetCDF Classic file or a CSV file. | True |
| Scale | Integer to define the level of output:  0 = all; 1 = Basin; 2 = Country; 3 = Region | True |

## Time Series Plotting
Config tag:
```ini
[TimeSeriesPlot]
```

Create time series plots of runoff.

| Variable | Description | Required |
| -------- | ----------- | -------- |
| MapID | Define the level of plotting:  999 = plot all; 0 = only global total; 1 = only by basin id | True |

## Drought
Config tag:
```ini
[Drought]
```

| Variable | Description | Required |
| -------- | ----------- | -------- |
| drought_var | Which output variable to use, either 'q' (runoff) or 'soilmoisture' (soil moisture). | True |
| drought_thresholds | Path to NumPy array of quantile-based drought thresholds. If provided, Severity, Intensity, and Duration statistics are calculated. If not provided, the following three parameters are required, and thresholds are calculated. | False |
| threshold_start_year | First year for which to calculate drought thresholds. | False |
| threshold_end_year | Last year for which to calculate drought thresholds. | False |
| threshold_nper | Number of periods to calculate thresholds for.  Generally either 1 (single threshold for all periods), or 12 (thresholds by month) | False |

## Accessible Water
Config tag:
```ini
[AccessibleWater]
```

| Variable | Description | Required |
| -------- | ----------- | -------- |
| ResCapacityFile | File name with extension of the reservior storage capacity file as a CSV. | True |
| BfiFile | File name with extension of the baseflow index (BFI) file as a CSV. | True |
| HistEndYear | Four digit year YYYY for the historical end year. | True |
| GCAM_StartYear | The start year YYYY of the Global Change Assessment Model (GCAM) to provide data for. | True |
| GCAM_EndYear | The end year YYYY of GCAM to provide data for. | True |
| GCAM_YearStep | The time step in years as an integer to process | True |
| MovingMeanWindow | The window size of the rolling mean as an integer. | True |
| Env_FlowPercent | A fractional decimal from 0 to 1.0 to calculate environmental flow requirements per basin. | True |

## Hydropower Potential
Config tag:
```ini
[HydropowerPotential]
```

| Variable | Description | Required |
| -------- | ----------- | -------- |
| hpot_start_date | Start date in "M/YYYY" format. | True |
| q_ex | Quantile of monthly flow above which additional power is unavailable; values from 0.0 to 1.0. | True |
| ef | Plant efficiency; values from 0.0 to 1.0. | True |

## Hydropower Actual
Config tag:
```ini
[HydropowerActual]
```

| Variable | Description | Required |
| -------- | ----------- | -------- |
| hact_start_date | Start date in "M/YYYY" format. | True |


## Calibration for the ABCD model
Config tag:
```ini
[Calibrate]
```

| Variable | Description | Required |
| -------- | ----------- | -------- |
| set_calibrate | Integer; 0 = calibrate runoff against observed runoff for the ABCD model; 1 = calibrate ABCD runoff parameters against observed streamflow using outputs from the MRTM model. | True |
| observed | Full path with filename and extension of observed runoff or streamflow dataset; Must be in `basin,year,month,value` CSV format with header. | True |
| obs_unit | Observed runoff or streamflow units (runoff options:  km3_per_mth, mm_per_mth; streamflow options: m3_per_sec). | True |
| calib_out_dir | Full path with to directory where output files for KGE and the ABCDM parameters will be saved. | True |
| calibration_basins | Specific basins to calibrate, specified as a comma-separated list (e.g. 1,5-9,12). | False |
