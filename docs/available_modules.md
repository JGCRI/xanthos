# Available Modules
This documentation details the parameters associated with each module available for use in Xanthos and how to represent them in the configuration file.

## Potential Evapotranspiration (PET)
| Variable | Description | Required |
| -------- | ----------- | -------- |
| pet_module | The reference name of the PET model to be used.  Currently supported:  pm (for Panman-Monteith), hargreaves (for Hargreaves), and none (for cases where the user will provide their own PET input file). | True |
| pet_file | Only used if `pet_module` is set to `none`.  This will be a full path file with extension to the input PET data in mm/month.  Array must be 2-D, where (gridcell_idx, month_idx).  File must be saved as a NumPy array file (.npy). | False |


### Penman-Monetith
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


### Hargreaves
| Variable | Description | Required |
| ------- | ------- | ------- |
| pet_dir | The directory name of the PET module chosen where its input files are stored | True |
| TemperatureFile | File name and extension to input mean monthly daily mean surface air temperature file in degree Celsius.  Array must be 2-D, where (gridcell_idx, month_idx).  File can be a NumPy array file (.npy) or a NetCDF Classic file (.nc). | True |
| TempVarName | If file is a NetCDF file, this is the variable name of temperature. | False |
| DailyTemperatureRangeFile | File name and extension to input mean monthly daily temperature range file in degree Celsius.  Array must be 2-D, where (gridcell_idx, month_idx).  File can be a NumPy array file (.npy) or a NetCDF Classic file (.nc). | True |
| DTRVarName | If file is a NetCDF file, this is the variable name of DTR. | False |


## Runoff
| Variable | Description | Required |
| -------- | ----------- | -------- |
| runoff_module | The name of the runoff model to be used.  Currently supported:  abcd (for ABCD), gwam (for GWAM), and none (for cases where the user will provide their own runoff input file). | True |
| runoff_dir | The directory name of the runoff module chosen where its input files are stored | True |
| runoff_file | Only used if `runoff_module` is set to `none`.  This will be a full path file with extension to the input runoff data in mm/month.  Array must be 2-D, where (gridcell_idx, month_idx).  File must be saved as a NumPy array file (.npy). | False |

### ABCD
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
| Variable | Description | Required |
| -------- | ----------- | -------- |
| runoff_dir | The directory name of the runoff module chosen where its input files are stored | True |
| PrecipitationFile | File name and extension to input mean monthly preciptation file in mm/month.  Array must be 2-D, where (gridcell_idx, month_idx).  File can be a NumPy array file (.npy) or a NetCDF Classic file (.nc). | True |
| PrecipVarName | If file is a NetCDF file, this is the variable name of precipitation. | False |


## Routing
| Variable | Description | Required |
| -------- | ----------- | -------- |
| routing_module | The reference name of the routing model to be used.  Currently supported:  mrtm (for MRTM), and none (for cases where the user does not which to run the routing module). | True |

### Modified River Transport Model (MRTM)
| Variable | Description | Required |
| -------- | ----------- | -------- |
| routing_dir | The directory name of the routing module chosen where its input files are stored | True |
| routing_spinup | The number of months from the start of the data that will be used to spin-up the data as an integer. | True |
