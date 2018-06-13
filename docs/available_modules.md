# Available Modules
This documentation details the parameters associated with each module available for use in Xanthos.

## Potential Evapotranspiration (PET)



### Hargreaves
| Variable | Description | Required |
| ------- | ------- | ------- |
| TemperatureFile | Full path with file name and extension to input mean temperature file | True |

### Penman-Monetith
| Variable | Description | Required |
| ------- | ------- | ------- |
| TemperatureFile | Full path with file name and extension to input minimum temperature file | True |




#### STEP 1:  Provide the observed runoff data in the required format
The format must be CSV and have a header row of ```basin,year,month,value```.
Each row must have the corresponding data.

#### STEP 2:  Set up your configuration file
Turn on the calibration setting:  `Calibration = 1`.
Add the `[Calibrate]` section to your config file if it is not already present.  The section should contain the following:
```
[Calibrate]
set_calibrate = <set to 0 for calibrating to observed runoff>
observed = <full path with file name and extension to your observational data>
obs_unit = <the unit of the input data; either km3_per_mth or mm_per_mth for runoff>
```