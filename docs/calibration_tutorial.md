# Calibrate ABCD Model
These tutorials demonstrate how to calibrate the ABCD model independently in combination with the streamflow.

## Calibrate the ABCD runoff module against observed runoff data

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

#### STEP 3:  Run it!
This may take a while depending on how many function evaluations it takes to optimize each basin.

## Calibrate the ABCD runoff module coupled with the MRTM routing module against observed streamflow

#### STEP 1:  Provide the observed streamflow data in the required format
The format must be CSV and have a header row of `basin,year,month,value`.
Each row must have the corresponding data.

#### STEP 2:  Set up your configuration file
Turn on the calibration setting:  `Calibration = 1`.
Add the `[Calibrate]` section to your config file if it is not already present.  The section should contain the following:
```
[Calibrate]
set_calibrate = <set to 1 for calibrating to observed streamflow>
observed = <full path with file name and extension to your observational data>
obs_unit = <the unit of the input data; m3_per_sec for streamflow>
```

#### STEP 3:  Run it!
This may take a while depending on how many function evaluations it takes to optimize each basin.
