# Adding Xanthos components
These tutorials demonstrate how to add a component to Xanthos using the accessible water module as an example.


## STEP 1:  Ensure compatibility with what your model ingests/outputs with what Xanthos outputs/ingests
Functionality for converting to Xanthos IO spatiotemporal scale and units may have to be made in your module to ensure compatibility before integration.

## STEP 2:  Set up your configuration file
Create a new, unique configuration tag for your module that will be represented in the configuration.ini file and read in `xanthos/data_reader/ini_reader.py`.  This section in the configuration.ini file will contain the variables and their associated values in a key: value pair similar to the following:

```ini
[AccessibleWater]

# Reservoir capacity at basin level
ResCapacityFile             = total_reservoir_storage_capacity_BM3.csv

# Baseflow index (BFI) file
BfiFile                     = bfi_per_basin.csv

# HistEndYear: End year of historical data, e.g. 2005
HistEndYear                 = 2001
GCAM_StartYear              = 1971
GCAM_EndYear                = 2001
GCAM_YearStep               = 1

# A parameter for moving average
MovingMeanWindow            = 9

# Used to calculate Environmental Flow Requirements (EFR) per basin, for example, use 10% of historical mean
Env_FlowPercent             = 0.1
```

Create a variable and value under the `[Project]` tag that will be used to determine the name of your directory that will contain any inputs that your module requires:

```ini
# directory name where the accessible water input file directory is contained
AccWatDir                   = accessible_water
```

Create a variable and value under the `[Project]` tag that will be used to determine whether or not your module will be executed during runtime:

```ini
# calculate accessible water; Default is 0 for False, 1 for True
CalculateAccessibleWater    = 1
```

## STEP 3:  Prepare `xanthos/data_reader/ini_reader.py`

Create a new module component check in `COMPONENT MODULE CHECK` section:

```python
try:
    a = c['AccessibleWater']
    self.AccWatDir = os.path.join(self.InputFolder, p['AccWatDir'])
except KeyError:
    a = False
```

Create module run expectation (0 or 1) in `PROJECT LEVEL SETTINGS` section:

```python
self.CalculateAccessibleWater = int(p['CalculateAccessibleWater'])
```

Create conditional statement to load config settings in `OPTIONAL POST-PROCESSING MODULES` section:

```python
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
            raise ValidationException("Accessible water range of GCAM years are outside the range of years in climate data.")
```

## STEP 4:  Create a module for your code containing a `__init__.py` file and the other code to run your module.  The `__init__.py` file can be empty; its presence lets Python know that your directory is a Python module so that it may be called accordingly.

## STEP 5:  Create a method in the `Components` class in the `xanthos\components.py` file that calls your module:

```python
def accessible_water(self):
    """
    Run accessible water module

    :@param s:          Settings object from configuration
    :@param data:       Data object from imported ancillary data
    :@param Q:          Monthly runoff NumPy array in mm/month; [ngrids, nmonths]
    """
    if self.s.CalculateAccessibleWater:
        AccessibleWater(self.s, self.data, self.Q)
```

## STEP 5:  Add your method call to an existing or new configuration function in `xanthos\configurations.py` after the simulation calls:

```python
def pm_abcd_mrtm(config):
    """
    Model configuration for the following:

    PET:                Penman-Monteith
    RUNOFF:             ABCD
    ROUTING:            Modified River Transport Model (MRTM)

    Spin-up is built in to the ABCD model and does not need to be calculated separately.

    :param config:      Configuration object generated from user-defined config.ini file
    :return:            Object containing all return values.
    """
    # instantiate hydro class
    c = Components(config)

    # run model
    c.simulation(pet=True,
                 pet_step=None,
                 runoff=True,
                 runoff_step=None,
                 routing=True,
                 routing_num_steps=config.nmonths,
                 routing_step='month',
                 notify='Simulation')

    # accessible water module
    c.accessible_water()

    # hydropower potential
    c.hydropower_potential()

    # hydropower actual
    c.hydropower_actual()

    # diagnostics
    c.diagnostics()

    # output simulation data
    c.output_simulation()

    # aggregate outputs
    c.aggregate_outputs()

    # create time series plots
    c.plots()

    return c
```

## STEP 6: Run it!
