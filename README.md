[![license](https://img.shields.io/pypi/l/pyunfold.svg 'license')](https://github.com/JGCRI/xanthos/blob/master/LICENSE)

# Xanthos
Xanthos is an open-source hydrologic model, written in Python, designed to quantify and analyze global water availability. Xanthos simulates historical and future global water availability on a monthly time step at a spatial resolution of 0.5 geographic degrees. Xanthos was designed to be extensible and used by scientists that study global water supply and work with the Global Change Assessment Model (GCAM). Xanthos uses a user-defined configuration file to specify model inputs, outputs and parameters. Xanthos has been tested using actual global data sets and the model is able to provide historical observations and future estimates of renewable freshwater resources in the form of total runoff, average streamflow, potential evapotranspiration, actual evapotranspiration, accessible water, hydropower potential, and more.

# Contact Us
For questions, technical supporting and user contribution, please contact:

Vernon, Chris R <Chris.Vernon@pnnl.gov>

Li, Xinya <Xinya.Li@pnnl.gov>

Link, Robert P <Robert.Link@pnnl.gov>

Hejazi, Mohamad I <Mohamad.Hejazi@pnnl.gov>


# Notice
This repository uses the Git Large File Storage (LFS) extension (see https://git-lfs.github.com/ for details).  Please run the following command before cloning if you do not already have Git LFS installed:
`git lfs install`.

# Xanthos 2.0 - Upgrades
With the ability to simulate historical and future global water availability on a monthly time step at a spatial resolution of 0.5 geographic degrees, Xanthos version 1.0 provided a solid foundation for continued advancements in global water dynamics science.  The goal of Xanthos version 2.0 was to build upon previous investments by creating an accessible computing environment where core components of the model (potential evapotranspiration (PET), runoff generation, and river routing) could be interchanged or added to without having to start from scratch.  Xanthos 2.0 utilizes a component-style architecture which enables researchers to quickly incorporate and test cutting-edge research into a stable modeling environment prebuilt with a diagnostics module.  Major advancements for Xanthos 2.0 were also achieved by creating a more robust default configuration for the model that is available to the scientific community.  These advancements include the addition of:  the Penman-Monteith PET module to capture the impacts of evolving land cover, the ABCD water balance module to account for groundwater recharge and discharge in runoff projections, improved water velocity considerations for the routing module, a built-in differential evolution optimization module to calibrate ABCD parameters to modeled global runoff, hydropower production assessment and potential capacity modules.  The figure below demonstrates the optimization module’s ability to calibrate Xanthos 2.0 runoff to the complex Variable Infiltration Capacity (VIC) model runoff projections when forced by the same climate data. Xanthos can be calibrated against other land surface models and Earth system models.

# Get Started (Xanthos 2.0)
Set up Xanthos using the following steps:
1.  This repository uses the Git Large File Storage (LFS) extension (see https://git-lfs.github.com/ for details).  Please run the following command before cloning if you do not already have Git LFS installed:
`git lfs install`.
2.  Clone Xanthos into your desired location `git clone https://github.com/JGCRI/xanthos.git`.  Some Windows users have had better luck with `git lfs clone https://github.com/JGCRI/xanthos.git`
3.  From the directory you cloned Xanthos into run `python setup.py install` .  This will install Xanthos as a Python package on your machine and install of the needed dependencies.
4.  Setup your configuration file (.ini).  Examples are located in the "example" directory.  Be sure to change the root directory to the directory that holds your data (use the 'xanthos/example' directory as an example).
5. If running Xanthos from an IDE:  Be sure to include the path to your config file.  See the "xanthos/example/example.py" script as a reference.
6. If running Xanthos from terminal:  Run model.py found in xanthos/xanthos/model.py passing the full path to the config file as the only argument. (e.g., `python model.py <dirpath>/config.ini`).

## Set up configuration file
An example configuration file set to the default Penman-Monteith/ABCD/MRTM configuration is included here: `examples/pm_abcd_mrtm.ini`.  The configuration parameters associated with each PET, runoff, and routing model option can be reviewed here: [Available modules](https://github.com/JGCRI/xanthos/blob/version-two-setup/docs/available_modules.md).  To get started with the default configuration, set the following configuration variables to your local preferences:
```ini
In [PROJECTS]...
RootDir = <the directory holding your input and output directories>

In [Runoff][abcd]...
TempMinFile = <the directory holding your example directory>/example/input/climate/tasmin_watch_monthly_degc_1991_2001.npy
PrecipitationFile = <the directory holding your example directory>/example/input/climate/pr_gpcc_watch_monthly_mmpermth_1991_2001.npy

In [Calibrate]...
observed =  <the directory holding your example directory>/example/input/calibration/vic_watch_basin_km3_1971_2001_monthly.csv
calib_out_dir = <a local directory of your choosing>
```

## Set up example.py
Change path to configuration file to your location and run!

```python
from xanthos.model import Xanthos

def run(ini):

    # instantiate model
    xth = Xanthos(ini)

    # run intended configuration
    xth.execute()

    return xth


if __name__ == "__main__":

    # full path to parameterized config file
    ini = '<the directory path to your example directory>/pm_abcd_mrtm.ini'

    # run the model
    xth = run(ini)
```
