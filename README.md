[![DOI](https://zenodo.org/badge/88797535.svg)](https://zenodo.org/badge/latestdoi/88797535) [![Build Status](https://travis-ci.org/JGCRI/xanthos.svg?branch=master)](https://travis-ci.org/JGCRI/xanthos)

# Xanthos
Xanthos is an open-source hydrologic model, written in Python, designed to quantify and analyze global water availability. Xanthos simulates historical and future global water availability on a monthly time step at a spatial resolution of 0.5 geographic degrees. Xanthos was designed to be extensible and used by scientists that study global water supply and work with the Global Change Analysis Model (GCAM). Xanthos uses a user-defined configuration file to specify model inputs, outputs and parameters. Xanthos has been tested using actual global data sets and the model is able to provide historical observations and future estimates of renewable freshwater resources in the form of total runoff, average streamflow, potential evapotranspiration, actual evapotranspiration, accessible water, hydropower potential, and more.

# Contact Us
For questions, technical supporting and user contribution, please contact:

Chris Vernon <chris.vernon@pnnl.gov>

Zarrar Khan <zarrar.khan.pnnl.gov>

# Xanthos 2 - Upgrades
With the ability to simulate historical and future global water availability on a monthly time step at a spatial resolution of 0.5 geographic degrees, Xanthos version 1.0 provided a solid foundation for continued advancements in global water dynamics science.  The goal of Xanthos version 2 was to build upon previous investments by creating an accessible computing environment where core components of the model (potential evapotranspiration (PET), runoff generation, and river routing) could be interchanged or added to without having to start from scratch.  Xanthos 2 utilizes a component-style architecture which enables researchers to quickly incorporate and test cutting-edge research in a stable modeling environment prebuilt with a diagnostics module.  Major advancements for Xanthos 2.0 were also achieved by creating a more robust default configuration for the model that is now available to the scientific community.  These advancements include the addition of:  the Penman-Monteith PET module to capture the impacts of evolving land cover, the ABCD water balance module to account for groundwater recharge and discharge in runoff projections, improved water velocity considerations for the Modified River Transport Model (MRTM) routing module, a built-in differential evolution optimization module to calibrate ABCD parameters to modeled global runoff, and hydropower production assessment and potential capacity modules.  The figure below demonstrates the optimization module’s ability to calibrate Xanthos 2 runoff to the complex Variable Infiltration Capacity (VIC) model runoff projections when forced by the same climate data. Xanthos can be calibrated against other land surface models and Earth system models.

![Xanthos to VIC](https://github.com/JGCRI/xanthos/blob/master/docs/xanthos2_to_vic_watch_basins.png)
Figure:  Xanthos 2.0 performance when calibrated to the VIC model forced by WATCH observational climate data.  Each point represents the mean annual runoff for each of the 235 river basins in GCAM.

# Get Started with Xanthos
Set up Xanthos using the following steps:
1.  Clone Xanthos into your desired location `git clone https://github.com/JGCRI/xanthos.git`.
2.  From the directory you cloned Xanthos into run `python setup.py install`.
3.  Download the example data using the following in a Python prompt:
    ```python
    import xanthos
    
    # the directory that you want to download and extract the example data to
    data_dir = "<my data download location>"
    
    # download and unzip the package data to your local machine
    xanthos.get_package_data(data_dir)```
4.  Setup your configuration file (.ini).  Examples are located in the "example" directory.  Be sure to change the root directory to the directory that holds your data (use the `xanthos/example` directory as an example).
5. If running Xanthos from an IDE:  Be sure to include the path to your config file.  See the "xanthos/example/example.py" script as a reference.
6. If running Xanthos from terminal:  Run model.py found in xanthos/xanthos/model.py passing the full path to the config file as the only argument. (e.g., `python model.py <dirpath>/config.ini`).

# Setting up a Xanthos run
A detailed Wiki set up to describe how to set up a Xanthos run can be viewed here:  https://github.com/JGCRI/xanthos/wiki/Tutorial-1:--Setting-up-a-Xanthos-run

# Available Modules
A detailed Wiki set up to describe available modules, as well as their associated configuration settings, can be viewed here: https://github.com/JGCRI/xanthos/wiki/Available-modules
