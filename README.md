[![DOI](https://zenodo.org/badge/88797535.svg)](https://zenodo.org/badge/latestdoi/88797535)

# Xanthos
Xanthos is an open-source hydrologic model, written in Python, designed to quantify and analyze global water availability. Xanthos simulates historical and future global water availability on a monthly time step at a spatial resolution of 0.5 geographic degrees. Xanthos was designed to be extensible and used by scientists that study global water supply and work with the Global Change Analysis Model (GCAM). Xanthos uses a user-defined configuration file to specify model inputs, outputs and parameters. Xanthos has been tested using actual global data sets and the model is able to provide historical observations and future estimates of renewable freshwater resources in the form of total runoff, average streamflow, potential evapotranspiration, actual evapotranspiration, accessible water, hydropower potential, and more.

# Get Started with Xanthos
Set up Xanthos using the following steps:
1. Install Xanthos 

    * install released version from PyPi using:
    ```bash
    pip install xanthos
    ```
   
    * or install the latest version hosted on Github using:
    ```bash
    python -m pip install git+https://github.com/JGCRI/xanthos.git
    ```
2. Download the example data using the following in a Python prompt:
    ```python
    import xanthos
    
    # the directory that you want to download and extract the example data to
    data_dir = "<my data download location>"
    
    # download and unzip the package data to your local machine
    xanthos.get_package_data(data_dir)
    ```
3. Setup your configuration file (.ini).  Examples are located in the "example" directory that you just downloaded.  Be sure to change the following variables to represent the local path to your example data:  `RootDir`, `TempMinFile`, `PrecipitationFile`.
4. To run Xanthos:

    ```python
    import xanthos
    
    # the path and file name that my example configuration (.ini) file was downloaded to
    config_file = '<path to my example config file>/pm_abcd_mrtm.ini'
    
    # run Xanthos 
    xanthos.run_model(config_file)
    ```

# Setting up a Xanthos run
A detailed Wiki set up to describe how to set up a Xanthos run can be viewed here:  https://github.com/JGCRI/xanthos/wiki/Tutorial-1:--Setting-up-a-Xanthos-run

# Available Modules
A detailed Wiki set up to describe available modules, as well as their associated configuration settings, can be viewed here: https://github.com/JGCRI/xanthos/wiki/Available-modules

# Xanthos 2 - Upgrades
With the ability to simulate historical and future global water availability on a monthly time step at a spatial resolution of 0.5 geographic degrees, Xanthos version 1.0 provided a solid foundation for continued advancements in global water dynamics science.  The goal of Xanthos version 2 was to build upon previous investments by creating an accessible computing environment where core components of the model (potential evapotranspiration (PET), runoff generation, and river routing) could be interchanged or added to without having to start from scratch.  Xanthos 2 utilizes a component-style architecture which enables researchers to quickly incorporate and test cutting-edge research in a stable modeling environment prebuilt with a diagnostics module.  Major advancements for Xanthos 2.0 were also achieved by creating a more robust default configuration for the model that is now available to the scientific community.  These advancements include the addition of:  the Penman-Monteith PET module to capture the impacts of evolving land cover, the ABCD water balance module to account for groundwater recharge and discharge in runoff projections, improved water velocity considerations for the Modified River Transport Model (MRTM) routing module, a built-in differential evolution optimization module to calibrate ABCD parameters to modeled global runoff, and hydropower production assessment and potential capacity modules.  The figure below demonstrates the optimization module’s ability to calibrate Xanthos 2 runoff to the complex Variable Infiltration Capacity (VIC) model runoff projections when forced by the same climate data. Xanthos can be calibrated against other land surface models and Earth system models.

![Xanthos to VIC](docs/xanthos2_to_vic_watch_basins.png)

Figure:  Xanthos 2.0 performance when calibrated to the VIC model forced by WATCH observational climate data.  Each point represents the mean annual runoff for each of the 235 river basins in GCAM.

# Contact Us
For questions, technical supporting and user contribution, please contact:

Mengqi Zhao <mengqi.zhao@pnnl.gov>

Chris Vernon <chris.vernon@pnnl.gov>

# Citations

Vernon, C.R., Hejazi, M.I., Turner, S.W.D., Liu, Y., Braun, C.J., Li, X. and Link, R.P. (2019). A Global Hydrologic Framework to Accelerate Scientific Discovery.  Journal of Open Research Software,  7(1), p.1. DOI: https://doi.org/10.5334/jors.245.

# Related Publications

* Abeshu, G. W., Tian, F., Wild, T., Zhao, M., Turner, S., Chowdhury, A. F. M. K., Vernon, C. R., Hu, H., Zhuang, Y., Hejazi, M., and Li, H.-Y. (2023). Enhancing the representation of water management in global hydrological models, Geosci. Model Dev. Discuss. https://doi.org/10.5194/gmd-2023-12.
* Chowdhury, K., Wessel, J., Wild, T.B., Lamontagne, J., Kanyako, F. (2023). Exploring Sustainable Electricity System Development Pathways in South America’s MERCOSUR Sub-Region. Energy Strategy Reviews. https://doi.org/10.1016/j.esr.2023.101150.
* Khan, Z., Thompson, I., Vernon, C., Graham, N., Chen, M., Wild, T.B. (2023). A global 0.5° gridded multi-sector monthly water withdrawal dataset for 2015-2100 under alternative futures. Nature Scientific Data. https://doi.org/10.1038/s41597-023-02086-2.
* Kanyako, F., Baker, E., Lamontagne, J., Turner, S., Wild, T.B. (2022). Seasonality and Trade in Hydro-heavy Electricity Markets: The West Africa Power Pool (WAPP). Applied Energy. DOI: https://doi.org/10.1016/j.apenergy.2022.120214.
* Birnbaum, A.N., Lamontagne, J.R., Wild, T.B., Dolan, F., Yarlagadda, B. (2022). Drivers of Physical and Economic Water Scarcity in Latin America and the Caribbean. Earth’s Future. DOI: https://doi.org/10.1029/2022EF002764.
* Liu, Y., Hejazi, M., Li, H., Zhang, X., and Leng, G. (2018). A hydrological emulator for global applications – HE v1.0.0, Geosci. Model Dev., 11, 1077–1092, https://doi.org/10.5194/gmd-11-1077-2018.
* Turner, S.W., Ng, J.Y. and Galelli, S. (2017). Examining global electricity supply vulnerability to climate change using a high-fidelity hydropower dam model. Science of the Total Environment, 590: 663–675. DOI: https://doi.org/10.1016/j.scitotenv.2017.03.022.
* Turner, S.W., Hejazi, M., Kim, S.H., Clarke, L., Edmonds, J. (2017). Climate impacts on hydropower and consequences for global electricity supply investment needs. Energy, 141: 2081–2090. DOI: https://doi.org/10.1016/j.energy.2017.11.089.
