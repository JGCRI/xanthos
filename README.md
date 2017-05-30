# xanthos
Xanthos is a Python model designed to quantify and analyse global water availability historically and in the future at 0.5° × 0.5° spatial resolution and a monthly time step.  Its performance and functionality was tested through real-world applications. It is open-source, extensible and accessible for researchers who work on long-term climate data for studies of global water supply, and the Global Change Assessment Model (GCAM). This model integrates inherent global gridded data maps, I/O modules, hydrologic processes and diagnostics modules parameterized by a user-defined configuration file.

# Notice
This repository uses the Git Large File Storage (LFS) extension (see https://git-lfs.github.com/ for details).  Please run the following command before cloning if you do not already have Git LFS installed:
`git lfs install`

# Introduction
Xanthos is a model for calculating global water availability, developed at the Joint Global Change Research Institute of the Pacific Northwest National laboratory, USA. The initial objective of this model is to quantify changes in future freshwater availability under various climate change regimes, and to serve as the freshwater supply component of the Global Change Assessment Model (GCAM) [1-3]. The model has been used in previous publications to explore different climate and socioeconomic scenarios over the 21st century and assess their implications on water scarcity regionally and globally [4-5].

The core algorithm of Xanthos is the hydrology model, which includes modules for calculating potential evapotranspiration (PET), runoff generation, and stream routing – shown in Figure 1. The underlying equations and algorithms for the hydrology model were described in details in Hejazi et al [4], and Zhou et al. [6]. The model requires gridded monthly temperature and precipitation, a maximum soil water storage capacity map (assumed static over time in this study) and calculates gridded monthly runoff, PET, actual evapotranspiration, water storage in the soil column, channel water storage and average channel flow. The PET module uses the Hargreaves method [5], and requires monthly temperature data (average temperature and average daily temperature range in Celsius). The stream routing module includes a cell-to-cell river routing scheme adopted from the modified river transport model (RTM), which uses a linear advection formula [6]. The aggregation module of Xanthos provides estimates of maximum naturally-available water fluxes (mm) or volumes (billion m3) monthly or annually (i.e., total monthly runoff) for each river basin. Xanthos also offers the user an option to obtain the annual volumes (billion m3/year) of renewable water resources that are accessible for human use in each basin [7], which is calculated by the accessible-water module.

Xanthos was designed for independent applications and it can also serve as the water supply component in GCAM. The primary goals are as follows:
1.	Be distributed as an open-source model that can be utilized and extended by researchers who study global hydrology
2.	Estimate gridded (0.5 degree resolution) monthly global results of runoff, average channel flow, actual evapotranspiration, and water storage in the soil column
3.	Employ a standardized configuration and simplified input structure that is easy for the user to update
4.	Direct diagnostics on the runoff results.

To achieve the above goals, Python was selected as the programming language because of its easy-to-use and easy-to-install-with-libraries and it flexibility.

![alt-text](https://github.com/JGCRI/xanthos/blob/master/docs/workflow.png)
Figure 1: Main inputs and outputs of the hydrology model (Runoff Generation & Stream Routing) which is the core module of Xanthos

# Implementation and architecture
The grid used for global map is a 0.5° resolution, which results in an original matrix of 360 × 720 to present the data. The total number of valid grids (land grids) is 67420. Those 67420 cells define a “gridded” map according to the coordinates listed in a separated coordinate file. Certain commonly used global data maps such as IDs of basins/countries/regions are harmonized into the gridded format.

Xanthos consists of four modules. The flow chart is illustrated in Figure 2. As described in the introduction, PET calculation, runoff generation, and river routing are the core components of the hydrology model. Xanthos requires comprehensive input data files, which is divided into two classes. The climate data files are the forcing input to the model. The global data maps are the static input that are kept the same for each simulation. The direct outputs from the hydrology model are monthly gridded runoff, average channel flow, potential evapotranspiration, actual evapotranspiration, channel storage, and water storage in soil column. The user can choose if they need: 1) aggregated runoff results by basin, country or region; 2) time series plots of runoff and average channel flow by basin, country or region; 3) preform comparison with other models (diagnostics); and 4) the calculation of accessible water.

Corresponding to the flowchart in Figure 2, Figure 3 shows the overall architecture of Xanthos.  In the source code, “GCAM_Hydro.py” defines the main function that integrates all the modules shown in Figure 2. “Test.py” is the executable Python file connecting “config.ini” and “GCAM_Hydro.py”. Thus, a simulation can be executed with a simple command:

`$ python test.py`

To simplify the how-to-run process for the users, the full model also has separated the input package from the code package. The input package provides:
1.	Climate data files, which constitute the major user controlled inputs, include three netcdf files for precipitation, temperature and daily maximum temperature range. The data is stored in a matrix format with the dimensions of 67420 rows and number-of-months columns.
2.	Harmonized gridded global data maps and other parameter data files required, which can be modified or replaced by the user easily.
3.	Example data files of other models for comparison in diagnostics.
4.	Example data files for accessible water.

The details (data source, format, related pre-processing, etc.) of the input files are listed in a document included in the input folder.

The interface between the user and Xanthos is a configuration file, i.e., “config.ini”. One simulation requires a unique configuration file. Inside this configuration file, sections are organized and each section contains a few of name-value pairs for model settings (Figure 4). Six sections are defined in the configuration file:
1.	Project (Required): defines the major settings such as the paths of input and output folders, the output formatting, options for aggregated runoff results, an option for diagnostics, an option for time series plot, and an option for calculation of accessible water
2.	Climate (Required): defines the historical/future mode, paths of input climate forcing data - precipitation, temperature and daily temperature range files, and an option for initialization
3.	GriddedMap (Required): defines the required data maps, such as coordinates of the grids, maximum soil moisture of each gird, ID of basin/country/region for each gird
4.	Diagnostics (Required only if “PerformDiagnostics = 1” in section 1): defines the paths of the data files from other models, and the options for comparison at basin, country or region scale.
5.	TimeSeriesPlot (Required only if “CreateTimeSeriesPlot = 1” in section 1): defines the scale for plots and the selected basin/county/region to be used for generating the plots.
6.	AccessibleWater (Required only if “CalculateAccessibleWater = 1” in section 1): defines all the related parameters to the calculation of accessible water, such as the total reservoir capacity at basin level, and baseflow index file.

Section 4, 5 and 6 are optional responding to the flags defined in Section 1. An output folder will be automatically generated at the path defined in the configuration file.

Xanthos can be executed in either historical or future mode. For the historical mode, the climate forcing data is taken from a historical time period (e.g. data before year 2005).  As for the future mode, initial soil moisture and channel storage conditions can be taken either from a historical case calculated by Xanthos or from other models, and the climate forcing data is for future years (e.g. 2006-2100). Xanthos allows the spin-up (model initiation by using the results from the first few years) phase by setting the years for the spin-up time.

Xanthos was developed using Python1 (version 2.7) and related libraries. NumPy2 and SciPy3 are the fundamental packages for scientific computing data processing. Matplotlib4 is a library used for 2D plotting of timer series plots and scatter plots from diagnostics. Pandas5 was used in accessible water module for data analysis. The results are able to be saved into two formats: CSV (comma-separated values) and NetCDF. NetCDF6 is a self-describing, machine-independent data formats for array-oriented scientific data which requires much smaller storage space. By default, the input climate data files should be stored as NetCDF files, but they can also be stored as mat (MATLAB7 formatted) files if needed. The mat files can be read and write by input/output modules of SciPy.

![alt-text](https://github.com/JGCRI/xanthos/blob/master/docs/flowchart.png)
Figure 2: Flowchart of Xanthos

![alt-text](https://github.com/JGCRI/xanthos/blob/master/docs/structure.png)

Figure 3: The architecture of Xanthos in application

![alt-text](https://github.com/JGCRI/xanthos/blob/master/docs/config.png)
Figure 4: Example of a configuration file

# Quality control
This model was intensively tested by cases of different climate forcing on Linux and Windows.  Two example cases were included with the input climate data files, to help the user get familiar with the features and functions of Xanthos:
  A. Case for a historical run, and with options for aggregation, diagnostics, accessible water and time series plots, turned on.
  B. Case for a future run, with the initial channel storage and soil moisture data taken from the ending values of case A, and with basic options.

A direct indicator of performance is the automatically generated log file. The log file lists model settings, the processing steps, and CPU cost and warnings if applicable. The diagnostics module provides another option for the user to examine the quality of the model results of runoff. It compares the runoff results by basin scale, country scale or region scale to data sets of other models, such as the Variable Infiltration Capacity (VIC) model [8], Water Balance Model (WBM) and Water Balance Model Composite (WBMc) [9], UNH/GRDC [10-11]. Figure 5 is a scatter plot that shows an example of the plots created by the diagnostics module. If the scatters are blew the “x=y” line, the model overestimates the runoff compared to four models (e.g., VIC, WBM, WBMc and UNH/GRDC), and if they are above then the model underestimates. The gridded streamflow outputs are not regulated streamflow results. They are not available for direct comparison with other researches (for example, Dai’s results [12]), and hence are not included in diagnostics.

![alt-text](https://github.com/JGCRI/xanthos/blob/master/docs/diag_basin.png)
Figure 5: Example of diagnostics plot for comparison of averaged annual runoff with results from other models at basin scale.

# References

[1]   Edmonds, J., and Reilly J. M., 1985. Global Energy: Assessing the Future. Oxford University Press, New York, pp.317.

[2]   Edmonds, J., Wise, M., Pitcher, H., Richels, R., Wigley, T. and Maccracken, C., 1997. An integrated assessment of climate change and the accelerated introduction of advanced energy technologies-an application of MiniCAM 1.0. Mitigation and adaptation strategies for global change 1(4):.311-339. DOI: http://dx.doi.org/10.1023/B:MITI.0000027386.34214.60

[3]   Kim, S.H., Edmonds, J., Lurz, J., Smith, S.J. and Wise, M., 2006. The objECTS framework for integrated assessment: Hybrid modeling of transportation. The Energy Journal Special Issue #2: 51-80.

[4]   Hejazi, M.I., Edmonds, J., Clarke, L., Kyle, P., Davies, E., Chaturvedi, V., Wise, M., Patel, P., Eom, J. and Calvin, K., 2014. Integrated assessment of global water scarcity over the 21st century under multiple climate change mitigation policies. Hydrology and Earth System Sciences 18: 2859-2883. DOI: http://dx.doi.org/10.5194/hess-18-2859-2014

[5]   Hejazi, M.I., Edmonds, J., Chaturvedi, V., Davies, E. and Eom, J., 2013. Scenarios of global municipal water-use demand projections over the 21st century. Hydrological Sciences Journal 58(3): 519-538. DOI: http://dx.doi.org/10.1080/02626667.2013.772301

[5]   Hargreaves, G. L., Hargreaves, G. H. and Riley, J. P., 1985. Irrigation water requirements for Senegal River basin. Journal of Irrigation and Drainage Engineering 111(3): 265-275. DOI: http://dx.doi.org/10.1061/(ASCE)0733-9437(1985)111:3(265)

[6]   Zhou, Y., Hejazi, M., Smith, S., Edmonds, J., Li, H., Clarke, L., Calvin, K. and Thomson, A., 2015. A comprehensive view of global potential for hydro-generated electricity. Energy & Environmental Science 8: 2622-2633. DOI: http://dx.doi.org/10.1039/C5EE00888C

[7]   Kim, S.H., Hejazi, M., Liu, L., Calvin, K., Clarke, L., Edmonds, J., Kyle, P., Patel, P., Wise, M. and Davies, E., 2016. Balancing global water availability and use at basin scale in an integrated assessment model. Climatic Change 136(2): 217-231. DOI: http://dx.doi.org/10.1007/s10584-016-1604-6

[8]   Leng, G., Tang, Q. and Rayburg, S., 2015. Climate change impacts on meteorological, agricultural and hydrological droughts in China. Global and Planetary Change 126:23-34. DOI:        http://doi.org/10.1016/j.gloplacha.2015.01.003

[9]   Fekete, B.M., Vörösmarty, C.J. and Grabs, W., 1999. Global, composite runoff fields based on observed river discharge and simulated water balances. Volume 22 of GRDC-Report. Global Runoff Data Centre, Federal Institute of Hydrology, Koblenz, Germany.

[10] Fekete, B.M., Vörösmarty, C.J. and Grabs, W., 2002. High-resolution fields of global runoff combining observed river discharge and simulated water balances. Global Biogeochemical Cycles 16(3):15-1-15-10. DOI: http://dx.doi.org/10.1029/1999GB001254

[11] Fekete, B.M. and Vörösmarty, C.J., 2011. ISLSCP II UNH/GRDC Composite Monthly Runoff. In Hall, Forrest G., G. Collatz, B. Meeson, S. Los, E. Brown de Colstoun, and D. Landis (eds.). ISLSCP Initiative II Collection. Data set. Available on-line [http://daac.ornl.gov/] from Oak Ridge National Laboratory Distributed Active Archive Center, Oak Ridge, Tennessee, USA. DOI: http://dx.doi.org/10.3334/ORNLDAAC/994

[12] Dai, A., Qian, T., Trenberth, K.E. and Milliman, J.D., 2009. Changes in continental freshwater discharge from 1948 to 2004. Journal of Climate 22(10): 2773-2792. DOI: http://dx.doi.org/10.1175/2008JCLI2592.1
