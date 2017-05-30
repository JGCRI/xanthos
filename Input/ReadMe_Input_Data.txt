#Date: 05/05/2017
#author: Xinya Li (xinya.li@pnl.gov)
#Project: Xanthos V1.0


<Climate>

NMonth: Number of Months (12 X NYear)
These files are prepared by the user.

1. Precipitation
Monthly average precipitation, positive values
Unit: mm/month
Dimension: 67420*NMonth
Format: MATLAB file (*.mat) or NETCDF file (*.nc), variable name is needed in configuration file

2. Temperature
Monthly average temperature
Unit: Celsius
Dimension: 67420*NMonth
Format: MATLAB file (*.mat) or NETCDF file (*.nc), variable name is needed in configuration file

3. Daily temperature range
Monthly average of range between maximum and minimum daily temperatures, positive values
Unit: Celsius
Dimension: 67420*NMonth
Format: MATLAB file (*.mat) or NETCDF file (*.nc), variable name is needed in configuration file


<GriddedMaps>
1. <harmonized_inputs>
Uniformed dimension and format in this folder: one-row header csv file, 67420*1, 0 means no assignment
   Basin ID for each grid: 235 Basins

   Country ID for each grid: 249 Countries

   Region ID for each grid: 32 regions
   
   The above three data files are GCAM related.

   Maximum Soil Moisture (MSM) (mm/month)
   References:
   FAO: Digital soil of the world and derived soil properties, Rome, 1998; 
   FAO: Digital soil map of the world and derived soil properties, [Version 3.5] Edn., FAO, Rome, Italy, 2003

e.g. basin.csv, country.csv, region32_grids.csv, soil_moisture.csv

2. Grid Area
Area value for each land grid cell: 67420 x 1, unit is ha, convert to km2 by *0.01
e.g. Grid_Areas_ID.csv

3. Coordinates
Coordinates for flattened grid:  67420 x 5, the columns are ID#, longitude, latitude, ilon, ilat
e.g. coordinates.csv

4. Name of basins
Corresponding with Basin ID Map, 235 Basin Names: 1D String Array
e.g. BasinNames235.txt

5. Name of countries
Corresponding with Country ID Map, 249 Country Names: 2D String Array
e.g. country-names.csv

6. Name of regions
Corresponding with Region ID Map, 32 Region Names: 2D String Array
e.g. Rgn32Names.csv

7. Cell ID of lakes
Water Bodies: assign MSM = 999, 306 x 2, Col 1 is the cell number in 67420
e.g. Lakes_wo_casp.csv

8. Cell ID of additional water bodies
Additional water bodies: assign MSM = 999, 421 x 2, Col 1 is the cell number in 67420
e.g. Addit_water421.csv

9. DRT (Hierarchical dominant river tracing algorithm) data for river routing
DRT data, 280 x 720, -9999 for missing values
e.g. DRT_half_FDR_globe_bystr50.txt

DRT data, 280 x 720, -9999 for missing values
e.g. DRT_half_FDISTANCE_globe.txt

Reference:
Wu, H., J. S. Kimball, H. Li, M. Huang, L. R. Leung, and R. F. Adler (2012), A new global river network database for macroscale hydrologic modeling, Water Resour. Res., 48, W09701, doi:10.1029/2012WR012313.


<Diagnostics>

1. VIC model
Annual gridded runoff of 30 years from 1971-2000
Unit: km^3/year
Dimension: 67420*30
Format: MATLAB file (*.mat), variable name is "q"
Reference: Leng, G., Tang, Q. and Rayburg, S., 2015. Climate change impacts on meteorological, agricultural and hydrological droughts in China. Global and Planetary Change 126:23-34. DOI: http://doi.org/10.1016/j.gloplacha.2015.01.003
e.g. vic_watch_hist_nosoc_co2_qtot_global_annual_1971_2000.mat
Note: Xanthos will use the averaged runoff according to the time period of climate forcing, e.g. climate forcing is in 1996-2005, the VIC data in 1996-2000 will be averaged for comparison.

2. WBM model
Gridded runoff
Unit: km^3/year
Dimension: 67420*2
Format: csv file (*.csv), no header
Pre-processing:
Each csv file has two columns: 1st - index of the grid in 67420 grids; 2nd - runoff. Xanthos will read in both column and assign girds not listed to zeros.
Reference: Fekete, B.M., Vörösmarty, C.J. and Grabs, W., 1999. Global, composite runoff fields based on observed river discharge and simulated water balances. Volume 22 of GRDC-Report. Global Runoff Data Centre, Federal Institute of Hydrology, Koblenz, Germany.
e.g. wbm_qestimates.csv, wbmc_qestimates.csv

3. UNH_GRDC model
Averaged gridded runoff of Year 1986-1995
Unit: km^3/year
Dimension: 67420*1
Format: MATLAB file (*.mat), variable name is "q"
e.g. UNH_GRDC_average_annual_1986_1995.mat
Pre-processing: 
Original Data can be accessed at http://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=994
120+2 columns: longitude, latitude, monthly runoff (mm/mon) from 1986-1995.
Reference: 
Fekete, B.M., Vörösmarty, C.J. and Grabs, W., 2002. High-resolution fields of global runoff combining observed river discharge and simulated water balances. Global Biogeochemical Cycles 16(3):15-1-15-10. DOI: http://dx.doi.org/10.1029/1999GB001254
Fekete, B.M. and Vörösmarty, C.J., 2011. ISLSCP II UNH/GRDC Composite Monthly Runoff. In Hall, Forrest G., G. Collatz, B. Meeson, S. Los, E. Brown de Colstoun, and D. Landis (eds.). ISLSCP Initiative II Collection. Data set. Available on-line [http://daac.ornl.gov/] from Oak Ridge National Laboratory Distributed Active Archive Center, Oak Ridge, Tennessee, USA. DOI: http://dx.doi.org/10.3334/ORNLDAAC/994


<AccessibleWater>
The following two files are corresponding with the basin map and basin names.
1. Baseflow index (BFI)
Dimension: 235*2
Format: csv file (*.csv), one-row header
e.g. bfi_per_basin.csv
Reference: 
Beck HE, van Dijk AIJM, Miralles DG et al. (2013) Global patterns in base flow index and recession based on streamflow observations from 3394 catchments. Water Resources Research 49(12):7843-7863.
Note: Xanthos will only read in the second column 

2. Reservoir capacity at basin level
Dimension: 235*1
Format: csv file (*.csv), no header
e.g. total_reservoir_storage_capacity_BM3.csv
Reference:
Liu L., S. Parkinson, M. Gidden, E. Byers, Y. Satoh, K. Riahi, and B. Forman (2017), Impacts of climate change and competing land-uses on the potential for reservoirs to secure reliable freshwater supply in the world’s large river basins, Environmental Research Letters (in progress).
