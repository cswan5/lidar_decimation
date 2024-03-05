# README
## Effects of point density on interpretability of lidar-derived forest structure metrics in two temperate forests

This R code was developed for lidar processing and analysis related to the manuscript entitled "Effects of point density on interpretability of lidar-derived forest structure metrics in two temperate forests" by A.C. Swanson et al. 

### Included files and descriptions
#### `00_lidar_decimation_functions.R`
This file contains the base functions used for decimating the lidar and quantifying forest structure metrics. It is called in `02_sensitivity analysis.R` and `03_lidar_decimation_processing.R`. Some of the functions used within this code were pulled from the `leafR` and `lidR` packages (see citations). 

#### `01_site_select.R`
This file contains code to select plots from the lidar data based on point density.

#### `02_sensitivity_analysis.R`
This file is used to run the sensitivity analysis to determine the ideal number of random decimation iterations to minimize standard error of mean.

#### `03_lidar_decimation_processing.R`
This is code for parallel processing the random lidar decimations and quantifying forest structure metrics from the decimated lidar.

#### `04_statistical_analysis.R`
This file contains the code to run all of the statistical analyses that were performed for the manuscript.

#### `05_figures.R`
This file contains the code to generate the figures contained in the manuscript.

### Citations
de Almeida, D.R.A., S.C. Stark, C.A. Silva, C. Hamamura, and R. Valbuena. 2021. leafR: Calculates the Leaf Area Index (LAD) and Other Related Functions. R package version 0.3.5. https://CRAN.R-project.org/package=leafR

Roussel, J.R. and D. Auty. 2023. Airborne LiDAR Data Manipulation and Visualization for Forestry Applications. R package version 4.0.4. https://cran.r-project.org/package=lidR
