This package includes files and MATLAB scripts to calculate morphologic metrics and stream power law erosion quantities for the submitted paper citable as:

Higa, Justin T., Rowland, Scott K. (in review) Volcanic and climatic controls on million-year valley incision, West Maui Volcano, Hawaiʻi.

Higa_and_Rowland_ValleyTopo_Zenodo.m and Higa_and_Rowland_DivideMigration_Zenodo.m require the TopoToolbox 2 (Scherler & Schwanghart, 2020; Schwanghart & Scherler, 2014) and DivideTools (Forte & Whipple, 2018) packages. For Topotoolbox 2, the refmat2XY.m function must be taken out of the private folder in the @GRIDobj folder and placed in a findable location (e.g., directly in the @GRIDobj folder). RGB triple of color name, version 2 (https://www.mathworks.com/matlabcentral/fileexchange/24497-rgb-triple-of-color-name-version-2; last accessed 15 May 2025) is also required for plotting.

Several MATLAB scripts are provided:
1) meanbin.m was used for simple plotting of binned data.
2) AcrossDivideJH.m is modified from AcrossDivideJH.m in DivideTools to include the 'ultralax' option for the 'divide_buffer' optional input.

Several rasters are provided:
1) The distance from Bouguer gravity high is from Strange et al. (1965) as strange_grav_cald_dist.asc.
2) The reconstructed surface is Wmaui_thinplatespline.asc, zipped as Wmaui_thinplatespline.zip. Points used to interpolate this surface is Wmaui_thinplatespline_points.csv.

Several shapefiles are provided:
1) The West Maui Volcano caldera outline is adapted from Sherrod et al. (2021a,b) as wailuku_caldera_full.shp.
2) West Maui valleys analyzed are WMV_Valleys.shp, valleys named in the paper are WMV_Valleys_main.shp.

Several tables are provided:
1) maui_valley_aspect.csv holds the valley ID number, valley widths [m], and corresponding valley lengths [m]. There is one row for each width measured, where the ID and valley length is constant for each valley.
2) Wmaui_thinplatespline_points.csv holds the Easting, Northing, elevation, and point type for points used in the reconstructed surface thin-plate spline.

Several rasters must be downloaded from other sources:
1) The digital elevation model is the USGS EDNA 10-m resolution digital elevation model of Hawaiʻi (downloadable from https://earthexplorer.usgs.gov/; last accessed 15 May 2025). See http://edna.usgs.gov for more information.
2) The rainfall map and standard deviation is from the Hawaiʻi Rainfall Atlas (Giambelluca et al., 2013; downloadable from https://www.hawaii.edu/climate-data-portal/rainfall-atlas/; last accessed 15 May 2025). These are the "Rainfall Grids" and "Uncertainty" maps in ASCII Grid Format, millimeters.
*The MATLAB code has a line to clip all rasters, including these that must be downloaded, to the same extent as that in the paper. The exact extent of self-downloaded rasters does not have to match that of the rasters provided in order to work. However, the extent of the self-downloaded rasters, particularly the EDNA elevation, must not be smaller than the extent of the provided ones. Also, if self-downloaded rasters are too large, this may affect run time.

Citations:

Forte, A. M., & Whipple, K. X. (2018). Criteria and tools for determining drainage divide stability. Earth and Planetary Science Letters, 493, 102–117. https://doi.org/10.1016/j.epsl.2018.04.026

Giambelluca, T. W., Chen, Q., Frazier, A. G., Price, J. P., Chen, Y.-L., Chu, P.-S., et al. (2013). Online Rainfall Atlas of Hawai‘i. Bulletin of the American Meteorological Society, 94(3), 313–316. https://doi.org/10.1175/BAMS-D-11-00228.1

Jonasson, K. (2025). RGB triple of color name, version 2 (https://www.mathworks.com/matlabcentral/fileexchange/24497-rgb-triple-of-color-name-version-2), MATLAB Central File Exchange. Retrieved May 15, 2025. 

Scherler, D., & Schwanghart, W. (2020). Drainage divide networks – Part 1: Identification and ordering in digital elevation models. Earth Surface Dynamics, 8(2), 245–259. https://doi.org/10.5194/esurf-8-245-2020

Schwanghart, W., & Scherler, D. (2014). Short Communication: TopoToolbox 2 – MATLAB-based software for topographic analysis and modeling in Earth surface sciences. Earth Surface Dynamics, 2(1), 1–7. https://doi.org/10.5194/esurf-2-1-2014

Sherrod, D. R., Robinson, J., Sinton, J., Watkins, S., & Brunt, K. (2021a). Geologic map database to accompany geologic map of the State of Hawai ‘i: US Geological Survey data release. https://doi.org/10.5066/P9YWXT41

Sherrod, D. R., Sinton, J. M., Watkins, S. E., & Brunt, K. M. (2021b). Geologic Map of the State of Hawaiʻi (No. 3143). U.S. Geological Survey. https://doi.org/10.3133/sim3143

Strange, W. E., Woollard, G. P., & Rose, J. C. (1965). An Analysis of the Gravity Field Over the Hawaiian Islands in Terms of Crustal Structure. Pacific Science, 19(3), 381–389.
