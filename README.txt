This package includes files and MATLAB scripts to calculate morphologic metrics and stream power law erosion quantities for the submitted paper citeable as:

Higa, Justin T., Rowland, Scott K. (in review) Volcanic and climatic controls on million-year valley incision, West Maui Volcano, Hawaiʻi.

Higa_and_Rowland_ValleyTopo_Zenodo.m and Higa_and_Rowland_DivideMigration_Zenodo.m require the TopoToolbox 2 and DivideTools packages.

In addition, bin1.m copyrighted by Taylor Perron was used for simple plotting, provided here.

AcrossDivide.m from DivideTools was modified as AcrossDivideJH.m and is provided here, as well.

Several rasters are provided
The example digital elevation is from the USGS EDNA 10-m resolution digital elevation model of West Maui Volcano as usgs_dem_10m_Wmaui_UTM.asc.
The rainfall map and standard deviation is from the Hawaii Rainfall Atlas (Giambelluca et al., 2013) as rfgrid_mm_maui_ann.asc and uncertaintygrid_mm_state_ann_UTM.asc, respectively.
The distance from Bouguer gravity high is from Strange et al. (1965) as strange_grav_cald_dist.asc.
The reconstructed surface is Wmaui_thinplatespline.asc. Points used to interpolate this surface is Wmaui_thinplatespline_points.csv.

Several shapefiles are provided
The West Maui Volcano caldera outline is adapted from Sherrod et al. (2021) as wailuku_caldera_full.shp.
West Maui valleys analyzed are WMV_Valleys.shp, valleys named in the paper are WMV_Valleys_main.shp.
