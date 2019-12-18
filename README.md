# Precipitation-Data-GHCN-station-observations

Analysis of GHCN (www.ncdc.noaa.gov/ghcn-daily-description/) station observation data of daily precipitation to calculate trends in Mean and Extreme precipitation over the period of 1951-2010, and their response to global warming.

Download GHCN-daily data from ( ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/ ) and place them in the “ghcn_all” directory.

Code for reading data: GHCH_Data_Initialization.m
This codes reads GHCN-daily station observation data, separates the daily precipitation fields (the data files include other meteorological data as well), excludes stations that have less than 20 years worth of data, and and saves the station data in .mat files to be analyzed by a different code. The matlab data are saved in “Matlab Data 20years” directory.


Code for analyzing data: GHCN_PrecipitationTrend_Mean_Max_1950_2010.m
This code reads .mat files created by the previous code and calculates average mean and extreme precipitation for each station and plots them on a Geo-referenced global map. It also calculates trends in mean and extreme precipitation over a 60 year period in absolute terms, relative terms, and their response to 1 degree of global warming. NASA-GISS global mean near surface temperature are used.

* NOTE: most calculations here are coded extensively and simple matlab functions are not used. For example, the code itself does linear regression calculation on each station and handle NaN values and saves them into an array. Extracting Year/Month/Day values from Time variables are done manually. Percentiles are calculated manually. All of these calculations probably can be done easier and faster using other available functions.


Results are partially published in: https://onlinelibrary.wiley.com/doi/full/10.1111/1752-1688.12472


