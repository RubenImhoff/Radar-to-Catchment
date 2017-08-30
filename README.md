# Radar to Catchment
Radar to Catchment aims to provide the user with an easy way to extract radar and satellite data for the area within a shapefile. The main design reason is for river catchments, but this tool could also be used for other purposes than hydrology.

The Radar to Catchment tool gives the user the following options:
 Extract data from a radar or satellite data set for a given shapefile of an area. The
resulting file is a GeoTiff or an ASCII file, with the following additional possibilities:
  { Choose the pixel resolution of the resulting file.
  { Choose the spatial reference system of the resulting file.
  { Change dBZ into millimeters of precipitation.
  { Let the tool calculate zonal statistics per time step, such as the minimum precipitation value, the maximum precipitation value and the mean value for the area.
 Obtain a .csv-file with, per time step, precipitation data for a certain pixel within the data set.

Radar to Catchment supports the following data sets (more information is available in the user manual):
 NL-Radar NL25 rac mfbs (hdf5 format with adjusted precipitation accumulations from the two KNMI - Royal Netherlands Meteorological Institute - weather radars for the Netherlands. Available at www.climate4impact.eu).
 NL-Radar NL21 rac mfbs (hdf5 format with adjusted precipitation accumulations from the two KNMI - Royal Netherlands Meteorological Institute - weather radars for the Netherlands. Available at www.climate4impact.eu).
 DE-Radar Radolan (ASCII format with adjusted hourly precipitation accumulations from the seventeen DWD (German Weather Service) weather radars for Germany. Available at ftp://ftp-cdc.dwd.de/pub/CDC/grids_germany/hourly/radolan/).
 NL-Radar KNMI Data Centre: Average Monthly Precipitation (NETCDF format with long term average of monthly precipitation for the period 1981-2010, based on an interpolation by kriging of precipitation data from the KNMI - Royal Netherlands Meteorological Institute - manual rain gauge network in the Netherlands. Also applicable for most other KNMI Data Centre data sets. Available at https://data.knmi.nl/datasets/Rd3/4?q=neerslag&dtend=2017-07-13T23:59Z). 
 EU-Radar EUMETNET OPERA (hdf5 format with fifteen minutes accumulations of weather radar data from most European countries. Not available as open data, but a data request is possible at http://eumetnet.eu/activities/observations-programme/current-activities/opera/).
 NASA GPM IMERG HQ Precipitation (hdf5 format with a thirty minutes accumulation of estimated precipitation in mm/h. Available at https://pmm.nasa.gov/data-access/downloads/gpm).

# Code availability
Radar to Catchment is available as python-script in this repository for both Linux and Windows. It is recommended to first read the available user manual before making use of the Radar2Catchment tool. 

# Copyright
Radar to Catchment is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your opinion) any later version. This program is redistributed in the hope that it will be useful, but without any warranty and without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licences/. The producer of this tool is by
no means to be held responsible for the use of and the results of this tool.

Suggestions for changes to increase functionality or to overcome possible errors, are highly
appreciated. Please contact:

Ruben Imhoff: ruben.imhoff@wur.nl or r.o.imhoff@live.nl

# Acknowledgments
The design of Radar to Catchment was part of a Capita Selecta during my study program for a Master's degree in Earth and Environment at the Wageningen University and Research. I would like to take this chance to thank Aart Overeem (KNMI) for supervising me during
these weeks, for his ideas to improve the tool and the possibility to design this tool at the KNMI office in De Bilt, the Netherlands.

Furthermore, I would like to thank Claudia Brauer for being the starting point of this possibility and for her suggested improvements of the tool. And finally, I would like to thank the KNMI (Royal Netherlands Meteorological Institute) and my colleagues in the
office for the provision of data sets and a great place to work.
