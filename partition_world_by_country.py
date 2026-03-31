#Use Geopandas environment
#Makes a single netcdf with masks for different countries given on a country dimension with ISO3 tags

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point
import argparse
import xarray as xr
from datetime import datetime
from simple_utils import *

testing = True
debugMode = False

if testing:
	file_out = '/hpc/group/shindell/ap851/mask_test.nc'
	gridlabel = "4.0x5.0"
	countrystring=""
	latlon_name = "lat lon"
	grid2Agg = "USA,CAN,KOR"
else:
	parser = argparse.ArgumentParser(description='Make mask for a given country.')
	parser.add_argument('-grid', '--grid_label', type=str, help='Label of GEOS-Chem grid (e.g. AS_MERRA2, 2.0x2.5), or a netcdf file with latitude/longitude as dimensions (will make mask with same dimensions). Curvilinear is fine.')
	parser.add_argument('-grid2Agg', '--grid_label_for_agg', default=None, type=str, help='If supplied, will aggregate the grid above to this resolution by area; mask will then represent percent of grid cell in boundary.')
	parser.add_argument('-latlon_name', '--latlon_name', type=str, default="lat lon", help='If grid supplied is a netcdf file, name of latitude/longitude dimensions respectively, space delimited. Example: "lat lon".')
	parser.add_argument('-country', '--country', type=str, default='', help='Comma separated three letter ISO codes for countries to make masks of. Each country will be placed in the country dimension')
	parser.add_argument('-o', '--file_out', type=str, help='Where to save mask as netcdf.')
	args = parser.parse_args()
	file_out = args.file_out
	gridlabel = args.grid_label
	grid2Agg = args.grid_label_for_agg
	countrystring=args.country
	latlon_name = args.latlon_name

countries = countrystring.split(',')
countrieswcode = {t:i for i,t in enumerate(test)}

world=gpd.read_file('WB_countries_Admin0_10m/WB_countries_Admin0_10m.shp')
geoms = {country:world[world['ISO_A3'] == country].geometry for country in countries}

try: #if a grid label, this will work.
	lon,lat = getLonLatFromLabel(gridlabel)
	lon2d, lat2d = np.meshgrid(lon,lat)
except:
	ds = xr.load_dataset(gridlabel)
	if grid2Agg is not None:
		raise ValueError("Can only aggregate non-custom lat/lon grids.")
	latlab,lonlab = latlon_name.split()
	lon = ds[lonlab].values
	lat = ds[latlab].values
	if len(lon.shape)==2:
		lon2d=lon
		lat2d=lat
	else:
		lon2d, lat2d = np.meshgrid(lon,lat)

if grid2Agg is not None:
	coarse_lon_edge,coarse_lat_edge = getLonLatEdgeFromLabel(grid2Agg)
	coarse_lon,coarse_lat = getLonLatFromLabel(grid2Agg)

# Reshape to 1D for easier iteration.
lon2 = lon2d.reshape(-1)
lat2 = lat2d.reshape(-1)
mask = []
# Iterate through all horizontal points in cube, and
# check for containment within the specified geometry.
for latval, lonval in zip(lat2, lon2):
	this_point = Point(lonval, latval)
	for country in countries:
		cs = geoms[country].contains(this_point)
		if (len(cs) > 0) and cs.values[0]:
			mask.append(countrieswcode[country])
			break

mask = np.array(mask).reshape(lon2d.shape)
to_return = np.zeros((len(countries),mask.shape[0],mask.shape[1]))
for country in countries:
	to_return[countrieswcode[country],:,:] = mask==countrieswcode[country]


if grid2Agg is not None:
	if debugMode:
		np.save(f'{file_out}_fine_onegrid',mask)
		np.save(f'{file_out}_fine_multidim',to_return)
	coarse_lon2d, coarse_lat2d = np.meshgrid(coarse_lon,coarse_lat)
	mask2return = np.zeros(to_return.shape)
	for i in range(len(coarse_lon)):
		lonstart,lonstop = coarse_lon_edge[i:i+2]
		if lonstart<lonstop: # normal case
			validlon = (lon2d>=lonstart)&(lon2d<lonstop)
		else:
			validlon = (lon2d>=lonstart)|(lon2d<lonstop) #date line
		for j in range(len(coarse_lat)):
			latstart,latstop = coarse_lat_edge[j:j+2]
			validlat = (lat2d>=latstart)&(lat2d<latstop)
			valid = validlon&validlat
			for country in countries:
				temp = to_return[countrieswcode[country],:,:]
				percent = np.sum(temp[valid])/np.sum(valid)
				mask2return[countrieswcode[country],j,i] = percent
	to_return = mask2return

if grid2Agg is not None:
	lon2save = coarse_lon
	lat2save = coarse_lat
else:
	lon2save = lon
	lat2save = lat
ds = xr.Dataset(
	{"mask": (("country","lat","lon"), to_return,{"long_name": "Mask", "units":"1"})},
	coords={
		"country": (["country"], countries,{"long_name": "Country ISO3 code", "units":"1"}),
		"lat": (["lat"], lat2save,{"long_name": "Latitude", "units":"degrees_north"}),
		"lon": (["lon"], lon2save,{"long_name": "Longitude", "units":"degrees_east"})
	},
	attrs={
		"Title":'Mask',
		"Conventions":"COARDS",
		"Format":"NetCDF-4",
		"Model":"GENERIC",
		"NLayers":"1",
		"History":f"Originally calculated by the country mask utility on {str(datetime.today())}",
		"Start_Date":"2020-01-01",
		"Start_Time":"0",
		"End_Date":f"2020-12-01",
		"End_Time":"0"
	}
)
ds.to_netcdf(file_out)


