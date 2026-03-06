#Use Geopandas environment
#Makes a mask for a given country (3 letter ISO code)
#ex (Korea): python make_country_landmask.py -grid AS_MERRA2 -country KOR -o /n/holylfs05/LABS/jacob_lab/Users/drewpendergrass/korea_comp_data/GC/korea_mask.npy
#ex (NO2 ratio masks): python make_country_landmask.py -grid '/n/holylfs05/LABS/jacob_lab/Users/drewpendergrass/korea_comp_data/sat/monthly_hcho_no2_ratio.nc' -country KOR -o "/n/holylfs05/LABS/jacob_lab/Users/drewpendergrass/korea_comp_data/sat/no2_satratio_mask.npy"
#ex population masks: python make_country_landmask.py -grid '0.5x0.5' -aggByPop True -continent "South America" -country "FRA" -lon "m180,m30" -custom "SUBREGION:Caribbean,Central America" -o /hpc/group/shindell/ap851/masks/south_and_central_america_0p5x0p5_pop_weighted_mask.nc4
#ex population masks, south korea test: python make_country_landmask.py -grid '0.5x0.5' -aggByPop True -country "KOR" -lon "123,132" -lat "32,40" -o /hpc/group/shindell/ap851/masks/skorea_0p5x0p5_pop_weighted_mask.npy

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point
import argparse
import xarray as xr
from datetime import datetime
from simple_utils import *

testing = False
debugMode=False

if testing:
	file_out = '/n/holylfs05/LABS/jacob_lab/Users/drewpendergrass/korea_comp_data/sat/testmask.npy'
	gridlabel = "2.0x2.5"
	countrystring=""
	custom_grouping=""
	continent = ""
	latbounds = "m90 0".split()
	lonbounds = "".split()
	latlon_name = "lat lon"
	grid2Agg = None
	aggByPop = False
	path2pop=None
	path2poplandmask=None
else:
	parser = argparse.ArgumentParser(description='Make mask for a given country.')
	parser.add_argument('-grid', '--grid_label', type=str, help='Label of GEOS-Chem grid (e.g. AS_MERRA2, 2.0x2.5), or a netcdf file with latitude/longitude as dimensions (will make mask with same dimensions). Curvilinear is fine.')
	parser.add_argument('-grid2Agg', '--grid_label_for_agg', default=None, type=str, help='If supplied, will aggregate the grid above to this resolution by area; mask will then represent percent of grid cell in boundary.')
	parser.add_argument('-aggByPop', '--grid_percent_by_population', default=False, type=str2bool, help='If True, will aggregate to the grid_label resolution, with boundary cells weighted between 0 and 1 by population; mask will then represent percent of grid cell in boundary weighted by population.')
	parser.add_argument('-path2pop', '--path_to_population_gpw4_file', default="/hpc/group/shindell/ap851/masks/gpwv4/gpw_v4_une_atotpopbt_cntm_2pt5_min.nc", type=str, help='If aggByPop is True, use this file of total population counts from GPWv4 at 2.5 arcmin to calculate population weights.')
	parser.add_argument('-path2poplandmask', '--path_to_population_gpw4_landmask_file', default="/hpc/group/shindell/ap851/masks/gpwv4/gpw_landmask.nc", type=str, help='If aggByPop is True, use this file of to apply a landmaks for the population counts to avoid coastline errors.')
	parser.add_argument('-latlon_name', '--latlon_name', type=str, default="lat lon", help='If grid supplied is a netcdf file, name of latitude/longitude dimensions respectively, space delimited. Example: "lat lon".')
	parser.add_argument('-country', '--country', type=str, default='', help='Three letter ISO code for country to make mask of. If multiple, comma separated. Use tilde to exclude a country. To only get country north of a latitude value, use colon and then that value (e.g. USA:50,CAN,~MEX)')
	parser.add_argument('-continent', '--continent', type=str, default='', help='Continent to make mask of. Countries can be appended if present')
	parser.add_argument('-custom', '--custom_grouping', type=str, default='', help='Custom grouping -- give name of column (e.g. REGION_WB) colon group (e.g. Sub-Saharan Africa) with groups comma separated, e.g. REGION_WB:Sub-Saharan Africa,North Africa & Middle East')
	parser.add_argument('-lat', '--lat_bounds', type=str, default='', help='Latitude bounds of mask. If provided, comma delimited. Use "m" in place of minus sign (optional)')
	parser.add_argument('-lon', '--lon_bounds', type=str, default='', help='Longitude bounds of mask. If provided, comma delimited. Use "m" in place of minus sign (optional)')
	parser.add_argument('-o', '--file_out', type=str, help='Where to save mask. If ends in csv, saves as csv; if npy, saves as NumPy object; if .nc or .nc4 it saves as a netcdf.')
	args = parser.parse_args()
	file_out = args.file_out
	gridlabel = args.grid_label
	grid2Agg = args.grid_label_for_agg
	aggByPop=args.grid_percent_by_population
	path2pop=args.path_to_population_gpw4_file
	path2poplandmask=args.path_to_population_gpw4_landmask_file
	countrystring=args.country
	custom_grouping=args.custom_grouping
	continent = args.continent
	latbounds = args.lat_bounds.split(',')
	lonbounds = args.lon_bounds.split(',')
	latlon_name = args.latlon_name

if (grid2Agg is not None) and aggByPop:
	raise ValueError('grid2Agg is supplied and aggByPop is true, but have to choose between areal or population weighting for divided grid cells')

if len(custom_grouping)>0:
	custom_grouping = custom_grouping.split(':')
	custom_grouping[1] = custom_grouping[1].split(',')

if len(countrystring)>0:
	countries_raw = countrystring.split(',')
else:
	countries_raw = []

excl_countries = []
countries_north = {}
countries = []
for c in countries_raw:
	if (':' in c) and ('~' in c):
		raise ValueError(f'{c} is not accepted string')
	if c.startswith('~'):
		excl_countries.append(c.split('~')[1])
	elif ':' in c:
		countries_north[c.split(':')[0]] = float(c.split(':')[1])
	else:
		countries.append(c)

world=gpd.read_file('WB_countries_Admin0_10m/WB_countries_Admin0_10m.shp')
geoms = []
excl_geoms = []
north_geoms = {}
excl_shapemask=False
north_shapemask = False
shapemask=True

if len(continent)>0:
	ids = list(world[world['CONTINENT'] == continent].OBJECTID)
	geoms.extend([world[world['OBJECTID'] == i].geometry for i in ids])

if len(custom_grouping)>0:
	ids = list(world[np.isin(world[custom_grouping[0]], custom_grouping[1])].OBJECTID)
	geoms.extend([world[world['OBJECTID'] == i].geometry for i in ids])

if len(countries)>0:
	geoms.extend([world[world['ISO_A3'] == country].geometry for country in countries])

if len(excl_countries)>0:
	excl_geoms.extend([world[world['ISO_A3'] == country].geometry for country in excl_countries])
	excl_shapemask=True

if len(countries_north)>0:
	north_shapemask=True
	for c in countries_north:
		north_geoms[c] = {}
		north_geoms[c]['geom'] = world[world['ISO_A3'] == c].geometry
		north_geoms[c]['lat'] = countries_north[c]

if (len(geoms)==0) and (len(north_geoms)==0):
	shapemask=False

if len(latbounds)==2:
	latbounds = [l.replace('m','-') for l in latbounds]
	latbounds = [float(l) for l in latbounds]
else:
	latbounds = None

if len(lonbounds)==2:
	lonbounds = [l.replace('m','-') for l in lonbounds]
	lonbounds = [float(l) for l in lonbounds]
else:
	lonbounds = None

if not aggByPop:
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
elif aggByPop:
	popdata = xr.open_dataset(path2pop)
	poplandmask = xr.open_dataset(path2poplandmask)
	poplandmask = poplandmask.rename({'lat':'latitude','lon':'longitude'})
	popdata = popdata.fillna(0)
	popdata = popdata.sortby("latitude", ascending=True) #fix lat dimension
	#Select 2020 data
	popdata = popdata['UN-Adjusted Population Count, v4.10 (2000, 2005, 2010, 2015, 2020): 2.5 arc-minutes'].sel(raster=5)
	#Apply pop landmask to avoid coast problems
	popdata = popdata*poplandmask['mask']
	coarse_lon_edge,coarse_lat_edge = getLonLatEdgeFromLabel(gridlabel)
	coarse_lon,coarse_lat = getLonLatFromLabel(gridlabel)
	lon = popdata.longitude.values
	lat = popdata.latitude.values
	lon2d, lat2d = np.meshgrid(lon,lat)



# Reshape to 1D for easier iteration.
lon2 = lon2d.reshape(-1)
lat2 = lat2d.reshape(-1)
mask = []
# Iterate through all horizontal points in cube, and
# check for containment within the specified geometry.
for latval, lonval in zip(lat2, lon2):
	this_point = Point(lonval, latval)
	validLatLon = True
	validCountry = True
	maskadd = 0
	if (latbounds is not None) and validLatLon:
		if (latval<latbounds[0]) or (latval>latbounds[1]):
			validLatLon=False
	if (lonbounds is not None) and validLatLon:
		if (lonval<lonbounds[0]) or (lonval>lonbounds[1]):
			validLatLon=False
	if validLatLon and excl_shapemask:
		for eg in excl_geoms:
			cs = eg.contains(this_point)
			if (len(cs) > 0) and cs.values[0]:
				validCountry = False
				break
	if validLatLon and shapemask and validCountry:
		for geom in geoms:
			cs = geom.contains(this_point)
			if (len(cs) > 0) and cs.values[0]:
				maskadd = 1
				break
		if (maskadd == 0) and north_shapemask:
			for c in north_geoms:
				if (latval>north_geoms[c]['lat']) and north_geoms[c]['geom'].contains(this_point).values[0]:
					maskadd = 1
					break
	elif validLatLon and validCountry and ~shapemask:
		maskadd = 1
	mask.append(maskadd)

mask = np.array(mask).reshape(lon2d.shape)

if (grid2Agg is not None) or aggByPop:
	if debugMode:
		np.save(f'{file_out}_fine',mask)
	coarse_lon2d, coarse_lat2d = np.meshgrid(coarse_lon,coarse_lat)
	mask2return = np.zeros(coarse_lon2d.shape)
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
			if grid2Agg is not None:
				percent = np.sum(mask[valid])/np.sum(valid)
			elif aggByPop:
				denom = np.sum(popdata.values[valid])
				if denom==0:
					percent = 0 #No population in gridcell
				else:
					percent = np.sum(popdata.values[valid]*mask[valid])/np.sum(popdata.values[valid])
			mask2return[j,i] = percent
	mask = mask2return

if file_out.split('.')[-1].lower() == "csv":
	np.savetxt(file_out,mask,delimiter=',')
elif file_out.split('.')[-1].lower() == "npy":
	np.save(file_out,mask)
elif file_out.split('.')[-1].lower() in ["nc","nc4"]:
	if (grid2Agg is not None) or aggByPop:
		lon2save = coarse_lon
		lat2save = coarse_lat
	else:
		lon2save = lon
		lat2save = lat
	ds = xr.Dataset(
		{"mask": (("lat","lon"), mask,{"long_name": "Mask", "units":"1"})},
		coords={
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


