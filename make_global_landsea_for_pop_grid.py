#Use Geopandas environment
#Make a global land/sea mask based on the GPWv4 grid for consistent calculations
#Parallelized for rapid computation


import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point
import argparse
import xarray as xr
from datetime import datetime
from simple_utils import *

parser = argparse.ArgumentParser(description='Make global land/sea mask on GPWv4 grid.')
parser.add_argument('-path2pop', '--path_to_population_gpw4_file', default="/hpc/group/shindell/ap851/masks/gpwv4/gpw_v4_une_atotpopbt_cntm_2pt5_min.nc", type=str, help='If aggByPop is True, use this file of total population counts from GPWv4 at 2.5 arcmin to calculate population weights.')
parser.add_argument('-n', '--number_parallel_elements', type=int, help='Number of parallel elements to split the lat/lon grid into. 12 divides lat/lon into 12 each, for 144 squares.')
parser.add_argument('-i', '--parallel_index', type=int, help='Parallelization index, ranging from 0 to one less than the square of number_parallel_elements')
parser.add_argument('-o', '--file_out', type=str, help='Where to save mask. Automatically saves as npy and adds index details.')
args = parser.parse_args()
npar = args.number_parallel_elements
ind = args.parallel_index
file_out = args.file_out
path2pop=args.path_to_population_gpw4_file

world=gpd.read_file('WB_countries_Admin0_10m/WB_countries_Admin0_10m.shp')
geoms = list(world['geometry'].values)

popdata = xr.open_dataset(path2pop)
popdata = popdata.sortby("latitude", ascending=True) #fix lat dimension
#Select 2020 data
popdata = popdata['UN-Adjusted Population Count, v4.10 (2000, 2005, 2010, 2015, 2020): 2.5 arc-minutes'].sel(raster=5)
globallon = popdata.longitude.values
globallat = popdata.latitude.values

#calculate latbounds two numbers
lat_chunks = np.array_split(globallat, npar)
lon_chunks = np.array_split(globallon, npar)
latind = ind%npar 
lonind = int(np.floor(ind/npar))
lat = lat_chunks[latind]
lon = lon_chunks[lonind]

lon2d, lat2d = np.meshgrid(lon,lat)



# Reshape to 1D for easier iteration.
lon2 = lon2d.reshape(-1)
lat2 = lat2d.reshape(-1)
mask = []
# Iterate through all horizontal points in cube, and
# check for containment within the specified geometry.
for latval, lonval in zip(lat2, lon2):
	this_point = Point(lonval, latval)
	maskadd = 0
	for geom in geoms:
		cs = geom.contains(this_point)
		if (len(cs) > 0) and cs.values[0]:
			maskadd = 1
			break
	mask.append(maskadd)

mask = np.array(mask).reshape(lon2d.shape)

np.save(f'{file_out}_lat_{latind}_lon_{lonind}_npar_{npar}.npy',mask)
