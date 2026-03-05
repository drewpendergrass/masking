#Use Geopandas environment
#Make a global land/sea mask based on the GPWv4 grid for consistent calculations
#Example call:
#python merge_landsea_par.py -n 12 -o /hpc/group/shindell/ap851/masks/gpwv4/gpw_landmask.nc -p /hpc/group/shindell/ap851/masks/landmask_scratch/temp


import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point
import argparse
import xarray as xr
from datetime import datetime
from glob import glob
from simple_utils import *

parser = argparse.ArgumentParser(description='Make global land/sea mask on GPWv4 grid.')
parser.add_argument('-path2pop', '--path_to_population_gpw4_file', default="/hpc/group/shindell/ap851/masks/gpwv4/gpw_v4_une_atotpopbt_cntm_2pt5_min.nc", type=str, help='If aggByPop is True, use this file of total population counts from GPWv4 at 2.5 arcmin to calculate population weights.')
parser.add_argument('-n', '--number_parallel_elements', type=int, help='Number of parallel elements to split the lat/lon grid into, same as used in make_global_landsea_for_pop_grid.py call.')
parser.add_argument('-p', '--chunk_pattern', type=str, help='Pattern of chunk file names, same as used for file_out in make_global_landsea_for_pop_grid.py call.')
parser.add_argument('-o', '--file_out', type=str, help='Where to save merged mask, ending in .nc or .nc4.')

args = parser.parse_args()
npar = args.number_parallel_elements
chunk_pattern = args.chunk_pattern
file_out = args.file_out
path2pop=args.path_to_population_gpw4_file


popdata = xr.open_dataset(path2pop)
popdata = popdata.sortby("latitude", ascending=True) #fix lat dimension
#Select 2020 data
popdata = popdata['UN-Adjusted Population Count, v4.10 (2000, 2005, 2010, 2015, 2020): 2.5 arc-minutes'].sel(raster=5)
globallon = popdata.longitude.values
globallat = popdata.latitude.values

# --- helper: boundaries for each chunk in index space ---
def chunk_bounds(chunks):
	sizes = [len(c) for c in chunks]
	starts = np.cumsum([0] + sizes[:-1])
	ends = starts + np.array(sizes)
	return list(zip(starts, ends))  # [(s0,e0), (s1,e1), ...]

lat_chunks = np.array_split(globallat, npar)
lon_chunks = np.array_split(globallon, npar)

lat_bounds = chunk_bounds(lat_chunks)
lon_bounds = chunk_bounds(lon_chunks)

global_arr = np.zeros((len(globallat), len(globallon)))  # or zeros if appropriate

for latind in range(npar):
	lat_s, lat_e = lat_bounds[latind]
	for lonind in range(npar):
		lon_s, lon_e = lon_bounds[lonind]
		tile = np.load(f'{chunk_pattern}_lat_{latind}_lon_{lonind}_npar_{npar}.npy')
		global_arr[lat_s:lat_e, lon_s:lon_e] = tile

ds = xr.Dataset(
	{"mask": (("lat","lon"), global_arr,{"long_name": "Mask", "units":"1"})},
	coords={
		"lat": (["lat"], globallat,{"long_name": "Latitude", "units":"degrees_north"}),
		"lon": (["lon"], globallon,{"long_name": "Longitude", "units":"degrees_east"})
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
