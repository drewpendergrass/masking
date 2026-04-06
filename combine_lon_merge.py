import numpy as np
import argparse
import xarray as xr

parser = argparse.ArgumentParser(description='Make mask for a given country.')
parser.add_argument('-n_par_combine', '--n_par_combine', type=int, help='Number of longitude values, which are used to do the partition_world_by_country parallelization.')
parser.add_argument('-o', '--file_out', type=str, help='Where to save mask as netcdf. Will use this to recognize files to combine.')
args = parser.parse_args()

file_out = args.file_out 
n_par_combine = args.n_par_combine 

for i in range(n_par_combine):
	ds = xr.open_dataset(f'{file_out}_mergepar_{i}.nc')
	if i == 0:
		mask_to_save = ds.mask.values
	else:
		mask_to_save += ds.mask.values
	ds.close()

ds = xr.open_dataset(f'{file_out}_mergepar_0.nc')
ds.mask.values = mask_to_save
ds.to_netcdf(file_out)
ds.close()
