#Use Geopandas environment
#Makes a single netcdf with masks for different countries given on a country dimension with ISO3 tags
#python partition_world_by_country.py -grid '1x1' -grid2Agg '4.0x5.0' -o '/hpc/group/shindell/ap851/masks/partition.nc'
#To partition by population, use (parallelization example)
#python partition_world_by_country.py -grid '4.0x5.0' -aggByPop 'True' -c 0 -n 10000 -country 'USA,CAN,JPN,CHN' -o '/hpc/group/shindell/ap851/masks/partition_by_pop.nc'
#To merge the 100 group from the sbatch file use,
#python partition_world_by_country.py -grid '4.0x5.0' -aggByPop 'True' -n 100 -o '/hpc/group/shindell/ap851/masks/partition_by_pop.nc'

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point
import argparse
import xarray as xr
from datetime import datetime
from simple_utils import *

testing = False
debugMode = True

if testing:
	file_out = '/hpc/group/shindell/ap851/mask_test.nc'
	gridlabel = "4.0x5.0"
	countrystring="USA,CAN,JPN,CHN"
	latlon_name = "lat lon"
	grid2Agg = None
	aggByPop = False
	path2pop=None
	path2poplandmask=None
	parallelization_code = None
	parallelization_n_tasks = None
else:
	parser = argparse.ArgumentParser(description='Make mask for a given country.')
	parser.add_argument('-grid', '--grid_label', type=str, help='Label of GEOS-Chem grid (e.g. AS_MERRA2, 2.0x2.5), or a netcdf file with latitude/longitude as dimensions (will make mask with same dimensions). Curvilinear is fine.')
	parser.add_argument('-grid2Agg', '--grid_label_for_agg', default=None, type=str, help='If supplied, will aggregate the grid above to this resolution by area; mask will then represent percent of grid cell in boundary.')
	parser.add_argument('-aggByPop', '--grid_percent_by_population', default=False, type=str2bool, help='If True, will aggregate to the grid_label resolution, with boundary cells weighted between 0 and 1 by population; mask will then represent percent of grid cell in boundary weighted by population.')
	parser.add_argument('-path2pop', '--path_to_population_gpw4_file', default="/hpc/group/shindell/ap851/masks/gpwv4/gpw_v4_une_atotpopbt_cntm_2pt5_min.nc", type=str, help='If aggByPop is True, use this file of total population counts from GPWv4 at 2.5 arcmin to calculate population weights.')
	parser.add_argument('-path2poplandmask', '--path_to_population_gpw4_landmask_file', default="/hpc/group/shindell/ap851/masks/gpwv4/gpw_landmask.nc", type=str, help='If aggByPop is True, use this file of to apply a landmaks for the population counts to avoid coastline errors.')
	parser.add_argument('-latlon_name', '--latlon_name', type=str, default="lat lon", help='If grid supplied is a netcdf file, name of latitude/longitude dimensions respectively, space delimited. Example: "lat lon".')
	parser.add_argument('-c', '--parallelization_code', type=int, default=None, help='If parallelizing and in step one (generating high-res mask), code for this worker portion of the task. If left as none but parallelizing, assumes we are in step 2')
	parser.add_argument('-n', '--parallelization_n_tasks', type=int, default=None, help='If parallelizing, number of workers. If a value is supplied here, we assume that parallization is happening.')
	parser.add_argument('-country', '--country', type=str, default='ABW,AFG,AGO,AIA,ALB,AND,ARE,ARG,ARM,ASM,ATF,ATG,AUS,AUT,AZE,BDI,BEL,BEN,BFA,BGD,BGR,BHR,BHS,BIH,BLM,BLR,BLZ,BMU,BOL,BRA,BRB,BRN,BTN,BWA,CAF,CAN,CHE,CHL,CHN,CIV,CMR,COD,COG,COK,COL,COM,CPV,CRI,CUB,CUW,CYM,CYP,CZE,DEU,DJI,DMA,DNK,DOM,DZA,ECU,EGY,ERI,ESP,EST,ETH,FIN,FJI,FLK,FRO,FSM,GAB,GBR,GEO,GGY,GHA,GIB,GIN,GMB,GNB,GNQ,GRC,GRD,GRL,GTM,GUM,GUY,HKG,HMD,HND,HRV,HTI,HUN,IDN,IMN,IND,IOT,IRL,IRN,IRQ,ISL,ISR,ITA,JAM,JEY,JOR,JPN,KAZ,KEN,KGZ,KHM,KIR,KNA,KOR,KWT,LAO,LBN,LBR,LBY,LCA,LIE,LKA,LSO,LTU,LUX,LVA,MAC,MAF,MAR,MCO,MDA,MDG,MDV,MEX,MHL,MKD,MLI,MLT,MMR,MNE,MNG,MNP,MOZ,MRT,MSR,MUS,MWI,MYS,NAM,NCL,NER,NFK,NGA,NIC,NIU,NLD,NPL,NRU,NZL,OMN,PAK,PAN,PCN,PER,PHL,PLW,PNG,POL,PRI,PRK,PRT,PRY,PSE,PYF,QAT,ROU,RUS,RWA,SAU,SDN,SEN,SGP,SGS,SHN,SLB,SLE,SLV,SMR,SOM,SPM,SRB,SSD,STP,SUR,SVK,SVN,SWE,SWZ,SXM,SYC,SYR,TCA,TCD,TGO,THA,TJK,TKM,TLS,TON,TTO,TUN,TUR,TUV,TZA,UGA,UKR,UMI,URY,USA,UZB,VAT,VCT,VEN,VGB,VIR,VNM,VUT,WLF,WSM,YEM,ZAF,ZMB,ZWE', help='Comma separated three letter ISO codes for countries to make masks of. Each country will be placed in the country dimension')
	parser.add_argument('-o', '--file_out', type=str, help='Where to save mask as netcdf.')
	args = parser.parse_args()
	file_out = args.file_out
	gridlabel = args.grid_label
	grid2Agg = args.grid_label_for_agg
	aggByPop=args.grid_percent_by_population
	path2pop=args.path_to_population_gpw4_file
	path2poplandmask=args.path_to_population_gpw4_landmask_file
	countrystring=args.country
	latlon_name = args.latlon_name
	parallelization_code = args.parallelization_code
	parallelization_n_tasks = args.parallelization_n_tasks


if (grid2Agg is not None) and aggByPop:
	raise ValueError('grid2Agg is supplied and aggByPop is true, but have to choose between areal or population weighting for divided grid cells')

if parallelization_n_tasks is not None:
	do_par = True 
	if parallelization_code is not None:
		step1 = True
	else:
		step1 = False 
else:
	do_par=False
	step1=False

countries = countrystring.split(',')
countrieswcode = {c:i+1 for i,c in enumerate(countries)}

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

mask = []
# Iterate through all horizontal points in cube, and
# check for containment within the specified geometry.
if do_par and step1:
	# Reshape to 1D for easier iteration.
	lon2full = lon2d.reshape(-1)
	lat2full = lat2d.reshape(-1)
	# Split to portion of the work assigend.
	lon2 = np.array_split(lon2full, parallelization_n_tasks)[parallelization_code]
	lat2 = np.array_split(lat2full, parallelization_n_tasks)[parallelization_code]
else:
	# Reshape to 1D for easier iteration.
	lon2 = lon2d.reshape(-1)
	lat2 = lat2d.reshape(-1)


if (not do_par) or (do_par and step1):
	for latval, lonval in zip(lat2, lon2):
		this_point = Point(lonval, latval)
		to_add = 0
		for country in countries:
			cs = geoms[country].contains(this_point)
			if (len(cs) > 0) and cs.values[0]:
				to_add = countrieswcode[country]
				break
		mask.append(to_add)
else:
	#We've already done the parallelization, and now are just completing everything
	mask_parts = []
	for i in range(parallelization_n_tasks):
		mask_parts.append(np.load(f'{file_out}_parworker_{i}'))
	mask = np.concatenate(mask_parts)

if do_par and step1:
	#save out intermediate parallelization data
	mask = np.save(f'{file_out}_parworker_{parallelization_code}',np.array(mask))
else:
	#Save out proper mask
	mask = np.array(mask).reshape(lon2d.shape)
	to_return = np.zeros((len(countries),mask.shape[0],mask.shape[1]))
	for country in countries:
		to_return[countrieswcode[country]-1,:,:] = mask==countrieswcode[country] #Subtract 1 to 0 index

	if (grid2Agg is not None) or aggByPop:
		if debugMode:
			np.save(f'{file_out}_fine_onegrid',mask)
			np.save(f'{file_out}_fine_multidim',to_return)
		coarse_lon2d, coarse_lat2d = np.meshgrid(coarse_lon,coarse_lat)
		mask2return = np.zeros((len(countries),len(coarse_lat),len(coarse_lon)))
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
					temp = to_return[countrieswcode[country]-1,:,:]
					if grid2Agg is not None:
						percent = np.sum(temp[valid])/np.sum(valid)
					elif aggByPop:
						denom = np.sum(popdata.values[valid])
						if denom==0:
							percent = 0 #No population in gridcell
						else:
							percent = np.sum(popdata.values[valid]*temp[valid])/denom
					mask2return[countrieswcode[country]-1,j,i] = percent
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


