import numpy as np
import argparse

def getLonLatFromLabel(gridlabel):
	#Create traditional GEOS-Chem longitude and latitude centers, as specified by the settings in ens_config.json
	if gridlabel == '4.0x5.0':
		lon = np.arange(-180.0,176.0, 5.0)
		lat = np.concatenate([[-89.0],np.arange(-86.0,87.0, 4.0), [89.0]])
	elif gridlabel == '2.0x2.5':
		lon = np.arange(-180.0,178.0, 2.5)
		lat = np.concatenate([[-89.5],np.arange(-88.0,89.0, 2.0), [89.5]])
	elif gridlabel == '1x1':
		lon = np.arange(-179.5,180.0, 1.0)
		lat = np.arange(-89.5,90, 1.0)
	elif gridlabel == '0.5x0.5':
		lon = np.arange(-179.75,180.0, 0.5)
		lat = np.arange(-89.75,90, 0.5)
	elif (gridlabel == '0.5x0.625') | (gridlabel == 'MERRA2'): #MERRA2 NATIVE GRID
		lon = -180.0 + (0.625*np.arange(0.0,576.0,1.0))
		lat = -90.0 + (0.5*np.arange(0.0,361.0,1.0))
	elif (gridlabel == 'AS_MERRA2'): #ASIA NESTED GRID FOR MERRA2
		lon = np.arange(60.0,150.01, 0.625)
		lat = np.arange(-11.0,55.01, 0.5)
	elif (gridlabel == 'EU_MERRA2'): #EUROPE NESTED GRID FOR MERRA2
		lon = np.arange(-30.0,50.01, 0.625)
		lat = np.arange(30.0,70.01, 0.5)
	elif (gridlabel == 'NA_MERRA2'): #NORTH AMERICA NESTED GRID FOR MERRA2
		lon = np.arange(-140.0,-39.99, 0.625)
		lat = np.arange(10.0,70.01, 0.5)
	elif (gridlabel == '0.25x0.3125') | (gridlabel == 'GEOSFP'):
		lon = -180.0 + (0.3125*np.arange(0.0,1152.0,1.0))
		lat = -90.0 + (0.25*np.arange(0.0,721.0,1.0))
	elif (gridlabel == 'CH_GEOSFP'): #CHINA NESTED GRID FOR GEOS-FP
		lon = np.arange(70.0,140.01, 0.3125)
		lat = np.arange(15.0,55.01, 0.25)
	elif (gridlabel == 'EU_GEOSFP'): #EU NESTED GRID FOR GEOS-FP
		lon = np.arange(-15.0,40.01, 0.3125)
		lat = np.arange(32.75,61.26, 0.25)
	elif (gridlabel == 'NA_GEOSFP'): #NA NESTED GRID FOR GEOS-FP
		lon = np.arange(-130.0,-59.99, 0.3125)
		lat = np.arange(9.75,60.01, 0.25)
	else:
		raise ValueError('Scaling factor initialization utility does not recognize grid specification.')
	return [lon,lat]


def getLonLatEdgeFromLabel(gridlabel):
	if gridlabel == '4.0x5.0':
		lon = calc_rectilinear_lon_edge(5,True)
		lat = calc_rectilinear_lat_edge(4,True)
	elif gridlabel == '2.0x2.5':
		lon = calc_rectilinear_lon_edge(2.5,True)
		lat = calc_rectilinear_lat_edge(2,True)
	elif gridlabel == '1x1':
		lon = calc_rectilinear_lon_edge(1,False)
		lat = calc_rectilinear_lat_edge(1,False)
	elif gridlabel == '0.5x0.5':
		lon = calc_rectilinear_lon_edge(0.5,False)
		lat = calc_rectilinear_lat_edge(0.5,False)
	elif gridlabel == '0.25x0.25':
		lon = calc_rectilinear_lon_edge(0.25,False)
		lat = calc_rectilinear_lat_edge(0.25,False)
	elif gridlabel == '0.1x0.1':
		lon = calc_rectilinear_lon_edge(0.1,False)
		lat = calc_rectilinear_lat_edge(0.1,False)
	else:
		raise ValueError(f'Scaling factor initialization utility does not recognize grid specification {gridlabel}.')
	return [lon,lat]

def calc_rectilinear_lon_edge(lon_stride, center_at_180):
	n_lon = np.round(360.0 / lon_stride)
	lon_edge = np.linspace(-180.0, 180.0, num=int(n_lon + 1))
	if center_at_180:
		lon_edge = lon_edge - (lon_stride / 2.0)
	lon_edge[lon_edge < -180.0] = lon_edge[lon_edge < -180] + 360.0
	lon_edge[lon_edge > 180.0] = lon_edge[lon_edge > 180.0] - 360.0
	return lon_edge

def calc_rectilinear_lat_edge(lat_stride, half_polar_grid):
	if half_polar_grid:
		start_pt = 90.0 + (lat_stride / 2.0)
	else:
		start_pt = 90.0
	lat_edge = np.linspace(-1.0 * start_pt, start_pt,
						   num=1 + int(np.round(2.0 * start_pt / lat_stride)))
	# Force back onto +/- 90
	lat_edge[lat_edge > 90.0] = 90.0
	lat_edge[lat_edge < -90.0] = -90.0
	return lat_edge

def str2bool(v):
	if isinstance(v, bool):
		return v
	if v.lower() in ('yes', 'true', 't', 'y', '1'):
		return True
	elif v.lower() in ('no', 'false', 'f', 'n', '0'):
		return False
	else:
		raise argparse.ArgumentTypeError('Boolean value expected.')

