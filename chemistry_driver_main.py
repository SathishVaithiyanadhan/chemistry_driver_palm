## Script for creating chemistry driver from PALM using emission GeoTIFFs
import netCDF4 as nc
import numpy as np
from chemistry_driver_config import static_pth, static

## Open static driver of simulation
print('Reading static driver')
ncs = nc.Dataset(static_pth+static+'_static', "r", format="NETCDF4")
origin_time = ncs.getncattr('origin_time')
origin_lat = ncs.getncattr('origin_lat')
origin_lon = ncs.getncattr('origin_lon')
origin_x = ncs.getncattr('origin_x')
origin_y = ncs.getncattr('origin_y')
nx = ncs.dimensions['x'].size
ny = ncs.dimensions['y'].size
res = float(10.0)  # Hardcoded resolution to match GeoTIFF processing
ncs.close()

xmax = origin_x + (nx * res)
ymax = origin_y + (ny * res)
dy = res
dx = res

## Create netCDF chemistry driver
import chemistry_driver_nc
