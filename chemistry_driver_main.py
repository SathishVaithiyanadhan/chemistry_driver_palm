## Script for creating chemistry driver from PALM using GRETA inventory
import netCDF4 as nc
import numpy as np
import pandas as pd

## Read configuration file
from chemistry_driver_config import static_pth, static

## Open static driver of simulation
print('Reading static driver')
ncs         = nc.Dataset(static_pth+static+'_static', "r", format="NETCDF4")
origin_time = ncs.getncattr('origin_time')
origin_lat  = ncs.getncattr('origin_lat')
origin_lon  = ncs.getncattr('origin_lon')
origin_x    = ncs.getncattr('origin_x')
origin_y    = ncs.getncattr('origin_y')
nx          = ncs.dimensions['x'].size
ny          = ncs.dimensions['y'].size
#res         = float(ncs.getncattr('resolution'))
res =float(10.0)
ncs.close()
xmax = origin_x+(nx*res)
ymax = origin_y+(ny*res)
dy = res
dx = res

## Resample and resize emission inventory
## Improved scripts with sectors to follow
print('Resample and resize emission inventory')
#import chemistry_driver_sum
#import chemistry_driver_transport
#import chemistry_driver_sectors

## Create netCDF chenistry driver
import chemistry_driver_nc
