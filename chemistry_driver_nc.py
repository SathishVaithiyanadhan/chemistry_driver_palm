## Script for creating chemistry driver from PALM using GeoTIFF emissions
import os
import netCDF4 as nc
import numpy as np
import pandas as pd
from osgeo import gdal
from chemistry_driver_config import (static_pth, static, emis_geotiff_pth, 
                                   country, month, cat_name_str, cat_name, 
                                   spec_name, spec_name_str, default_monthly,
                                   default_daily, default_hourly)
from chemistry_driver_main import origin_time, origin_lat, origin_lon, origin_x, origin_y, nx, ny, res, dy, dx 

def read_geotiff_band(geotiff_path, band_num):
    """Read specific band from multi-band GeoTIFF"""
    ds = gdal.Open(geotiff_path)
    if ds is None:
        raise ValueError(f"Could not open GeoTIFF file: {geotiff_path}")
    band = ds.GetRasterBand(band_num)
    arr = band.ReadAsArray()
    ds = None  # Close the dataset
    return arr

print('Creating netCDF chemistry driver')

## NetCDF global attributes
ds = nc.Dataset(static_pth+static+'_chemistry', 'w', format='NETCDF4')
ds.lod = 1
ds.legacy_mode = "yes (z dimension enabled)"
ds.origin_time = origin_time
ds.origin_lat = origin_lat
ds.origin_lon = origin_lon
ds.origin_x = origin_x
ds.origin_y = origin_y
ds.resolution = res

## Add dimensions
ds.createDimension('z', 1)
ds.createDimension('y', ny)
ds.createDimension('x', nx)
ds.createDimension('nspecies', spec_name.size)
ds.createDimension('ncat', cat_name.size)
ds.createDimension('npm', 3)
ds.createDimension('nvoc', 2)
ds.createDimension('nmonthdayhour', 91)
ds.createDimension('max_string_length', 25)
ds.createDimension('nox_comp', 2)
ds.createDimension('sox_comp', 2)
ds.createDimension('pm_comp', 3)
ds.createDimension('voc_comp', 2)

## Add variables and Assign values
z = ds.createVariable('z', 'f8', ('z',))
z.units = 'm'
z.axis = 'Z'
z.long_name = 'distance to origin in z-direction'
z[:] = 1

x = ds.createVariable('x', 'f8', ('x',))
x.units = 'm'
x.axis = 'X'
x[:] = np.arange(10, (nx * dx + 1), dx)

y = ds.createVariable('y', 'f8', ('y',))
y.units = 'm'
y.axis = 'Y'
y[:] = np.arange(10, (ny * dy + 1), dy)

emission_name = ds.createVariable('emission_name', 'S1', ('nspecies', 'max_string_length'))
emission_name.long_name = "emission species name"
emission_name.standard_name = "emission name"
emission_name.units = ""
emission_name[:] = nc.stringtochar(spec_name)

emission_index = ds.createVariable('emission_index', 'u1', ('nspecies',), fill_value=-9999)
emission_index.long_name = "emission species index"
emission_index.standard_name = "emission index"
emission_index.units = ""
emission_index[:] = np.arange(1, len(emission_name) + 1)

emission_values = ds.createVariable('emission_values', 'f4', ('z', 'y', 'x', 'nspecies', 'ncat'), fill_value=-9999.9)
emission_values.long_name = "emission species values"
emission_values.standard_name = "emission values"
emission_values.lod = 1
emission_values.units = "kg/m2/year"
emission_values.coordinates = "E_UTM N_UTM lon lat"
emission_values.grid_mapping = "crsUTM: E_UTM N_UTM crsETRS: lon lat"

# Load emission data from GeoTIFFs
for nsp, spec in enumerate(spec_name_str):
    try:
        # Handle special case for PM2.5 filename
        if spec == 'pm2_5':
            geotiff_path = f"{emis_geotiff_pth}emission_{spec}_yearly.tif"
        else:
            geotiff_path = f"{emis_geotiff_pth}emission_{spec}_yearly.tif"
        
        print(f"Loading emissions from: {geotiff_path}")
        
        # Read each band (category) from the GeoTIFF
        for nct in range(len(cat_name_str)):
            try:
                # Band numbers start at 1
                chem_array = read_geotiff_band(geotiff_path, nct + 1)
                
                # Ensure array matches expected dimensions
                if chem_array.shape != (ny, nx):
                    raise ValueError(f"GeoTIFF dimensions {chem_array.shape} don't match expected ({ny}, {nx})")
                
                # Convert from g/m²/year to kg/m²/year and assign
                emission_values[:, :, :, nsp, nct] = (chem_array / 1000).reshape(1, ny, nx)
                
            except Exception as e:
                print(f"Warning: Could not read band {nct+1} ({cat_name_str[nct]}) from {geotiff_path}: {str(e)}")
                emission_values[:, :, :, nsp, nct] = 0.0
                
    except Exception as e:
        print(f"Error processing {spec} emissions: {str(e)}")
        emission_values[:, :, :, nsp, :] = 0.0

# Create simple time factors (equal distribution)
emission_time_factors = ds.createVariable('emission_time_factors', 'f4', ('nmonthdayhour', 'ncat'))
emission_time_factors.long_name = "emission_time_scaling_factors"
emission_time_factors.standard_name = "emission_time_scaling_factors"
emission_time_factors.lod = 1
emission_time_factors.units = ""

# Fill time factors with default values
time_factors = np.zeros((91, len(cat_name_str)))
time_factors[0:12, :] = default_monthly.reshape(12, 1)  # Monthly
time_factors[12:19, :] = default_daily.reshape(7, 1)    # Daily
time_factors[19:43, :] = default_hourly.reshape(24, 1)  # Weekday hours
time_factors[43:67, :] = default_hourly.reshape(24, 1)  # Saturday hours
time_factors[67:91, :] = default_hourly.reshape(24, 1)  # Sunday hours
emission_time_factors[:, :] = time_factors

# Emission category name
emission_category_name = ds.createVariable('emission_category_name', 'S1', ('ncat', 'max_string_length'))
emission_category_name.long_name = "emission category name"
emission_category_name.standard_name = "emission_cat_name"
emission_category_name.units = ""
emission_category_name[:] = nc.stringtochar(cat_name)

# Emission category index
emission_category_index = ds.createVariable('emission_cat_index', 'f4', ('ncat',), fill_value=-9999.9)
emission_category_index.long_name = "emission category index"
emission_category_index.standard_name = "emission_cat_index"
emission_category_index.units = ""
emission_category_index[:] = np.arange(1, len(cat_name_str) + 1)

# Composition variables (same as before)
composition_nox = ds.createVariable('composition_nox', 'f4', ('nox_comp', 'ncat'))
composition_nox.long_name = "composition of NOx"
composition_nox.standard_name = "composition_nox"
composition_nox.units = ""
composition_nox[0, :] = 0.8  # NO
composition_nox[1, :] = 0.2  # NO2

composition_sox = ds.createVariable('composition_sox', 'f4', ('sox_comp', 'ncat'))
composition_sox.long_name = "composition of SOx"
composition_sox.standard_name = "composition_sox"
composition_sox.units = ""
composition_sox[0, :] = 0.95  # SO2
composition_sox[1, :] = 0.05  # SO4

pm_name = np.array(['PM10', 'PM2.5', 'PM1'], dtype='S25')
emission_pm_name = ds.createVariable("emission_pm_name", 'S1', ('npm', 'max_string_length'))
emission_pm_name.long_name = "PM name"
emission_pm_name.standard_name = "pm_name"
emission_pm_name.units = ""
emission_pm_name[:] = nc.stringtochar(pm_name)

composition_pm = ds.createVariable("composition_pm", 'f4', ('pm_comp', 'npm'))
composition_pm.long_name = "composition of PM"
composition_pm.standard_name = "composition_PM"
composition_pm.units = ""
composition_pm[0, :] = [0.9, 0.09, 0.01]  # PM10, PM2.5, PM1 splits

voc_name = np.array(['NMVOC', 'MVOC'], dtype='S25')
emission_voc_name = ds.createVariable("emission_voc_name", 'S1', ('nvoc', 'max_string_length'))
emission_voc_name.long_name = "VOC name"
emission_voc_name.standard_name = "voc_name"
emission_voc_name.units = ""
emission_voc_name[:] = nc.stringtochar(voc_name)

composition_voc = ds.createVariable("composition_voc", 'f4', ('voc_comp', 'nvoc'))
composition_voc.long_name = "composition of VOC"
composition_voc.standard_name = "composition_VOC"
composition_voc.units = ""
composition_voc[0, :] = 0.9  # NMVOC
composition_voc[1, :] = 0.1  # MVOC

print('Chemistry driver created!')
ds.close()
##########
