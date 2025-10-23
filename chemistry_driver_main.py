import netCDF4 as nc
import numpy as np
from osgeo import gdal
import re
import os
from collections import defaultdict
from chemistry_driver_config import *

# Initialize global cache
_geotiff_cache = {}

def extract_static_parameters(static_file):
    """Extract domain parameters from static driver file"""
    with nc.Dataset(static_file, "r") as ncs:
        params = {}
        
        # Extract basic attributes
        params['origin_time'] = ncs.getncattr('origin_time')
        params['origin_lat'] = ncs.getncattr('origin_lat')
        params['origin_lon'] = ncs.getncattr('origin_lon')
        
        # Convert center point to UTM coordinates
        center_x, center_y = transformer_to_utm.transform(
            params['origin_lon'], params['origin_lat'])
        
        # Get grid dimensions and resolution
        params['nx'] = ncs.dimensions['x'].size
        params['ny'] = ncs.dimensions['y'].size
        params['nz'] = ncs.dimensions['z'].size if 'z' in ncs.dimensions else 1
        
        x_coords = ncs.variables['x'][:]
        y_coords = ncs.variables['y'][:]
        z_coords = ncs.variables['z'][:] if 'z' in ncs.variables else np.array([0])
        
        params['dx'] = float(x_coords[1] - x_coords[0]) if len(x_coords) > 1 else 1.0
        params['dy'] = float(abs(y_coords[1] - y_coords[0])) if len(y_coords) > 1 else 1.0
        params['dz'] = float(z_coords[1] - z_coords[0]) if len(z_coords) > 1 else 1.0
        params['z_origin'] = float(z_coords[0])
        
        # Read building height data
        if 'buildings_2d' in ncs.variables:
            params['building_height'] = ncs.variables['buildings_2d'][:, :].astype(np.float32)
        else:
            params['building_height'] = np.zeros((params['ny'], params['nx']), dtype=np.float32)
        
        # Calculate domain boundaries
        half_nx = (params['nx'] - 1) * params['dx'] / 2
        half_ny = (params['ny'] - 1) * params['dy'] / 2
        
        params['west'] = float(center_x - half_nx)
        params['east'] = float(center_x + half_nx)
        params['south'] = float(center_y - half_ny)
        params['north'] = float(center_y + half_ny)
        
        params['origin_x'] = params['west']
        params['origin_y'] = params['north']
        
        # Convert boundaries back to WGS84
        params['lon_w'], params['lat_s'] = transformer_to_wgs.transform(
            params['west'], params['south'])
        params['lon_e'], params['lat_n'] = transformer_to_wgs.transform(
            params['east'], params['north'])
            
    return params

def resample_entire_geotiff(geotiff_path, static_params):
    """Resample entire GeoTIFF once and cache it - MAJOR PERFORMANCE BOOST"""
    global _geotiff_cache
    
    if geotiff_path in _geotiff_cache:
        return _geotiff_cache[geotiff_path]
    
    if not os.path.exists(geotiff_path):
        raise FileNotFoundError(f"GeoTIFF file not found: {geotiff_path}")
    
    print(f"  Resampling entire file: {os.path.basename(geotiff_path)}")
    
    # Resample entire file to match static grid
    output_path = '/vsimem/resampled_entire.tif'
    
    # Single warp operation for entire file
    resampled_ds = gdal.Warp(
        output_path,
        geotiff_path,
        format='GTiff',
        outputBounds=[
            static_params['west'],
            static_params['south'],
            static_params['east'], 
            static_params['north']
        ],
        width=static_params['nx'],
        height=static_params['ny'],
        dstSRS=config_proj,
        resampleAlg=gdal.GRA_NearestNeighbour
    )
    
    if resampled_ds is None:
        raise RuntimeError(f"Failed to resample GeoTIFF: {geotiff_path}")
    
    # Read ALL bands at once into memory
    num_bands = resampled_ds.RasterCount
    all_bands_data = []
    
    for band_num in range(1, num_bands + 1):
        band = resampled_ds.GetRasterBand(band_num)
        arr = band.ReadAsArray().astype(np.float32)
        arr = np.flipud(arr)  # Match PALM's coordinate system
        all_bands_data.append(arr)
    
    # Cache the entire dataset
    _geotiff_cache[geotiff_path] = all_bands_data
    
    # Clean up
    resampled_ds = None
    
    return all_bands_data

def read_geotiff_band(geotiff_path, band_num, static_params):
    """reading bands - uses pre-resampled cached data"""
    global _geotiff_cache
    
    # Get or create cached resampled data
    if geotiff_path not in _geotiff_cache:
        all_bands = resample_entire_geotiff(geotiff_path, static_params)
    else:
        all_bands = _geotiff_cache[geotiff_path]
    
    # Return the specific band (band_num is 1-indexed in GDAL)
    if band_num - 1 < len(all_bands):
        return all_bands[band_num - 1].copy()
    else:
        raise ValueError(f"Band {band_num} not found in {geotiff_path}")

def clear_geotiff_cache():
    """Clear the cache to free memory"""
    global _geotiff_cache
    _geotiff_cache.clear()
    import gc
    gc.collect()

if __name__ == "__main__":
    print('\nExtracting static driver parameters')
    try:
        static_params = extract_static_parameters(static_pth + static + '_static')
    except Exception as e:
        raise RuntimeError(f"Failed to read static driver: {str(e)}")

    print("\nStatic Driver Configuration:")
    print(f"Grid: {static_params['nx']}x{static_params['ny']}x{static_params['nz']}")
    print(f"Resolution: {static_params['dx']}m x {static_params['dy']}m x {static_params['dz']}m")

    from chemistry_driver_nc import create_chemistry_driver
    create_chemistry_driver(static_params)