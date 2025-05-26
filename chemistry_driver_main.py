#LOD2
import netCDF4 as nc
import numpy as np
from osgeo import gdal
from chemistry_driver_config import *
import re
import os
from collections import defaultdict

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
        
        params['dx'] = x_coords[1] - x_coords[0] if len(x_coords) > 1 else 1.0
        params['dy'] = abs(y_coords[1] - y_coords[0]) if len(y_coords) > 1 else 1.0
        params['dz'] = z_coords[1] - z_coords[0] if len(z_coords) > 1 else 1.0
        params['z_origin'] = z_coords[0]
        
        # Read building height data (2D array, no time dimension)
        if 'buildings_2d' in ncs.variables:
            params['building_height'] = ncs.variables['buildings_2d'][:, :]  # Correct 2D access
        else:
            params['building_height'] = np.zeros((params['ny'], params['nx']))
        
        # Calculate domain boundaries
        half_nx = (params['nx'] - 1) * params['dx'] / 2
        half_ny = (params['ny'] - 1) * params['dy'] / 2
        
        params['west'] = center_x - half_nx
        params['east'] = center_x + half_nx
        params['south'] = center_y - half_ny
        params['north'] = center_y + half_ny
        
        params['origin_x'] = params['west']
        params['origin_y'] = params['north']
        
        # Convert boundaries back to WGS84
        params['lon_w'], params['lat_s'] = transformer_to_wgs.transform(
            params['west'], params['south'])
        params['lon_e'], params['lat_n'] = transformer_to_wgs.transform(
            params['east'], params['north'])
            
    return params

def read_geotiff_band(geotiff_path, band_num, static_params):
    """Read and properly orient a single band from GeoTIFF"""
    if not os.path.exists(geotiff_path):
        raise FileNotFoundError(f"GeoTIFF file not found: {geotiff_path}")
    
    resampled_ds = resample_geotiff(geotiff_path, static_params)
    band = resampled_ds.GetRasterBand(band_num)
    arr = band.ReadAsArray()
    arr = np.flipud(arr)  # Match PALM's coordinate system
    
    if arr.shape != (static_params['ny'], static_params['nx']):
        raise ValueError(f"Resampled dimensions {arr.shape} don't match expected "
                       f"({static_params['ny']}, {static_params['nx']})")
    
    resampled_ds = None
    return arr

def resample_geotiff(input_path, static_params, output_path=None):
    """Resample GeoTIFF to exactly match static driver grid"""
    src_ds = gdal.Open(input_path)
    if src_ds is None:
        raise ValueError(f"Could not open GeoTIFF file: {input_path}")
    
    print("\nInput GeoTIFF Information:")
    src_gt = src_ds.GetGeoTransform()
    print(f"Dimensions: {src_ds.RasterXSize} x {src_ds.RasterYSize}")
    print(f"Original resolution: {src_gt[1]}m x {abs(src_gt[5])}m")
    print(f"Original CRS: {src_ds.GetProjection()}")
    
    if output_path is None:
        output_path = '/vsimem/temp_resampled.tif'
    
    resampled_ds = gdal.Warp(
        output_path,
        src_ds,
        format='GTiff',
        outputBounds=[
            static_params['west'],
            static_params['south'],
            static_params['east'],
            static_params['north']
        ],
        xRes=static_params['dx'],
        yRes=static_params['dy'],
        resampleAlg=gdal.GRA_NearestNeighbour,
        dstSRS=config_proj,
        targetAlignedPixels=True
    )
    
    if resampled_ds is None:
        raise RuntimeError("Failed to resample GeoTIFF")
    
    if (resampled_ds.RasterXSize != static_params['nx'] or 
        resampled_ds.RasterYSize != static_params['ny']):
        resampled_ds = gdal.Warp(
            '/vsimem/temp_resized.tif',
            resampled_ds,
            format='GTiff',
            width=static_params['nx'],
            height=static_params['ny'],
            resampleAlg=gdal.GRA_NearestNeighbour
        )
    
    print("\nResampled GeoTIFF Information:")
    print(f"Dimensions: {resampled_ds.RasterXSize} x {resampled_ds.RasterYSize}")
    print(f"New resolution: {static_params['dx']}m x {static_params['dy']}m")
    print(f"Output extent (UTM):")
    print(f"  West-East: {static_params['west']:.2f}m - {static_params['east']:.2f}m")
    print(f"  South-North: {static_params['south']:.2f}m - {static_params['north']:.2f}m")
    
    return resampled_ds

if __name__ == "__main__":
    print('\nExtracting static driver parameters')
    try:
        static_params = extract_static_parameters(static_pth + static + '_static')
    except Exception as e:
        raise RuntimeError(f"Failed to read static driver: {str(e)}")

    print("\nStatic Driver Configuration:")
    print(f"Center (WGS84): {static_params['origin_lat']:.6f}°N, {static_params['origin_lon']:.6f}°E")
    print(f"Center (UTM): {static_params['origin_x']:.2f}m, {static_params['origin_y']:.2f}m")
    print(f"Grid: {static_params['nx']}x{static_params['ny']}x{static_params['nz']}")
    print(f"Resolution: {static_params['dx']}m x {static_params['dy']}m x {static_params['dz']}m")
    print(f"Z levels: {static_params['nz']} from {static_params['z_origin']}m")

    from chemistry_driver_nc import create_chemistry_driver
    create_chemistry_driver(static_params)