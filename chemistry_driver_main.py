import netCDF4 as nc
import numpy as np
from osgeo import gdal
import os
from chemistry_driver_config import *

# Initialize global cache with LRU behavior
_geotiff_cache = {}
_cache_order = []
_MAX_CACHE_SIZE = 5  # Keep only 5 files in memory

def filter_negative_and_zero_values(emission_array):
    """
    Set all values ≤ 0 to NaN in the emission array.
    This ensures only positive emission values are kept.
    """
    # Create a copy to avoid modifying the original
    filtered_array = emission_array.copy()
    
    # Set values ≤ 0 to NaN
    filtered_array[filtered_array <= 0] = np.float32(-9999.9)
    
    return filtered_array

def extract_static_parameters(static_file):
    """Ultra-fast static parameter extraction with minimal operations"""
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
    """Ultra-fast GeoTIFF resampling with optimized GDAL operations"""
    global _geotiff_cache, _cache_order
    
    if geotiff_path in _geotiff_cache:
        # Move to front of LRU
        _cache_order.remove(geotiff_path)
        _cache_order.insert(0, geotiff_path)
        return _geotiff_cache[geotiff_path]
    
    if not os.path.exists(geotiff_path):
        raise FileNotFoundError(f"GeoTIFF file not found: {geotiff_path}")
    
    print(f"  Resampling: {os.path.basename(geotiff_path)}")
    
    # Optimized warp options
    warp_options = gdal.WarpOptions(
        format='MEM',
        outputBounds=[static_params['west'], static_params['south'], 
                     static_params['east'], static_params['north']],
        width=static_params['nx'],
        height=static_params['ny'],
        dstSRS=config_proj,
        resampleAlg=gdal.GRA_NearestNeighbour
    )
    
    resampled_ds = gdal.Warp('', geotiff_path, options=warp_options)
    
    if resampled_ds is None:
        raise RuntimeError(f"Failed to resample GeoTIFF: {geotiff_path}")
    
    # Read all bands efficiently
    num_bands = resampled_ds.RasterCount
    all_bands_data = []
    
    for band_num in range(1, num_bands + 1):
        band = resampled_ds.GetRasterBand(band_num)
        arr = band.ReadAsArray().astype(np.float32)
        arr = np.flipud(arr)  # Match PALM's coordinate system
        all_bands_data.append(arr)
    
    # Cache management with LRU
    if len(_cache_order) >= _MAX_CACHE_SIZE:
        oldest = _cache_order.pop()
        del _geotiff_cache[oldest]
    
    _geotiff_cache[geotiff_path] = all_bands_data
    _cache_order.insert(0, geotiff_path)
    
    # Clean up
    resampled_ds = None
    
    return all_bands_data

def read_geotiff_band(geotiff_path, band_num, static_params):
    """Ultra-fast band reading with cache optimization"""
    global _geotiff_cache
    
    # Get or create cached resampled data
    if geotiff_path not in _geotiff_cache:
        all_bands = resample_entire_geotiff(geotiff_path, static_params)
    else:
        all_bands = _geotiff_cache[geotiff_path]
    
    # Return the specific band (band_num is 1-indexed in GDAL)
    if band_num - 1 < len(all_bands):
        return all_bands[band_num - 1]
    else:
        raise ValueError(f"Band {band_num} not found in {geotiff_path}")

def clear_geotiff_cache():
    """Clear the cache to free memory"""
    global _geotiff_cache, _cache_order
    _geotiff_cache.clear()
    _cache_order.clear()
    import gc
    gc.collect()

def process_species_emissions(spec, spec_idx, static_params, all_time_info, time_steps):
    """Process emissions for a single species with zero/negative filtering"""
    print(f"  Processing {spec}...")
    species_emissions = np.full((len(time_steps), static_params['ny'], static_params['nx']), 
                               np.float32(-9999.9))
    
    if spec not in all_time_info:
        return spec_idx, species_emissions
    
    time_info = all_time_info[spec]
    
    for ts_idx, ts in enumerate(time_steps):
        date_key = ts['date']
        hour_key = ts['hour']
        
        if date_key not in time_info or hour_key not in time_info[date_key]:
            continue
        
        # Initialize with fill values
        total_emission = np.full((static_params['ny'], static_params['nx']), 
                               np.float32(-9999.9))
        
        # Aggregate emissions
        bands = time_info[date_key][hour_key]
        for band in bands:
            arr = read_geotiff_band(
                f"{emis_geotiff_pth}emission_{spec}_temporal.tif",
                band['band_num'],
                static_params
            )
            
            # Vectorized aggregation with unit conversion
            valid_mask = ~np.isnan(arr)
            fill_mask = (total_emission == np.float32(-9999.9)) & valid_mask
            add_mask = valid_mask & ~fill_mask
            
            # Convert from kg/m2/hour to g/m2/s
            arr_g_per_sec = arr * (1000.0 / 3600.0)
            
            total_emission[fill_mask] = arr_g_per_sec[fill_mask]
            total_emission[add_mask] += arr_g_per_sec[add_mask]
        
        # Apply filtering: Set all values ≤ 0 to NaN
        total_emission = filter_negative_and_zero_values(total_emission)
        
        species_emissions[ts_idx] = total_emission
    
    return spec_idx, species_emissions

if __name__ == "__main__":
    print('\nExtracting static driver parameters')
    try:
        static_params = extract_static_parameters(static_pth + static + '_static')
    except Exception as e:
        raise RuntimeError(f"Failed to read static driver: {str(e)}")
    
    print("\nStatic Driver Configuration:")
    print(f"Grid: {static_params['nx']}x{static_params['ny']}x{static_params['nz']}")
    print(f"Resolution: {static_params['dx']}m x {static_params['dy']}m x {static_params['dz']}m")
    print(f"Date range: {start_date} to {end_date}")
    print(f"Traffic tag: {'ENABLED' if tag == 'traffic' else 'DISABLED'}")
    print(f"Zero/negative filtering: ENABLED (all values ≤ 0 will be set to NaN)")
    
    from chemistry_driver_nc import create_chemistry_driver
    create_chemistry_driver(static_params)