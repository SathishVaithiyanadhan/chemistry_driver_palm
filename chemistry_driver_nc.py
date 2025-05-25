import os
import re
import numpy as np
import netCDF4 as nc
from osgeo import gdal
from chemistry_driver_config import *
from chemistry_driver_main import resample_geotiff, read_geotiff_band
from collections import defaultdict

def parse_band_description(desc):
    """Parse band description into category, hour, and date components"""
    # Pattern to match: Category_hour_YYYYMMDD
    pattern = r'^([A-Za-z_]+)_h(\d+)_(\d{8})$'
    match = re.match(pattern, desc)
    if match:
        return {
            'category': match.group(1),
            'hour': f"h{match.group(2)}",
            'date': match.group(3),
            'hour_num': int(match.group(2))
        }
    return None

def get_time_bands_info(geotiff_path):
    """Extract complete time and category information from all bands"""
    ds = gdal.Open(geotiff_path)
    if ds is None:
        raise ValueError(f"Could not open GeoTIFF file: {geotiff_path}")
    
    time_info = defaultdict(lambda: defaultdict(list))
    
    for band_num in range(1, ds.RasterCount + 1):
        band = ds.GetRasterBand(band_num)
        desc = band.GetDescription()
        if not desc:
            continue
            
        parsed = parse_band_description(desc)
        if parsed:
            time_info[parsed['date']][parsed['hour']].append({
                'band_num': band_num,
                'category': parsed['category'],
                'desc': desc,
                'hour_num': parsed['hour_num']
            })
    
    ds = None
    return time_info

def create_time_dimensions(time_info):
    """Create sorted time dimensions from the collected time info"""
    # Get all unique dates and hours
    all_dates = sorted(time_info.keys())
    all_hours = set()
    
    for date in all_dates:
        all_hours.update(time_info[date].keys())
    
    # Sort hours numerically
    sorted_hours = sorted(all_hours, key=lambda x: int(x[1:]))
    
    # Create time steps combining date and hour
    time_steps = []
    for date in all_dates:
        for hour in sorted_hours:
            if hour in time_info[date]:
                time_steps.append({
                    'date': date,
                    'hour': hour,
                    'hour_num': int(hour[1:])
                })
    
    return time_steps

def create_chemistry_driver(static_params):
    """Main function to create the chemistry driver"""
    print("\nStarting chemistry driver creation process")
    
    # First collect all time information from all files
    all_time_info = {}
    
    for spec in spec_name_str:
        geotiff_path = f"{emis_geotiff_pth}emission_{spec}_temporal.tif"
        if not os.path.exists(geotiff_path):
            print(f"Warning: Emission file not found: {geotiff_path}")
            continue
            
        print(f"\nProcessing emissions for: {spec}")
        time_info = get_time_bands_info(geotiff_path)
        all_time_info[spec] = time_info
    
    if not all_time_info:
        raise RuntimeError("No valid time steps found in any input files")
    
    # Use the first species to determine time structure (assuming all have same structure)
    ref_spec = next(iter(all_time_info.keys()))
    time_steps = create_time_dimensions(all_time_info[ref_spec])
    
    print(f"\nFound {len(time_steps)} time steps across {len(all_time_info[ref_spec])} days")
    
    # Create NetCDF file
    ds = nc.Dataset(static_pth + static + '_chemistry', 'w', format='NETCDF4')
    
    # Global attributes
    ds.lod = 2
    ds.legacy_mode = "yes (z dimension enabled)"
    ds.origin_time = static_params['origin_time']
    ds.origin_lat = static_params['origin_lat']
    ds.origin_lon = static_params['origin_lon']
    ds.origin_x = static_params['origin_x']
    ds.origin_y = static_params['origin_y']
    ds.resolution = static_params['dx']
    
    # Create dimensions
    ds.createDimension('z', static_params['nz'])
    ds.createDimension('y', static_params['ny'])
    ds.createDimension('x', static_params['nx'])
    ds.createDimension('nspecies', len(spec_name_str))
    ds.createDimension('ncat', len(cat_name_str))
    ds.createDimension('ntime', len(time_steps))
    ds.createDimension('max_string_length', 25)
    
    # Coordinate variables
    z = ds.createVariable('z', 'f4', ('z',))
    z.units = 'm'
    z.axis = 'Z'
    z.long_name = 'height above origin'
    z[:] = np.linspace(
        static_params['z_origin'],
        static_params['z_origin'] + (static_params['nz'] - 1) * static_params['dz'],
        static_params['nz']
    )
    
    x = ds.createVariable('x', 'f4', ('x',))
    x.units = 'm'
    x.axis = 'X'
    x.long_name = 'distance to origin in x-direction'
    x[:] = np.linspace(
        static_params['origin_x'],
        static_params['origin_x'] + (static_params['nx'] - 1) * static_params['dx'],
        static_params['nx']
    )
    
    y = ds.createVariable('y', 'f4', ('y',))
    y.units = 'm'
    y.axis = 'Y'
    y.long_name = 'distance to origin in y-direction'
    y[:] = np.linspace(
        static_params['south'], 
        static_params['north'],  
        static_params['ny']
    )
    
    # Time variables
    time = ds.createVariable('time', 'f4', ('ntime',))
    time.units = 'hours since 00:00:00'
    time.long_name = 'time'
    
    time_name = ds.createVariable('time_name', 'S1', ('ntime', 'max_string_length'))
    time_name.long_name = "time step names"
    
    # Initialize time variables
    for i, ts in enumerate(time_steps):
        time[i] = ts['hour_num']  # Store hour number
        time_name[i] = nc.stringtochar(np.array([f"{ts['date']}_{ts['hour']}".ljust(25)], dtype='S25'))
    
    # Emission variables
    emission_name = ds.createVariable('emission_name', 'S1', ('nspecies', 'max_string_length'))
    emission_name.long_name = "emission species name"
    emission_name[:] = nc.stringtochar(np.array(spec_name_str, dtype='S25'))
    
    emission_values = ds.createVariable('emission_values', 'f4', 
                                      ('ntime', 'z', 'y', 'x', 'nspecies', 'ncat'), 
                                      fill_value=0.0)
    emission_values.long_name = "emission species values"
    emission_values.units = "kg/m2/hour"
    
    emission_category_name = ds.createVariable('emission_category_name', 'S1', 
                                             ('ncat', 'max_string_length'))
    emission_category_name.long_name = "emission category name"
    emission_category_name[:] = nc.stringtochar(np.array(cat_name_str, dtype='S25'))
    
    # Process all emission data
    for nsp, spec in enumerate(spec_name_str):
        if spec not in all_time_info:
            continue
            
        print(f"\nProcessing emissions for: {spec}")
        time_info = all_time_info[spec]
        
        # Process each time step
        for time_idx, ts in enumerate(time_steps):
            date = ts['date']
            hour = ts['hour']
            
            if date not in time_info or hour not in time_info[date]:
                continue
                
            # Get all bands for this time step
            bands = time_info[date][hour]
            
            # Process each category
            for band_info in bands:
                try:
                    cat = band_info['category']
                    if cat not in cat_name_str:
                        continue
                        
                    nct = list(cat_name_str).index(cat)
                    
                    # Read and store the data
                    chem_array = read_geotiff_band(
                        f"{emis_geotiff_pth}emission_{spec}_temporal.tif", 
                        band_info['band_num'], 
                        static_params
                    )
                    emission_values[time_idx, 0, :, :, nsp, nct] = chem_array
                    
                except Exception as e:
                    print(f"Error processing {spec} {date} {hour} {cat}: {str(e)}")
    
    print("\nChemistry driver created successfully!")
    ds.close()