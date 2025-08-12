import os
import re
import numpy as np
import netCDF4 as nc
from osgeo import gdal
from chemistry_driver_config import *
from chemistry_driver_main import resample_geotiff, read_geotiff_band
from collections import defaultdict
import datetime 

def parse_band_description(desc):
    pattern = r'^([A-Za-z_]+)_h(\d+)_(\d{8})$'
    match = re.match(pattern, desc)
    return match.groups() if match else None

def get_time_bands_info(geotiff_path):
    ds = gdal.Open(geotiff_path)
    if not ds: raise ValueError(f"Can't open {geotiff_path}")
    
    time_info = defaultdict(lambda: defaultdict(list))
    for band_num in range(1, ds.RasterCount + 1):
        band = ds.GetRasterBand(band_num)
        if desc := band.GetDescription():
            if (parsed := parse_band_description(desc)) and parsed[0] in active_categories:
                time_info[parsed[2]][f"h{parsed[1]}"].append({
                    'band_num': band_num,
                    'category': parsed[0],
                    'hour_num': int(parsed[1])
                })
    ds = None
    return time_info

def create_time_dimensions(time_info):
    dates = sorted(time_info.keys())
    hours = sorted({h for d in dates for h in time_info[d]}, key=lambda x: int(x[1:]))
    return [{'date': d, 'hour': h, 'hour_num': int(h[1:])} 
            for d in dates for h in hours if h in time_info[d]]

def create_chemistry_driver(static_params):
    print("\nCreating PALM LOD2 Chemistry Driver")
    
    # Define species name mapping to uppercase format
    species_mapping = {
        'n2o': 'N2O',
        'nox': 'NOX',
        'nmvoc': 'NMVOC',
        'so2': 'SO2',
        'co': 'CO',
        'pm10': 'PM10',
        'pm2_5': 'PM25',
        'nh3': 'NH3',
        'pb': 'PB',
        'cd': 'CD',
        'hg': 'HG',
        'as': 'AS',
        'ni': 'NI',
        'bc': 'BC',
        'co2': 'CO2',
        'ch4': 'CH4',
        'no': 'NO',
        'no2': 'NO2',
        'ec': 'EC',
        'na': 'NA',
        'so4': 'SO4',
        'oc': 'OC',
        'othmin': 'OTHMIN',
        'o3' : 'O3'
    }
    
    # Convert species names to uppercase format
    uppercase_spec_names = []
    for spec in spec_name_str:
        if spec in species_mapping:
            uppercase_spec_names.append(species_mapping[spec])
        else:
            uppercase_spec_names.append(spec.upper().replace('_', ''))
            
    # Collect temporal emission data
    all_time_info = {}
    for spec in spec_name_str:
        gt_path = f"{emis_geotiff_pth}emission_{spec}_temporal.tif"
        if os.path.exists(gt_path):
            all_time_info[spec] = get_time_bands_info(gt_path)
    
    if not all_time_info: 
        raise RuntimeError("No valid emission files found")

    # Collect all unique time steps to determine end_time
    all_time_steps = []
    for spec_info in all_time_info.values():
        for date_str in spec_info:
            for hour_str in spec_info[date_str]:
                hour_num = int(hour_str[1:])  # Convert 'hHH' to integer
                dt = datetime.datetime.strptime(date_str, "%Y%m%d") + datetime.timedelta(hours=hour_num)
                all_time_steps.append(dt)
    
    # Remove duplicates, sort, and get latest time
    unique_dts = sorted(list(set(all_time_steps))) if all_time_steps else []
    end_time_str = unique_dts[-1].strftime("%Y-%m-%d %H:%M:%S +000") if unique_dts else "N/A"
    current_time = datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S +000")

    # Create NetCDF dataset
    with nc.Dataset(f"{static_pth}{static}_chemistry", 'w') as ds:
        # Global attributes
        ds.setncatts({
            'description': 'Chemistry driver for PALM model simulation with LOD2 emissions from the Downscaled GRETA Emissions',
            'author': 'Sathish Kumar Vaithiyanadhan (sathish.vaithiyanadhan@med.uni-augsburg.de)',
            'Institution': 'Chair of Model-based Environmental Exposure Science (MBEES), Faculty of Medicine, University of Augsburg, Germany.',
            'lod': 2,
            'origin_time': static_params['origin_time'],
            'origin_lat': static_params['origin_lat'],
            'origin_lon': static_params['origin_lon'],
            'origin_x': static_params['origin_x'],
            'origin_y': static_params['origin_y'],
            'resolution': static_params['dx'],
            'history': f'Created on {current_time}',
            'source': 'multiple sources; netCDF5',
            'origin_z': static_params['z_origin'],
            'Conventions': 'CF-1.7'
        })

        # Dimensions
        ds.createDimension('time', None)
        ds.createDimension('z', static_params['nz'])
        ds.createDimension('y', static_params['ny'])
        ds.createDimension('x', static_params['nx'])
        ds.createDimension('nspecies', len(spec_name_str))
        ds.createDimension('field_length', 64)  # Standard PALM character length

        # Coordinate variables
        z = ds.createVariable('z', 'f4', ('z',))
        z[:] = np.linspace(static_params['z_origin'], 
                          static_params['z_origin'] + (static_params['nz']-1)*static_params['dz'],
                          static_params['nz'])
        z.setncatts({'units': 'm', 'axis': 'Z', 'long_name': 'height above origin'})

        x = ds.createVariable('x', 'f4', ('x',))
        x[:] = np.linspace(static_params['origin_x'], 
                          static_params['origin_x'] + (static_params['nx']-1)*static_params['dx'],
                          static_params['nx'])
        x.setncatts({'units': 'm', 'axis': 'X', 'long_name': 'x-distance from origin'})

        y = ds.createVariable('y', 'f4', ('y',))
        y[:] = np.linspace(static_params['south'], static_params['north'], static_params['ny'])
        y.setncatts({'units': 'm', 'axis': 'Y', 'long_name': 'y-distance from origin'})

        nspecies_var = ds.createVariable('nspecies', 'i4', ('nspecies',))
        nspecies_var[:] = np.arange(1, len(spec_name_str) + 1)  # Start from 1
        nspecies_var.long_name = "nspecies"

        # Time variables
        time = ds.createVariable('time', 'i4', ('time',))
        time.setncatts({
            'long_name': 'time',
            'standard_name': 'time',
            'units': 'hours since first timestamp'
        })

        timestamp = ds.createVariable('timestamp', 'S1', ('time', 'field_length'))
        timestamp.long_name = "time stamp"

        # Emission metadata variables
        emission_name = ds.createVariable('emission_name', 'S1', 
                                         ('nspecies', 'field_length'))
        emission_name.long_name = "emission species name"
        emission_name.standard_name = "emission_name"
        emission_name[:] = nc.stringtochar(np.array(uppercase_spec_names, dtype='S64'))

        emission_index = ds.createVariable('emission_index', 'f4', ('nspecies',),
                                         fill_value=-9999.9)
        emission_index.long_name = "emission species index"
        emission_index.standard_name = "emission_index"
        emission_index[:] = np.arange(1, len(spec_name_str) + 1)  # Start from 1

        # Main emission data with NaN handling
        emission_values = ds.createVariable('emission_values', 'f4', 
                                          ('time', 'z', 'y', 'x', 'nspecies'),
                                          fill_value=np.float32(-9999.9))
        emission_values.setncatts({
            'long_name': 'emission species values',
            'standard_name': 'emission_values',
            'units': 'kg/m2/hour',
            'coordinates': 'E_UTM N_UTM lon lat',
            'grid_mapping': 'crsUTM: E_UTM N_UTM crsETRS: lon lat',
            'missing_value': np.float32(-9999.9)
        })

        # Stack height (building-based)
        stack_height = ds.createVariable('emission_stack_height', 'f4', ('y', 'x'),
                                        fill_value=np.float32(-9999.9))
        stack_height.setncatts({
            'long_name': 'emission stack height',
            'standard_name': 'emission_stack_height',
            'units': 'm',
            'coordinates': 'E_UTM N_UTM lon lat',
            'grid_mapping': 'crsUTM: E_UTM N_UTM crsETRS: lon lat',
            'missing_value': np.float32(-9999.9)
        })
        stack_height[:] = np.where(static_params['building_height'] > 0,
                                 static_params['building_height'], np.float32(-9999.9))
        
        for var in ds.variables:
            if "_ChunkSizes" in ds.variables[var].__dict__:
               del ds.variables[var]._ChunkSizes

        # Process emissions data with NaN preservation
        for spec_idx, spec in enumerate(spec_name_str):
            if spec not in all_time_info: continue
            
            print(f"Processing {spec} emissions...")
            time_info = all_time_info[spec]
            time_steps = create_time_dimensions(time_info)
            
            for ts_idx, ts in enumerate(time_steps):
                # Time index starts at 1 (ts_idx + 1)
                time[ts_idx] = ts_idx + 1
                
                # Create timestamp - subtract 1 hour from the original time to shift from 01:00:00 to 00:00:00
                dt = datetime.datetime.strptime(ts['date'], "%Y%m%d") + datetime.timedelta(hours=ts['hour_num'] - 1)
                formatted_timestamp = dt.strftime("%Y-%m-%d %H:%M:%S +000")
                timestamp[ts_idx] = nc.stringtochar(np.array(formatted_timestamp.ljust(64), dtype='S64'))
                
                # Initialize output array with fill values
                total_emission = np.full((static_params['ny'], static_params['nx']), np.float32(-9999.9))
                
                # Aggregate emissions across categories
                for band in time_info[ts['date']][ts['hour']]:
                    arr = read_geotiff_band(
                        f"{emis_geotiff_pth}emission_{spec}_temporal.tif",
                        band['band_num'],
                        static_params
                    )
                    # Preserve NaN values from input
                    valid_mask = ~np.isnan(arr)
                    total_emission[valid_mask] = np.where(
                        total_emission[valid_mask] == np.float32(-9999.9),
                        arr[valid_mask],
                        total_emission[valid_mask] + arr[valid_mask]
                    )
                
                emission_values[ts_idx, 0, :, :, spec_idx] = total_emission

    print("\nSuccessfully created PALM chemistry driver with:")
    print(f"- {len(spec_name_str)} species")
    print(f"- {len(time_steps)} time steps")
    print(f"- {static_params['nx']}x{static_params['ny']} grid")