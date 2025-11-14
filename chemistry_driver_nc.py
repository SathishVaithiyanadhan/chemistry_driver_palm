import os
import re
import numpy as np
import netCDF4 as nc
from osgeo import gdal
from chemistry_driver_config import *
from chemistry_driver_main import read_geotiff_band, clear_geotiff_cache, process_species_emissions_sequential
from collections import defaultdict
import datetime 
from concurrent.futures import ThreadPoolExecutor
import multiprocessing as mp

# Pre-compiled regex for faster parsing - UPDATED for h00-h23 format
BAND_PATTERN = re.compile(r'^([A-Za-z_]+)_h(\d{2})_(\d{8})$')

def parse_band_description(desc):
    """Fast parsing with pre-compiled regex"""
    match = BAND_PATTERN.match(desc)
    return match.groups() if match else None

def get_time_bands_info_fast(geotiff_path):
    """Optimized band info extraction with date filtering - FIXED for h00-h23 format"""
    ds = gdal.Open(geotiff_path, gdal.GA_ReadOnly)
    if not ds: 
        raise ValueError(f"Can't open {geotiff_path}")
    
    # Use regular dict instead of defaultdict for serialization
    time_info = {}
    active_set = set(active_categories)
    
    for band_num in range(1, ds.RasterCount + 1):
        band = ds.GetRasterBand(band_num)
        desc = band.GetDescription()
        if desc:
            parsed = parse_band_description(desc)
            if parsed and parsed[0] in active_set:
                band_date_str = parsed[2]  # YYYYMMDD
                band_hour_str = parsed[1]  # hour as string '00' to '23'
                band_hour = int(band_hour_str)  # hour as integer 0-23
                
                # Convert band date to datetime for comparison
                band_dt = datetime.datetime.strptime(band_date_str, "%Y%m%d") + \
                         datetime.timedelta(hours=band_hour)
                
                # Check if band falls within the specified date range
                if start_dt <= band_dt <= end_dt:
                    # Initialize date dict if not exists
                    if band_date_str not in time_info:
                        time_info[band_date_str] = {}
                    # Initialize hour list if not exists
                    if f"h{band_hour_str}" not in time_info[band_date_str]:
                        time_info[band_date_str][f"h{band_hour_str}"] = []
                    
                    time_info[band_date_str][f"h{band_hour_str}"].append({
                        'band_num': band_num,
                        'category': parsed[0],
                        'hour_num': band_hour,  # 0-23
                        'hour_str': band_hour_str,  # '00'-'23'
                        'datetime': band_dt
                    })
    ds = None
    return time_info

def create_all_time_steps():
    """Create ALL expected time steps for the date range"""
    current_dt = start_dt
    time_steps = []
    
    while current_dt <= end_dt:
        date_str = current_dt.strftime("%Y%m%d")
        hour_num = current_dt.hour  # 0-23 (matches h00-h23 format)
        hour_str = f"h{hour_num:02d}"  # h00, h01, ..., h23
        
        time_steps.append({
            'date': date_str,
            'hour': hour_str,
            'hour_num': hour_num,  # 0-23
            'datetime': current_dt
        })
        
        # Move to next hour
        current_dt += datetime.timedelta(hours=1)
    
    return time_steps

def create_chemistry_driver(static_params):
    """Ultra-fast chemistry driver creation with optimized processing"""
    print("\nCreating PALM LOD2 Chemistry Driver")
    print(f"Date range: {start_date} to {end_date}")
    
    # Species mapping
    species_mapping = {
        'n2o': 'N2O', 'nox': 'NOX', 'nmvoc': 'RH', 'so2': 'SO2', 'co': 'CO',
        'pm10': 'PM10', 'pm2_5': 'PM25', 'nh3': 'NH3', 'pb': 'PB', 'cd': 'CD',
        'hg': 'HG', 'as': 'AS', 'ni': 'NI', 'bc': 'BC', 'co2': 'CO2', 'ch4': 'CH4',
        'no': 'NO', 'no2': 'NO2', 'ec': 'EC', 'na': 'NA', 'so4': 'SO4', 'oc': 'OC',
        'othmin': 'OTHMIN', 'o3': 'O3', 'hno3': 'HNO3', 'rcho': 'RCHO', 'ho2': 'HO2',
        'ro2': 'RO2', 'oh': 'OH', 'h2o': 'H2O'
    }
    
    # Pre-compute uppercase names
    uppercase_spec_names = []
    for spec in spec_name_str:
        if spec in species_mapping:
            uppercase_spec_names.append(species_mapping[spec])
        else:
            uppercase_spec_names.append(spec.upper().replace('_', ''))
    
    # Collect temporal emission data
    print("Scanning emission files...")
    all_time_info = {}
    valid_species = []
    
    for spec in spec_name_str:
        gt_path = f"{emis_geotiff_pth}emission_{spec}_temporal.tif"
        if os.path.exists(gt_path):
            print(f"  Found: {spec}")
            all_time_info[spec] = get_time_bands_info_fast(gt_path)
            valid_species.append(spec)
        else:
            print(f"  Missing: {spec}")
    
    if not all_time_info: 
        raise RuntimeError("No valid emission files found")

    # Create ALL expected time steps (24 hours for a full day)
    time_steps = create_all_time_steps()
    
    if not time_steps:
        raise RuntimeError(f"No time steps generated for the specified date range: {start_date} to {end_date}")

    print(f"Generated {len(time_steps)} time steps for date range")
    print(f"Time range: {time_steps[0]['datetime'].strftime('%Y-%m-%d %H:%M')} to {time_steps[-1]['datetime'].strftime('%Y-%m-%d %H:%M')}")

    current_time = datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S +000")
    output_file = f"{static_pth}{static}_chemistry"

    # Pre-compute timestamps
    print("Pre-computing timestamps...")
    timestamps_array = []
    for ts in time_steps:
        formatted_ts = ts['datetime'].strftime("%Y-%m-%d %H:%M:%S +000").ljust(64)
        timestamps_array.append(formatted_ts)

    # Create NetCDF file
    with nc.Dataset(output_file, 'w') as ds:
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
            'simulation_start': start_date,
            'simulation_end': end_date,
            'Conventions': 'CF-1.7'
        })

        # Dimensions
        n_time = len(time_steps)
        n_species = len(spec_name_str)
        
        ds.createDimension('z', static_params['nz'])
        ds.createDimension('field_length', 64)
        ds.createDimension('x', static_params['nx'])
        ds.createDimension('y', static_params['ny'])
        ds.createDimension('nspecies', n_species)
        ds.createDimension('time', None)

        # Coordinate variables
        z = ds.createVariable('z', 'f4', ('z',))
        z_data = np.linspace(static_params['z_origin'], 
                           static_params['z_origin'] + (static_params['nz']-1)*static_params['dz'],
                           static_params['nz'])
        z[:] = z_data
        z.units = "m"
        z.axis = "Z"
        z.long_name = "height above origin"

        x = ds.createVariable('x', 'f4', ('x',))
        x_data = np.linspace(static_params['origin_x'], 
                           static_params['origin_x'] + (static_params['nx']-1)*static_params['dx'],
                           static_params['nx'])
        x[:] = x_data
        x.units = "m"
        x.axis = "X"
        x.long_name = "x-distance from origin"

        y = ds.createVariable('y', 'f4', ('y',))
        y_data = np.linspace(static_params['south'], static_params['north'], static_params['ny'])
        y[:] = y_data
        y.axis = "Y"
        y.long_name = "y-distance from origin"
        y.units = "m"

        # nspecies variable
        nspecies_var = ds.createVariable('nspecies', 'i4', ('nspecies',))
        nspecies_var[:] = np.arange(1, n_species + 1, dtype=np.int32)
        nspecies_var.long_name = "nspecies"

        # Time variables
        time = ds.createVariable('time', 'i4', ('time',))
        time_data = np.arange(1, n_time + 1, dtype=np.int32)
        time[:] = time_data
        time.long_name = "time"
        time.standard_name = "time"
        time.units = "hours since first timestamp"

        timestamp = ds.createVariable('timestamp', 'S1', ('time', 'field_length'))
        timestamp_data = nc.stringtochar(np.array(timestamps_array, dtype='S64'))
        timestamp[:] = timestamp_data
        timestamp.long_name = "time stamp"

        # Emission metadata variables
        emission_name = ds.createVariable('emission_name', 'S1', 
                                         ('nspecies', 'field_length'))
        emission_name_data = nc.stringtochar(np.array(uppercase_spec_names, dtype='S64'))
        emission_name[:] = emission_name_data
        emission_name.long_name = "emission species name"
        emission_name.standard_name = "emission_name"

        emission_index = ds.createVariable('emission_index', 'f4', ('nspecies',),
                                         fill_value=np.float32(-9999.9))
        emission_index_data = np.arange(1, n_species + 1, dtype=np.float32)
        emission_index[:] = emission_index_data
        emission_index.long_name = "emission species index"
        emission_index.standard_name = "emission_index"

        # Main emission data
        emission_values = ds.createVariable('emission_values', 'f4', 
                                          ('time', 'z', 'y', 'x', 'nspecies'),
                                          fill_value=np.float32(-9999.9))
        emission_values.long_name = "emission species values"
        emission_values.standard_name = "emission_values"
        emission_values.units = "g/m2/s"
        emission_values.coordinates = "E_UTM N_UTM lon lat"
        emission_values.grid_mapping = "crsUTM: E_UTM N_UTM crsETRS: lon lat"
        emission_values.missing_value = np.float32(-9999.9)
        emission_values.lod = np.int32(2)

        # Stack height
        stack_height = ds.createVariable('emission_stack_height', 'f4', ('y', 'x'),
                                        fill_value=np.float32(-9999.9))
        building_mask = static_params['building_height'] > 0
        stack_height_data = np.full((static_params['ny'], static_params['nx']), 
                                  np.float32(-9999.9))
        stack_height_data[building_mask] = static_params['building_height'][building_mask]
        stack_height[:] = stack_height_data
        stack_height.long_name = "emission stack height"
        stack_height.standard_name = "emission_stack_height"
        stack_height.units = "m"
        stack_height.coordinates = "E_UTM N_UTM lon lat"
        stack_height.grid_mapping = "crsUTM: E_UTM N_UTM crsETRS: lon lat"
        stack_height.missing_value = np.float32(-9999.9)

        # Ultra-fast sequential emissions processing (still very fast with caching)
        print("Processing emissions data...")
        
        # Process species sequentially but with optimized caching
        for spec_idx, spec in enumerate(spec_name_str):
            spec_idx, species_emissions = process_species_emissions_sequential(
                spec, spec_idx, static_params, all_time_info, time_steps
            )
            emission_values[:, 0, :, :, spec_idx] = species_emissions

    # Clear cache
    clear_geotiff_cache()
    
    print(f"\nSUCCESS: Created {output_file}")
    print(f"- Species: {len(spec_name_str)} (nspecies starts from 1)")
    print(f"- Time steps: {len(time_steps)} (time starts from 1)")
    print(f"- Date range: {start_date} to {end_date}")
    print(f"- Emission indices: {list(range(1, len(spec_name_str) + 1))}")
    print(f"- Grid: {static_params['nx']}x{static_params['ny']}x{static_params['nz']}")
    print(f"- Units: g/m²/s")
    print(f"- Emission input for PALM LOD2 created successfully.")

# Replace the original function
create_chemistry_driver = create_chemistry_driver