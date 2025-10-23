import os
import re
import numpy as np
import netCDF4 as nc
from osgeo import gdal
from chemistry_driver_config import *
from chemistry_driver_main import read_geotiff_band, clear_geotiff_cache
from collections import defaultdict
import datetime 

# Pre-compiled regex for faster parsing
BAND_PATTERN = re.compile(r'^([A-Za-z_]+)_h(\d+)_(\d{8})$')

def parse_band_description(desc):
    """Fast parsing with pre-compiled regex"""
    match = BAND_PATTERN.match(desc)
    return match.groups() if match else None

def get_time_bands_info_fast(geotiff_path):
    """Optimized band info extraction"""
    ds = gdal.Open(geotiff_path, gdal.GA_ReadOnly)
    if not ds: 
        raise ValueError(f"Can't open {geotiff_path}")
    
    time_info = defaultdict(lambda: defaultdict(list))
    active_set = set(active_categories)
    
    for band_num in range(1, ds.RasterCount + 1):
        band = ds.GetRasterBand(band_num)
        desc = band.GetDescription()
        if desc:
            parsed = parse_band_description(desc)
            if parsed and parsed[0] in active_set:
                time_info[parsed[2]][f"h{parsed[1]}"].append({
                    'band_num': band_num,
                    'category': parsed[0],
                    'hour_num': int(parsed[1])
                })
    ds = None
    return time_info

def create_time_dimensions_fast(time_info):
    """Optimized time dimension creation"""
    dates = sorted(time_info.keys())
    if not dates:
        return []
    
    hours = set()
    for date_info in time_info.values():
        hours.update(date_info.keys())
    
    hours = sorted(hours, key=lambda x: int(x[1:]))
    
    time_steps = []
    for d in dates:
        for h in hours:
            if h in time_info[d]:
                time_steps.append({
                    'date': d, 
                    'hour': h, 
                    'hour_num': int(h[1:])
                })
    return time_steps

def create_chemistry_driver(static_params):
    """chemistry driver creation - EXACT MATCH TO SAMPLE"""
    print("\nCreating PALM LOD2 Chemistry Driver")
    
    # Species mapping - exactly as in sample
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

    # Determine time structure
    first_spec = valid_species[0]
    time_steps = create_time_dimensions_fast(all_time_info[first_spec])
    
    if not time_steps:
        raise RuntimeError("No valid time steps found in emission data")

    print(f"Found {len(time_steps)} time steps across {len(valid_species)} species")

    current_time = datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S +000")
    output_file = f"{static_pth}{static}_chemistry"

    # Pre-compute timestamps
    print("Pre-computing timestamps...")
    timestamps_array = []
    for ts in time_steps:
        dt = datetime.datetime.strptime(ts['date'], "%Y%m%d") + \
             datetime.timedelta(hours=ts['hour_num'] - 1)
        formatted_ts = dt.strftime("%Y-%m-%d %H:%M:%S +000").ljust(64)
        timestamps_array.append(formatted_ts)

    # Create NetCDF with EXACT structure matching sample
    with nc.Dataset(output_file, 'w') as ds:
        # Global attributes - EXACT MATCH
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

        # Dimensions - EXACT MATCH
        n_time = len(time_steps)
        n_species = len(spec_name_str)
        
        ds.createDimension('z', static_params['nz'])
        ds.createDimension('field_length', 64)
        ds.createDimension('x', static_params['nx'])
        ds.createDimension('y', static_params['ny'])
        ds.createDimension('nspecies', n_species)
        ds.createDimension('time', None)  # UNLIMITED as in sample

        # Coordinate variables - EXACT MATCH
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

        # nspecies variable - EXACT MATCH (starts from 1)
        nspecies_var = ds.createVariable('nspecies', 'i4', ('nspecies',))
        nspecies_var[:] = np.arange(1, n_species + 1, dtype=np.int32)  # Starts from 1
        nspecies_var.long_name = "nspecies"

        # Time variables - EXACT MATCH
        time = ds.createVariable('time', 'i4', ('time',))
        time_data = np.arange(1, n_time + 1, dtype=np.int32)  # Starts from 1
        time[:] = time_data
        time.long_name = "time"
        time.standard_name = "time"
        time.units = "hours since first timestamp"

        timestamp = ds.createVariable('timestamp', 'S1', ('time', 'field_length'))
        timestamp_data = nc.stringtochar(np.array(timestamps_array, dtype='S64'))
        timestamp[:] = timestamp_data
        timestamp.long_name = "time stamp"

        # Emission metadata variables - EXACT MATCH (with fill_value set at creation)
        emission_name = ds.createVariable('emission_name', 'S1', 
                                         ('nspecies', 'field_length'))
        emission_name_data = nc.stringtochar(np.array(uppercase_spec_names, dtype='S64'))
        emission_name[:] = emission_name_data
        emission_name.long_name = "emission species name"
        emission_name.standard_name = "emission_name"

        # emission_index with fill_value set at creation
        emission_index = ds.createVariable('emission_index', 'f4', ('nspecies',),
                                         fill_value=np.float32(-9999.9))
        emission_index_data = np.arange(1, n_species + 1, dtype=np.float32)  # Starts from 1, float as in sample
        emission_index[:] = emission_index_data
        emission_index.long_name = "emission species index"
        emission_index.standard_name = "emission_index"

        # Main emission data with fill_value set at creation
        emission_values = ds.createVariable('emission_values', 'f4', 
                                          ('time', 'z', 'y', 'x', 'nspecies'),
                                          fill_value=np.float32(-9999.9))
        emission_values.long_name = "emission species values"
        emission_values.standard_name = "emission_values"
        emission_values.units = "kg/m2/hour"
        emission_values.coordinates = "E_UTM N_UTM lon lat"
        emission_values.grid_mapping = "crsUTM: E_UTM N_UTM crsETRS: lon lat"
        emission_values.missing_value = np.float32(-9999.9)

        # Stack height with fill_value set at creation
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

        # emissions processing
        print("Processing emissions data...")
        
        for spec_idx, spec in enumerate(spec_name_str):
            if spec not in all_time_info: 
                # Fill with missing values for all time steps
                emission_values[:, 0, :, :, spec_idx] = np.float32(-9999.9)
                continue
                
            print(f"  Processing {spec}...")
            time_info = all_time_info[spec]
            
            for ts_idx, ts in enumerate(time_steps):
                date_key = ts['date']
                hour_key = ts['hour']
                
                if date_key not in time_info or hour_key not in time_info[date_key]:
                    # No data for this time step
                    emission_values[ts_idx, 0, :, :, spec_idx] = np.float32(-9999.9)
                    continue
                
                # Initialize with fill values
                total_emission = np.full((static_params['ny'], static_params['nx']), 
                                       np.float32(-9999.9))
                
                # Aggregate emissions
                bands = time_info[date_key][hour_key]
                for band in bands:
                    # band reading from cache
                    arr = read_geotiff_band(
                        f"{emis_geotiff_pth}emission_{spec}_temporal.tif",
                        band['band_num'],
                        static_params
                    )
                    
                    # Vectorized aggregation
                    valid_mask = ~np.isnan(arr)
                    fill_mask = (total_emission == np.float32(-9999.9)) & valid_mask
                    add_mask = valid_mask & ~fill_mask
                    
                    total_emission[fill_mask] = arr[fill_mask]
                    total_emission[add_mask] += arr[add_mask]
                
                emission_values[ts_idx, 0, :, :, spec_idx] = total_emission

    # Clear cache
    clear_geotiff_cache()
    
    print(f"\nSUCCESS: Created {output_file}")
    print(f"- Species: {len(spec_name_str)} (nspecies starts from 1)")
    print(f"- Time steps: {len(time_steps)} (time starts from 1)")
    print(f"- Emission indices: {list(range(1, len(spec_name_str) + 1))}")
    print(f"- Grid: {static_params['nx']}x{static_params['ny']}x{static_params['nz']}")

create_chemistry_driver = create_chemistry_driver