#species ordered
import os
import re
import numpy as np
import netCDF4 as nc
from osgeo import gdal
from chemistry_driver_config import *
from chemistry_driver_main import read_geotiff_band, clear_geotiff_cache
import datetime 

# Pre-compiled regex for faster parsing - UPDATED for h00-h23 format
BAND_PATTERN = re.compile(r'^([A-Za-z_]+)_h(\d{2})_(\d{8})$')

def parse_band_description(desc):
    """Fast parsing with pre-compiled regex"""
    match = BAND_PATTERN.match(desc)
    return match.groups() if match else None

def get_time_bands_info_fast(geotiff_path, filter_traffic=False):
    """Optimized band info extraction with optional traffic filtering"""
    ds = gdal.Open(geotiff_path, gdal.GA_ReadOnly)
    if not ds: 
        raise ValueError(f"Can't open {geotiff_path}")
    
    time_info = {}
    
    if filter_traffic:
        active_set = set(traffic_sectors)
    else:
        active_set = set(active_categories)
    
    for band_num in range(1, ds.RasterCount + 1):
        band = ds.GetRasterBand(band_num)
        desc = band.GetDescription()
        if desc:
            parsed = parse_band_description(desc)
            if parsed and parsed[0] in active_set:
                band_date_str = parsed[2]
                band_hour_str = parsed[1]
                band_hour = int(band_hour_str)
                
                band_dt = datetime.datetime.strptime(band_date_str, "%Y%m%d") + \
                         datetime.timedelta(hours=band_hour)
                
                if start_dt <= band_dt <= end_dt:
                    if band_date_str not in time_info:
                        time_info[band_date_str] = {}
                    if f"h{band_hour_str}" not in time_info[band_date_str]:
                        time_info[band_date_str][f"h{band_hour_str}"] = []
                    
                    time_info[band_date_str][f"h{band_hour_str}"].append({
                        'band_num': band_num,
                        'category': parsed[0],
                        'hour_num': band_hour,
                        'hour_str': band_hour_str,
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
        hour_num = current_dt.hour
        hour_str = f"h{hour_num:02d}"
        
        time_steps.append({
            'date': date_str,
            'hour': hour_str,
            'hour_num': hour_num,
            'datetime': current_dt
        })
        
        current_dt += datetime.timedelta(hours=1)
    
    return time_steps

def filter_negative_and_zero_values(emission_array):
    """Set all values < 0 to 0 in the emission array."""
    filtered_array = emission_array.copy()
    filtered_array[filtered_array < 0] = np.float32(0)
    return filtered_array

def create_chemistry_driver(static_params):
    """Chemistry driver creation with traffic species from spec_name_str"""
    print("\nCreating PALM LOD2 Chemistry Driver")
    print(f"Date range: {start_date} to {end_date}")

    # FORCE z dimension to 1 for surface-only emissions
    static_params['nz'] = 1
    print("Note: z dimension forced to 1 (surface emissions only)")
    
    if tag == "traffic":
        print("Traffic species separation: ENABLED")
        print(f"Traffic sectors: {', '.join(traffic_sectors)}")
        print("Zero/negative value filtering: ENABLED (all values < 0 will be set to 0)")
    else:
        print("Traffic species separation: DISABLED")
        print("Zero/negative value filtering: ENABLED (all values < 0 will be set to 0)")
    
    species_mapping = {
        'n2o': 'N2O', 'nox': 'NOX', 'nmvoc': 'RH', 'so2': 'H2SO4', 'co': 'CO',
        'pm10': 'PM10', 'pm2_5': 'PM25', 'nh3': 'NH3', 'pb': 'PB', 'cd': 'CD',
        'hg': 'HG', 'as': 'AS', 'ni': 'NI', 'bc': 'BC', 'co2': 'CO2', 'ch4': 'CH4',
        'no': 'NO', 'no2': 'NO2', 'ec': 'EC', 'na': 'NA', 'so4': 'SO4', 'ocnv': 'OCNV',
        'othmin': 'OTHMIN', 'o3': 'O3', 'hno3': 'HNO3', 'rcho': 'RCHO', 'ho2': 'HO2',
        'ro2': 'RO2', 'oh': 'OH', 'h2o': 'H2O', 'ocsv': 'OCSV'
    }
    
    all_species_to_process = list(spec_name_str)
    
    # Create uppercase versions for output
    uppercase_spec_names = []
    for spec in spec_name_str:
        is_traffic = spec.lower().endswith('_tra')
        
        # FIX: Use replace instead of slicing to correctly remove _tra suffix
        if is_traffic:
            base_spec = spec.lower().replace('_tra', '')  # "pm10_tra" -> "pm10"
        else:
            base_spec = spec.lower()
        
        if base_spec in species_mapping:
            base_upper = species_mapping[base_spec]
        else:
            base_upper = base_spec.upper().replace('_', '')
        
        if is_traffic:
            uppercase_spec_names.append(f"{base_upper}_tra")
        else:
            uppercase_spec_names.append(base_upper)
    
    print(f"Total species to process: {len(all_species_to_process)}")
    print(f"Species order (preserved from config):")
    for i, spec in enumerate(all_species_to_process):
        print(f"  {i+1}: {spec} -> {uppercase_spec_names[i]}")
    
    # Collect temporal emission data for ALL species
    print("Scanning emission files...")
    all_time_info = {}
    
    for spec in all_species_to_process:
        # FIX: Check for _tra suffix (not _traffic)
        is_traffic_species = spec.lower().endswith('_tra')
        
        # FIX: Use replace to correctly get base species name
        if is_traffic_species:
            base_spec = spec.lower().replace('_tra', '')  # "pm10_tra" -> "pm10", "pm2_5_tra" -> "pm2_5", "no2_tra" -> "no2"
        else:
            base_spec = spec
        
        gt_path = f"{emis_geotiff_pth}emission_{base_spec}_temporal.tif"
        if os.path.exists(gt_path):
            print(f"  Found: {spec} ({'traffic sectors only' if is_traffic_species else 'all sectors'}) -> file: emission_{base_spec}_temporal.tif")
            all_time_info[spec] = get_time_bands_info_fast(gt_path, filter_traffic=is_traffic_species)
        else:
            print(f"  Missing: {spec} (file: {os.path.basename(gt_path)})")
            all_time_info[spec] = {}
    
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
        attrs = {
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
            'traffic_species_enabled': 'yes' if tag == "traffic" else 'no',
            'zero_negative_filter': 'yes (all values < 0 set to 0)',
            'Conventions': 'CF-1.7'
        }
        
        if tag == "traffic":
            attrs['traffic_sectors'] = ', '.join(traffic_sectors)
        
        ds.setncatts(attrs)
        
        n_time = len(time_steps)
        n_species = len(all_species_to_process)
        
        ds.createDimension('z', static_params['nz'])
        ds.createDimension('field_length', 64)
        ds.createDimension('x', static_params['nx'])
        ds.createDimension('y', static_params['ny'])
        ds.createDimension('nspecies', n_species)
        ds.createDimension('time', None)
        
        z = ds.createVariable('z', 'f4', ('z',))
        z[:] = np.float32(static_params['z_origin'])  # Single value at surface height
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
        
        nspecies_var = ds.createVariable('nspecies', 'i4', ('nspecies',))
        nspecies_var[:] = np.arange(1, n_species + 1, dtype=np.int32)
        nspecies_var.long_name = "nspecies"
        
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
        
        print("Processing emissions data...")
        print("Note: All values < 0 will be filtered and set to 0")
        
        for spec_idx, spec_name in enumerate(all_species_to_process):
            # FIX: Check for _tra suffix
            is_traffic_species = spec_name.lower().endswith('_tra')
            
            # FIX: Use replace to correctly get base species name
            if is_traffic_species:
                base_spec = spec_name.lower().replace('_tra', '')  # "pm10_tra" -> "pm10", "pm2_5_tra" -> "pm2_5"
            else:
                base_spec = spec_name
            
            print(f"  Processing {spec_idx+1}/{len(all_species_to_process)}: {spec_name} "
                  f"({'traffic sectors only' if is_traffic_species else 'all sectors'}) "
                  f"-> reading from emission_{base_spec}_temporal.tif")
            
            species_emissions = np.full((n_time, static_params['ny'], static_params['nx']), 
                                       np.float32(-9999.9))
            
            if spec_name not in all_time_info or not all_time_info[spec_name]:
                print(f"    WARNING: No time info for {spec_name}, filling with -9999.9")
                emission_values[:, 0, :, :, spec_idx] = species_emissions
                continue
            
            time_info = all_time_info[spec_name]
            
            for ts_idx, ts in enumerate(time_steps):
                date_key = ts['date']
                hour_key = ts['hour']
                
                if date_key not in time_info or hour_key not in time_info[date_key]:
                    continue
                
                total_emission = np.full((static_params['ny'], static_params['nx']), 
                                       np.float32(-9999.9))
                
                bands = time_info[date_key][hour_key]
                for band in bands:
                    # FIX: Read from BASE species file
                    arr = read_geotiff_band(
                        f"{emis_geotiff_pth}emission_{base_spec}_temporal.tif",
                        band['band_num'],
                        static_params
                    )
                    
                    valid_mask = ~np.isnan(arr)
                    fill_mask = (total_emission == np.float32(-9999.9)) & valid_mask
                    add_mask = valid_mask & ~fill_mask
                    
                    arr_g_per_sec = arr * (1000.0 / 3600.0)
                    
                    total_emission[fill_mask] = arr_g_per_sec[fill_mask]
                    total_emission[add_mask] += arr_g_per_sec[add_mask]
                
                filtered_emission = filter_negative_and_zero_values(total_emission)
                species_emissions[ts_idx] = filtered_emission
            
            emission_values[:, 0, :, :, spec_idx] = species_emissions
    
    clear_geotiff_cache()
    
    print(f"\nSUCCESS: Created {output_file}")
    print(f"- Total species: {len(all_species_to_process)}")
    print(f"- Species list (in order):")
    for i, (orig, upper) in enumerate(zip(all_species_to_process, uppercase_spec_names)):
        traffic_indicator = " (traffic)" if orig.lower().endswith('_tra') else ""
        print(f"    {i+1}: {orig} -> {upper}{traffic_indicator}")
    print(f"- Time steps: {len(time_steps)}")
    print(f"- Date range: {start_date} to {end_date}")
    print(f"- Grid: {static_params['nx']}x{static_params['ny']}x{static_params['nz']}")
    print(f"- Units: g/m²/s")
    print(f"- Traffic species separation: {'ENABLED' if tag == 'traffic' else 'DISABLED'}")
    
    return output_file