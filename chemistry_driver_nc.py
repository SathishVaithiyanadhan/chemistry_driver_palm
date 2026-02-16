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
    
    # Use regular dict
    time_info = {}
    
    # Determine which categories to include
    if filter_traffic:
        # When filter_traffic=True, we only want traffic sectors
        active_set = set(traffic_sectors)
    else:
        # When filter_traffic=False, we want all sectors
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
    """
    Set all values ≤ 0 to NaN in the emission array.
    This ensures only positive emission values are kept.
    """
    # Create a copy to avoid modifying the original
    filtered_array = emission_array.copy()
    
    # Set values ≤ 0 to NaN
    filtered_array[filtered_array <= 0] = np.float32(-9999.9)
    
    return filtered_array

def create_chemistry_driver(static_params):
    """Chemistry driver creation with selective traffic species support and zero/negative value filtering"""
    print("\nCreating PALM LOD2 Chemistry Driver")
    print(f"Date range: {start_date} to {end_date}")
    
    # Check if traffic tag is enabled
    if tag == "traffic":
        print("Traffic species separation: ENABLED")
        print(f"Traffic sectors: {', '.join(traffic_sectors)}")
        print(f"Species with traffic versions: {', '.join(tag_spec_name_str)}")
        print("Zero/negative value filtering: ENABLED (all values ≤ 0 will be set to NaN)")
    else:
        print("Traffic species separation: DISABLED")
        print("Zero/negative value filtering: ENABLED (all values ≤ 0 will be set to NaN)")
    
    # Species mapping
    species_mapping = {
        'n2o': 'N2O', 'nox': 'NOX', 'nmvoc': 'RH', 'so2': 'H2SO4', 'co': 'CO',
        'pm10': 'PM10', 'pm2_5': 'PM25', 'nh3': 'NH3', 'pb': 'PB', 'cd': 'CD',
        'hg': 'HG', 'as': 'AS', 'ni': 'NI', 'bc': 'BC', 'co2': 'CO2', 'ch4': 'CH4',
        'no': 'NO', 'no2': 'NO2', 'ec': 'EC', 'na': 'NA', 'so4': 'SO4', 'oc': 'OCNV',
        'othmin': 'OTHMIN', 'o3': 'O3', 'hno3': 'HNO3', 'rcho': 'RCHO', 'ho2': 'HO2',
        'ro2': 'RO2', 'oh': 'OH', 'h2o': 'H2O'
    }
    
    # Prepare species names based on traffic tag
    all_species_to_process = []
    uppercase_spec_names = []
    
    # First, add all regular species (from spec_name_str)
    for spec in spec_name_str:
        base_name = species_mapping.get(spec, spec.upper().replace('_', ''))
        all_species_to_process.append(spec)
        uppercase_spec_names.append(base_name)
    
    # Then, if traffic tag is enabled, add traffic species ONLY for those in tag_spec_name_str
    if tag == "traffic":
        for spec in tag_spec_name_str:
            if spec in spec_name_str:  # Only create traffic versions for species that exist in regular list
                base_name = species_mapping.get(spec, spec.upper().replace('_', ''))
                traffic_spec_name = f"{spec}_traffic"
                all_species_to_process.append(traffic_spec_name)
                uppercase_spec_names.append(f"{base_name}_traffic")
            else:
                print(f"  Warning: {spec} is in tag_spec_name_str but not in spec_name_str. Skipping traffic version.")
    
    print(f"Total species to process: {len(all_species_to_process)}")
    print(f"Regular species: {len(spec_name_str)}")
    if tag == "traffic":
        print(f"Traffic species: {len(tag_spec_name_str)}")
    
    # Collect temporal emission data for ALL species (regular and traffic)
    print("Scanning emission files...")
    all_time_info = {}
    
    for spec in all_species_to_process:
        # Determine if this is a traffic species
        is_traffic_species = spec.endswith('_traffic')
        
        # Get base species name (remove _traffic suffix if present)
        base_spec = spec.replace('_traffic', '') if is_traffic_species else spec
        
        gt_path = f"{emis_geotiff_pth}emission_{base_spec}_temporal.tif"
        if os.path.exists(gt_path):
            print(f"  Found: {spec} ({'traffic sectors only' if is_traffic_species else 'all sectors'})")
            # For traffic species, filter to only traffic sectors
            all_time_info[spec] = get_time_bands_info_fast(gt_path, filter_traffic=is_traffic_species)
        else:
            print(f"  Missing: {spec} (file: {os.path.basename(gt_path)})")
            all_time_info[spec] = {}  # Empty dict for missing files
    
    # Create ALL expected time steps
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
            'zero_negative_filter': 'yes (all values ≤ 0 set to NaN)',
            'Conventions': 'CF-1.7'
        }
        
        # Add traffic-specific attributes only if traffic is enabled
        if tag == "traffic":
            attrs['traffic_sectors'] = ', '.join(traffic_sectors)
            attrs['species_with_traffic_versions'] = ', '.join(tag_spec_name_str)
        
        ds.setncatts(attrs)
        
        # Dimensions
        n_time = len(time_steps)
        n_species = len(all_species_to_process)
        
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
        
        # Process emissions for all species
        print("Processing emissions data...")
        print("Note: All values ≤ 0 will be filtered and set to NaN")
        
        for spec_idx, spec_name in enumerate(all_species_to_process):
            is_traffic_species = spec_name.endswith('_traffic')
            base_spec = spec_name.replace('_traffic', '') if is_traffic_species else spec_name
            
            print(f"  Processing {spec_idx+1}/{len(all_species_to_process)}: {spec_name} "
                  f"({'traffic sectors only' if is_traffic_species else 'all sectors'})")
            
            # Initialize emissions array for this species
            species_emissions = np.full((n_time, static_params['ny'], static_params['nx']), 
                                       np.float32(-9999.9))
            
            if spec_name not in all_time_info or not all_time_info[spec_name]:
                # No data for this species
                emission_values[:, 0, :, :, spec_idx] = species_emissions
                continue
            
            time_info = all_time_info[spec_name]
            
            # Process each time step
            for ts_idx, ts in enumerate(time_steps):
                date_key = ts['date']
                hour_key = ts['hour']
                
                if date_key not in time_info or hour_key not in time_info[date_key]:
                    # No data for this time step
                    continue
                
                # Initialize with fill values
                total_emission = np.full((static_params['ny'], static_params['nx']), 
                                       np.float32(-9999.9))
                
                # Aggregate emissions for this time step
                bands = time_info[date_key][hour_key]
                for band in bands:
                    # Read band data
                    arr = read_geotiff_band(
                        f"{emis_geotiff_pth}emission_{base_spec}_temporal.tif",
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
                
                # FILTER: Set all values ≤ 0 to NaN (fill value)
                total_emission = filter_negative_and_zero_values(total_emission)
                
                species_emissions[ts_idx] = total_emission
            
            # Store in NetCDF variable
            emission_values[:, 0, :, :, spec_idx] = species_emissions
    
    # Clear cache
    clear_geotiff_cache()
    
    print(f"\nSUCCESS: Created {output_file}")
    print(f"- Total species: {len(all_species_to_process)}")
    print(f"- Regular species: {len(spec_name_str)}")
    if tag == "traffic":
        print(f"- Traffic species: {len(tag_spec_name_str)}")
        print(f"- Species with traffic versions: {', '.join(tag_spec_name_str)}")
    print(f"- Time steps: {len(time_steps)}")
    print(f"- Date range: {start_date} to {end_date}")
    print(f"- Grid: {static_params['nx']}x{static_params['ny']}x{static_params['nz']}")
    print(f"- Units: g/m²/s")
    print(f"- Traffic species separation: {'ENABLED' if tag == 'traffic' else 'DISABLED'}")
    
    return output_file

# Replace the original function
create_chemistry_driver = create_chemistry_driver