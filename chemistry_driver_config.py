#species order
import numpy as np
from pyproj import Proj, Transformer
import warnings
import os
from datetime import datetime

# Suppress warnings
warnings.filterwarnings("ignore")
warnings.simplefilter("ignore")

print('Reading PALM chemistry configuration')

# Projection configurations
config_proj = "EPSG:25832"  # UTM Zone 32N
default_proj = "EPSG:4326"  # WGS84

# Coordinate transformers
transformer_to_utm = Transformer.from_crs(default_proj, config_proj, always_xy=True)
transformer_to_wgs = Transformer.from_crs(config_proj, default_proj, always_xy=True)

# Path configurations
emis_geotiff_pth = '/home/vaithisa/Downscale_Emissions_simple/downscale/'
static_pth = '/home/vaithisa/GEO4PALM-main/JOBS/Augs_Bourges_Platz/OUTPUT/'
static = 'Augs_Bourges_Platz_400'

# Date and time range configuration
start_date = "2024-08-11 00:00:00"  # Format: "YYYY-MM-DD HH:MM:SS"
end_date = "2024-08-11 23:00:00"    # Format: "YYYY-MM-DD HH:MM:SS"

# Convert to datetime objects for easier comparison
start_dt = datetime.strptime(start_date, "%Y-%m-%d %H:%M:%S")
end_dt = datetime.strptime(end_date, "%Y-%m-%d %H:%M:%S")

# Traffic tag configuration
# Set tag = "traffic" to enable traffic-specific species separation
# Set tag = "" or any other value to disable
tag = "traffic"  # Enable traffic species separation

# Traffic sectors (these will be separated when tag = "traffic")
traffic_sectors = ['F_RoadTransport']  # 'I_OffRoad' can be added if needed

# Species that should have traffic versions (when tag = "traffic")
# These are the base species names (without _traffic suffix)
tag_spec_name_str = ('no', 'no2')  # Create traffic versions for these species

# Active emission categories
active_categories = [
    'A_PublicPower', 
    'B_Industry', 
    'C_OtherStationaryComb', 
    'D_Fugitives',
    'E_Solvents', 
    'F_RoadTransport', 
    'G_Shipping', 
    'H_Aviation',
    'I_OffRoad', 
    'J_Waste', 
    'K_AgriLivestock', 
    'L_AgriOther',
]
cat_name_str = tuple(active_categories)
cat_name = np.array(cat_name_str, dtype='S64')


# This can contain both regular species and traffic species (with _traffic suffix)
#spec_name_str = ('pm10','pm2_5')  # passive mechanism
# spec_name_str = ('pm10','no', 'no2', 'o3')   # phstatp mechanism
#spec_name_str = ('hno3', 'rcho', 'nmvoc', 'ho2', 'no2', 'ro2', 'no2_traffic', 'no_traffic', 'oh', 'o3', 'no', 'h2o')  #simple_traffic
#spec_name_str = ('so2', 'nh3', 'oc', 'hno3','rcho','nmvoc', 'no2', 'no2_traffic', 'ho2','o3', 'no_traffic', 'no', 'oh', 'ro2', 'h2o') #salsa+simple_tra
spec_name_str = ('hno3', 'rcho', 'nmvoc', 'ho2', 'ro2', 'oh', 'no2',  'o3', 'no', 'h2o')  #simple

# Global cache for entire resampled GeoTIFFs
_geotiff_cache = {}