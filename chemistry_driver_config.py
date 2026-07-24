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
emis_geotiff_pth = '/hpc/gpfs2/home/u/vaithisa/UniA/Downscale_Emissions_simple/Downscale_10m_3days/' #'/mnt/t/PhD_data/downscale/UA_CLC_100m_3days/'  #'/home/vaithisa/Downscale_Emissions_simple/downscale/' #'/mnt/d/downscaled_emissions/UA_CLC_100m_3days/'
#static_pth = '/home/vaithisa/GEO4PALM-main/JOBS/Augs_Bourges_Platz/OUTPUT/'   
static_pth = '/hpc/gpfs2/scratch/u/vaithisa/palm_25.10/palm_mbees/palm/JOBS/Salsa_tra128/INPUT/'  #'/mnt/d/downscaled_emissions/UA_CLC_100m_3days/'
static = 'Salsa_tra128_static'

# Date and time range configuration
start_date = "2024-08-25 00:00:00"  # Format: "YYYY-MM-DD HH:MM:SS"
end_date = "2024-08-25 23:00:00"    # Format: "YYYY-MM-DD HH:MM:SS"

# Convert to datetime objects for easier comparison
start_dt = datetime.strptime(start_date, "%Y-%m-%d %H:%M:%S")
end_dt = datetime.strptime(end_date, "%Y-%m-%d %H:%M:%S")

# Traffic tag configuration
# Set tag = "traffic" to enable traffic-specific species separation
# Set tag = "" or any other value to disable
tag = ""  # Enable traffic species separation

# Traffic sectors (these will be separated when tag = "traffic")
traffic_sectors = ['F_RoadTransport']  # 'I_OffRoad' can be added if needed

# Species that should have traffic versions (when tag = "traffic")
# These are the base species names (without _traffic suffix)
#tag_spec_name_str = ('no', 'no2', 'pm10','pm2_5')  # Create traffic versions for these species
tag_spec_name_str = () 

# Active emission categories
active_categories = [
    'A_PublicPower', 
    'B_Industry', 
    'C_OtherStationaryComb', 
    'D_Fugitives',
    'E_Solvents', 
    #'F_RoadTransport', 
    'G_Shipping', 
    'H_Aviation',
    'I_OffRoad', 
    'J_Waste', 
    'K_AgriLivestock', 
    'L_AgriOther',
]
cat_name_str = tuple(active_categories)
cat_name = np.array(cat_name_str, dtype='S64')


# This can contain both regular species and traffic species (with _tra suffix)
# Full list matching salsa+simple_tra mechanism (32 species in mechanism order)
# spec_name_str = ('pm2_5',)  # passive mechanism
# spec_name_str = ('pm10','no', 'no2', 'o3')   # phstatp mechanism
# spec_name_str = ('hno3', 'rcho', 'nmvoc', 'ho2', 'no2', 'ro2', 'no2_traffic', 'no_traffic', 'oh', 'o3', 'no', 'h2o')  #simple_traffic
#spec_name_str = ('so2', 'nh3', 'oc', 'hno3','rcho','nmvoc', 'no2', 'no2_traffic', 'ho2','o3', 'no_traffic', 'no', 'oh', 'ro2', 'h2o') #salsa+simple_tra
#spec_name_str = ('hno3', 'rcho', 'nmvoc', 'ho2', 'ro2', 'oh', 'no2',  'o3', 'no', 'h2o')  #simple
#spec_name_str = ('so4', 'nh3', 'ocnv', 'ocsv', 'pm10','pm2_5', 'pm10_tra', 'pm2_5_tra', 'hno3','rcho','nmvoc', 'no2_tra', 'ho2','ro2', 'no_tra', 'oh', 'o3', 'no', 'no2', 'h2o') # simple+salsa_tra mechanism
#spec_name_str = ('so4', 'nh3', 'ocnv', 'ocsv', 'hno3','rcho','nmvoc','o3', 'oh', 'no2', 'no', 'ho2','ro2', 'h2o') # simple+salsa mechanism
spec_name_str = ('so4', 'nh3', 'ocnv', 'ocsv', 'so4_tra', 'nh3_tra', 'ocnv_tra', 'ocsv_tra', 'pm10', 'pm2_5', 'pm10_tra', 'pm2_5_tra', 'hno3', 'rcho', 
                 'h2o_tra', 'rcho_tra', 'hno3_tra', 'nmvoc', 'nmvoc_tra', 'ro2_tra', 'ro2', 'o3', 'no2', 'no', 
                 'oh', 'ho2', 'ho2_tra', 'oh_tra', 'no2_tra', 'no_tra', 'o3_tra', 'h2o')  # simple+salsa_tra mechanism (32 species, full mechanism order)

# Global cache for entire resampled GeoTIFFs
_geotiff_cache = {}