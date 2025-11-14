##Specific date and time
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
emis_geotiff_pth = '/mnt/t/PhD_data/Downscale_Augsburg_center_10m_03022025/'  #Small_area_Emissions_Yearly/'
static_pth = '/home/vaithisa/GEO4PALM-main/JOBS/Augs_Bourges_Platz/OUTPUT/'  #/home/vaithisa/palm_model_system-v25.04/palm_model_system-v25.04/JOBS/Augsburg_passive/INPUT/Augsburg_passive_static"
static = 'Augs_Bourges_Platz_120824'
#static_pth = '/home/vaithisa/palm_model_system-v25.04/palm_model_system-v25.04/JOBS/Augsburg_small_allsector/INPUT/'
#static = 'Augsburg_small_allsector'

# Date and time range configuration
start_date = "2025-02-03 00:00:00"  # Format: "YYYY-MM-DD HH:MM:SS"
end_date = "2025-02-03 23:00:00"    # Format: "YYYY-MM-DD HH:MM:SS"

# Convert to datetime objects for easier comparison
start_dt = datetime.strptime(start_date, "%Y-%m-%d %H:%M:%S")
end_dt = datetime.strptime(end_date, "%Y-%m-%d %H:%M:%S")

# Active emission categories (edit these to select sectors)
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
    #'SumAllSectors'
]
cat_name_str = tuple(active_categories)
cat_name = np.array(cat_name_str, dtype='S64')

# Chemical species configuration
#spec_name_str = ('pm10','pm2_5')  #passive mechanism
#spec_name_str = ('pm10','no', 'no2', 'o3')   #phstatp mechanism
#spec_name_str = ('hno3','rcho','nmvoc', 'ho2','ro2','oh','no2', 'o3','no', 'h2o')  #simple mechanism 
spec_name_str = ('n2o', 'nox', 'nmvoc', 'so2', 'co', 'pm10', 'pm2_5', 'nh3', 
                'pb', 'cd', 'hg', 'as', 'ni', 'bc', 'co2', 'ch4', 'no', 'no2', 'ec', 'oc', 'na', 'so4', 
                'othmin', 'o3', 'h2o', 'oh', 'ho2', 'ro2', 'rcho', 'hno3')
spec_name = np.array(spec_name_str, dtype='S64')

# Global cache for entire resampled GeoTIFFs
_geotiff_cache = {}