import numpy as np
from pyproj import Proj, Transformer
import warnings
import os

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
emis_geotiff_pth = '/home/vaithisa/Downscale_Emissions/Downscale_10m/'  #Small_area_Emissions_Yearly/'
static_pth = '/home/vaithisa/GEO4PALM-main/palm_model_system-v22.10/JOBS/Augsburg_konig/INPUT/'
static = 'Augsburg_konig'

# Active emission categories (edit these to select sectors)
active_categories = [
      #'A_PublicPower', 
    #'B_Industry', 
    #'C_OtherStationaryComb', 
    #'D_Fugitives',
    #'E_Solvents', 
    'F_RoadTransport', 
    #'G_Shipping', 
    #'H_Aviation',
    #'I_OffRoad', 
    #'J_Waste', 
    #'K_AgriLivestock', 
    #'L_AgriOther',
    #'SumAllSectors'
]
cat_name_str = tuple(active_categories)
cat_name = np.array(cat_name_str, dtype='S64')

# Chemical species configuration
spec_name_str = ('pm10',)
#spec_name_str = ('n2o', 'nox', 'nmvoc', 'so2', 'co', 'pm10', 'pm2_5', 'nh3', 
                #'pb', 'cd', 'hg', 'as', 'ni', 'bc', 'co2', 'ch4')
spec_name = np.array(spec_name_str, dtype='S64')