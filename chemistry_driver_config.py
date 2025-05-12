# Script for creating chemistry driver from PALM using emission GeoTIFFs
import numpy as np

print('Reading configuration')
## Directories and filenames
emis_geotiff_pth = '/home/vaithisa/Downscale_Emissions/Small_area_Emissions_Yearly/'
static_pth = '/home/vaithisa/GEO4PALM-main/palm_model_system-v22.10/JOBS/Augsburg_konig/INPUT/'
static = 'Augsburg_konig'
country = 'DEU'
month = 2  # February

## Define emission categories matching your GeoTIFF bands
cat_name_str = (
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
    'SumAllSectors'
)
cat_name = np.array(cat_name_str, dtype='S25')

## Define species (must match your GeoTIFF filenames)
spec_name_str = ('n2o', 'nox', 'nmvoc', 'so2', 'co', 'pm10', 'pm2_5', 'nh3', 'pb', 'cd', 'hg', 'as', 'ni', 'bc', 'co2', 'ch4')  # 'n2o', 'nox', 'nmvoc', 'so2', 'co', 'pm10', 'pm2_5', 'nh3', 'pb', 'cd', 'hg', 'as', 'ni', 'bc', 'co2', 'ch4')  # Lowercase to match filenames
spec_name = np.array(spec_name_str, dtype='S25')

# Simple time profile defaults (monthly, weekly, hourly factors)
default_monthly = np.array([1/12]*12)  # Equal across months
default_daily = np.array([1/7]*7)      # Equal across days
default_hourly = np.array([1/24]*24)   # Equal across hours