#species order
import numpy as np
from pyproj import Proj, Transformer
import warnings
import os
import yaml
import sys

from datetime import datetime

# Suppress warnings
warnings.filterwarnings("ignore")
warnings.simplefilter("ignore")

# =============================================================================
# Load configuration from YAML file passed as argument
# =============================================================================
# Usage: python chemistry_driver_main.py chem_config.yaml

if len(sys.argv) < 2:
    print("Usage: python chemistry_driver_main.py <config.yaml>")
    print("Example config files in : chem_driver_config.yaml")
    sys.exit(1)

config_file = sys.argv[1]
with open(config_file, "r") as f:
    _cfg = yaml.safe_load(f)

print(f'Reading PALM chemistry configuration from: {config_file}')

# Projection configurations
config_proj = "EPSG:25832"  # UTM Zone 32N
default_proj = "EPSG:4326"  # WGS84

# Coordinate transformers
transformer_to_utm = Transformer.from_crs(default_proj, config_proj, always_xy=True)
transformer_to_wgs = Transformer.from_crs(config_proj, default_proj, always_xy=True)

# Path configurations
emis_geotiff_pth = _cfg["paths"]["emis_geotiff_pth"]
static_pth       = _cfg["paths"]["static_pth"]
static           = _cfg["paths"]["static"]

# Date and time range configuration
start_date = _cfg["time"]["start_date"]
end_date   = _cfg["time"]["end_date"]

# Convert to datetime objects for easier comparison
start_dt = datetime.strptime(start_date, "%Y-%m-%d %H:%M:%S")
end_dt = datetime.strptime(end_date, "%Y-%m-%d %H:%M:%S")

# Traffic sectors (sectors treated as "traffic" for _tra species filtering)
traffic_sectors = _cfg["traffic_sectors"]

# Tag is set automatically based on whether any species ends with '_tra'
tag = "traffic" if any(s.endswith('_tra') for s in _cfg["species"]) else ""

# Active emission categories
active_categories = list(_cfg["active_categories"])
cat_name_str = tuple(active_categories)
cat_name = np.array(cat_name_str, dtype='S64')

# Species list
spec_name_str = tuple(_cfg["species"])

# Global cache for entire resampled GeoTIFFs
_geotiff_cache = {}