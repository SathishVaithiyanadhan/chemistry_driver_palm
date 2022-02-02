# Script for creating chemistry driver from PALM using GRETA inventory
## Files with configuration for chemistry driver - make changes here for each simulation
import numpy as np

print('Reading configuration')
## Directories and filenames
edgar_pth   = '/cfs/home/d/u/dupreeda/MBEES/GRETA/EDGAR/'
chem_dr_pth = '/cfs/home/d/u/dupreeda/MBEES/GRETA/chemistry_driver/processed/'
static_pth  = '/cfs/home/d/u/dupreeda/MBEES/PALM/palm_model_system-v21.10/JOBS/augs10/INPUT/'
static      = 'augs10'
country     = 'DEU'
region      = 22
month       = 6

## Define emission driver details (Emission categories and species)
cat_name_str = ('Transport',)
cat_name = np.array([cat_name_str],dtype= 'S25')       ## emission_category_name
spec_name_str = ('NOX',)                               ## emission_name string
spec_name = np.array([spec_name_str],dtype= 'S25')
