# chemistry_driver_palm

These scripts are used to create the chemistry driver for PALM from the GRETA emission inventory.

The process to downscale the GREAT emission inventory for a PALM simualtion follows the following process:

- The configuration are read from the `chemistry_driver_config. py` after being called by `chemistry_driver_main.py`. The configurations include the direcories to the processed chemistry, static layer as well as the emission category and emission name.
Note: Only use the emission categories in the `emission_time_factors.xlsx` as these are used to obtain the emission time factors for the netcdf driver.

- Next the static driver for the simualtion is read, necessary data is extracted.


- Next, the emission inventory is resampled and resized. NOTE! Improve on the existig scripts are in process.


- Lastly, the resampled emission inventory is written to a netcdf file (`chemistry driver netcdf creator.py`). Emission_time_factors are read from the EDGAR profiles using the emission category, month and region/country. EDGAR time profiles can found here: https://edgar.jrc.ec.europa.eu/dataset_temp_profile

TODO:
- Improve resampling and resizing of emission inventory.
- Update README file.
- Update document describing chemistry driver dimensions and variables.
