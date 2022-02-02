## Script for creating chemistry driver from PALM using GRETA inventory
## Create the chemistry driver (netCDF) from the resampled emission inventory
import os
import netCDF4 as nc
import numpy as np
import pandas as pd

from chemistry_driver_config import static_pth, static, edgar_pth, chem_dr_pth, country, region, month, cat_name_str, cat_name, spec_name, spec_name_str
from chemistry_driver_main import origin_time, origin_lat, origin_lon, origin_x, origin_y, nx,ny, res, dy, dx 

print('Creating netCDF chemistry driver')
## NetCDF global attributes
ds = nc.Dataset(static_pth+static+'_chemistry', 'w', format='NETCDF4')
ds.lod = 1
ds.legacy_mode = "yes (z dimension enabled)"
ds.origin_time = origin_time
ds.origin_lat  = origin_lat
ds.origin_lon  = origin_lon
ds.origin_x    = origin_x
ds.origin_y    = origin_y
ds.resolution  = res

## Add dimensions
ds.createDimension('z', 1)
ds.createDimension('y', ny)
ds.createDimension('x', nx)
ds.createDimension('nspecies',spec_name.size)          ## Number of speices
ds.createDimension('ncat',cat_name.size)      ## Numer of emission categories
ds.createDimension('npm',3)               ## Number of PM species
ds.createDimension('nvoc',2)              ## Number of voc species
ds.createDimension('nmonthdayhour',91)    ## Time scaling factors
ds.createDimension('max_string_length',25)
ds.createDimension('nox_comp',2)          ## Composition of NOX
ds.createDimension('sox_comp',2)          ## Composition of SOX
ds.createDimension('pm_comp',3)           ## Composition of PM
ds.createDimension('voc_comp',2)          ## Composition of VOC

## Add variables and Assign values
z = ds.createVariable('z','f8',('z',))
z.units = 'm'
z.axis = 'Z'
z.long_name='distance to origin in z-direction'
z[:] = 1
## x
x = ds.createVariable('x','f8',('x',))
x.units = 'm'
x.axis = 'X'
x [:] = np.arange(10,(nx*dx+1),dx)
## y
y = ds.createVariable('y','f8',('y',))
y.units = 'm'
y.axis = 'Y'
y[:] = np.arange(10,(ny*dy+1),dy)
## emission_name
emission_name = ds.createVariable('emission_name','S1',('nspecies','max_string_length'))
emission_name.long_name = "emission species name"
emission_name.standard_name = "emission name"
emission_name.units =""
emission_name[:] = nc.stringtochar(spec_name)
## emission_index
emission_index = ds.createVariable('emission_index', 'u1',('nspecies',),fill_value=-9999)
emission_index.long_name = "emission species index"
emission_index.standard_name ="emission index"
emission_index.units =""
emission_index[:] = np.arange(1,len(emission_name)+1,1.0)
## emission_values
emission_values = ds.createVariable('emission_values', 'f4',('z','y','x','nspecies','ncat'),fill_value=-9999.9)
emission_values.long_name = "emission species values"
emission_values.standard_name ="emission values"
emission_values.lod = 1
emission_values.units = "kg/m2/year"
emission_values.coordinates = "E_UTM N_UTM lon lat"
emission_values.grid_mapping = "crsUTM: E_UTM N_UTM crsETRS: lon lat"
for nct in range(0,(ds.dimensions['ncat'].size)):
    for nsp in range(0,(ds.dimensions['nspecies'].size)):
        chem_array = np.loadtxt(chem_dr_pth+'rasters/interp/aoi/e_sum_'+spec_name_str[nsp].casefold()+'_interp_aoi.txt')
        emission_values[:,:,:,0,0] = (chem_array/1000000) ## convert to kg/m2/year

## stack height (include when palm does this)
#emission_stack_height = ds.createVariable('emission_stack_height','f4',('y','x') ,fill_value=-9999.9)
#emission_stack_height.long_name = "emission stack height"
#emission_stack_height.standard_name = "emission_stack_height"
#emission_stack_height.units = "m"
#emission_stack_height.coordinates = "E_UTM N_UTM lon lat"
#emission_stack_height.grid_mapping = "crsUTM: E_UTM N_UTM crsETRS: lon lat"
#psource_array = np.loadtxt('/cfs/home/d/u/dupreeda/MBEES/GRETA/processed/rasters/interp/aoi/point_sources.txt', delimiter=',')
#emission_stack_height[:,:] = psource_array

## emission time factor (based on EDGAR time profiles)
## months:1-12, days:13-19, work_hours:20-43, saturday_hours:44-67, sunday_hours:68-91
emission_time_factors = ds.createVariable('emission_time_factors','f4',('nmonthdayhour', 'ncat'))
emission_time_factors.long_name = "emission_time_scaling_factors"
emission_time_factors.standard_name = "emission_time_scaling_factors"
emission_time_factors.lod = 1
emission_time_factors.units = ""
## Find EDGAR time profiles for emission categories and region
tf_df = pd.read_excel(edgar_pth+'emission_time_factors.xlsx', 0)
xs = cat_name.size
emt = np.zeros((91,xs))
## Find index value of EDGAR and IPCC name
for n in range(0,xs):
    em_idx = tf_df.index[tf_df['Emission_category']==cat_name_str[n]].tolist()
    for i in em_idx:
        em_idxval = i
        edg_nm = tf_df.iloc[em_idxval]['EDGAR_sector']
        ipcc_nm = tf_df.iloc[em_idxval]['IPCC_2006']
        ## Find monthly profile (by region, IPCC, year)
        mth_df = pd.read_excel(edgar_pth+'EDGAR_temporal_profiles_r1.xlsx', 1)
        mth_idx = mth_df.index[mth_df['Region/country'] == region].tolist()
        mth_df2 = mth_df.loc[(mth_df['Region/country']==region) & (mth_df['IPCC_2006_source_category']== ipcc_nm) & (mth_df['Year']==2017)]
        if mth_df2.empty:
            mth_df2 = mth_df.loc[(mth_df['Region/country']==region) & (mth_df['IPCC_2006_source_category']== ipcc_nm) & (mth_df['Year']==0)]
        elif mth_df2.size > 17:
            mth_df2 = mth_df2.iloc[0]
        emt[0:12,n] = mth_df2[['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']].to_numpy()
        ## Find daily profile (by country, EDGAR sector, day number)
        week_df = pd.read_csv(edgar_pth+'auxiliary_tables/weekly_profiles.csv')
        for dn in range(1,8):
            week_df2 = week_df.loc[(week_df['Country_code_A3']== country) & (week_df['activity_code']== edg_nm) & (week_df['Weekday_id']== dn)]
            emt[11+dn,n] = week_df2[['daily_factor']].to_numpy()
        ## Find hourly weekday profile (by country, EDGAR sector, month, daynumber)
        hour_df = pd.read_csv(edgar_pth+'auxiliary_tables/hourly_profiles.csv')
        hour_df2 = hour_df.loc[(hour_df['Country_code_A3']== country) & (hour_df['activity_code']== edg_nm) &
                (hour_df['month_id']== month) & (hour_df['Daytype_id']==1)]
        emt[19:43,n] = hour_df2[['h1','h2','h3','h4','h5','h6','h7','h8','h9','h10','h11','h12','h13','h14','h15',
            'h16','h17','h18','h19','h20','h21','h22','h23','h24']].to_numpy()
        ## Find hourly saturday profile (by country, EDGAR sector, month)
        hour_df2 = hour_df.loc[(hour_df['Country_code_A3']== country) & (hour_df['activity_code']== edg_nm) &
                (hour_df['month_id']== month) & (hour_df['Daytype_id']==2)]
        emt[43:67,n] = hour_df2[['h1','h2','h3','h4','h5','h6','h7','h8','h9','h10','h11','h12','h13','h14','h15',
            'h16','h17','h18','h19','h20','h21','h22','h23','h24']].to_numpy()
        ## Find hourly sunday profile (by country, EDGAR sector, month)
        hour_df2 = hour_df.loc[(hour_df['Country_code_A3']== country) & (hour_df['activity_code']== edg_nm) &
                (hour_df['month_id']== month) & (hour_df['Daytype_id']==3)]
        emt[67:91,n] = hour_df2[['h1','h2','h3','h4','h5','h6','h7','h8','h9','h10','h11','h12','h13','h14','h15',
            'h16','h17','h18','h19','h20','h21','h22','h23','h24']].to_numpy()
emission_time_factors[:,:] = emt
#if (np.sum(emission_time_factors[0:12,:],axis=0)!=1 or np.sum(emission_time_factors[12:19,:],axis=0)!=1 or (np.sum(emission_time_factors[19:43,:],axis=0)!=1) or (np.sum(emission_time_factors[43:67,:],axis=0)!=1) or (np.sum(emission_time_factors[67:91,:],axis=0)!=1)):
#        print('ERROR!! Check emission time factor')
#        os.remove(static_pth+static+'_chemistry')
#        exit()

## emission category name
emission_category_name = ds.createVariable('emission_category_name','S1',('ncat','max_string_length'))
emission_category_name.long_name =  "emission category name"
emission_category_name.standard_name = "emission_cat_name"
emission_category_name.units = ""
#cat_name = np.array(['traffic exhaust'],dtype= 'S25')
emission_category_name[:] = nc.stringtochar(cat_name)
if ds.dimensions['ncat'].size != len(emission_category_name):
    print('ERROR!! Check ncat names and size!')
    os.remove(static_pth+static+'_chemistry')
    exit()
## emission category index
emission_category_index = ds.createVariable('emission_cat_index','f4',('ncat',),fill_value=-9999.9)
emission_category_index.long_name =  "emission category index"
emission_category_index.standard_name = "emission_cat_index"
emission_category_index.units = ""
emission_category_index[:] = np.arange(1,len(emission_category_name)+1,1.0)
## composition_nox (must sum to 1)
composition_nox = ds.createVariable('composition_nox','f4',('nox_comp','ncat'))
composition_nox.long_name = "composition of NOx"
composition_nox.standard_name = "composition_nox"
composition_nox.units = ""
composition_nox[0,:] = 0.8
composition_nox[1,:] = 0.2
## composition_sox (must sum to 1)
composition_sox = ds.createVariable('composition_sox','f4',('sox_comp','ncat'))
composition_sox.long_name = "composition of SOx"
composition_sox.standard_name = "composition_sox"
composition_sox.units = ""
composition_sox[0,:] = 0.95
composition_sox[1,:] = 0.05
## emission_pm_name
emission_pm_name = ds.createVariable("emission_pm_name",'S1',('npm','max_string_length'))
emission_pm_name.long_name = "PM name"
emission_pm_name.standard_name = "pm_name"
emission_pm_name.units = ""
pm_name = np.array(['PM10',' PM2.5','PM1'],dtype= 'S25')
emission_pm_name[:] = nc.stringtochar(pm_name)
if ds.dimensions['npm'].size != len(emission_pm_name):
    print('ERROR!! Check npm names and size!')
    os.remove(static_pth+static+'_chemistry')
    exit()
## composition_pm (must sum to 1)
composition_pm = ds.createVariable("composition_pm",'f4',('pm_comp','npm'))
composition_pm.long_name = "composition of PM"
composition_pm.standard_name = "composition_PM"
composition_pm.units = ""
composition_pm[0,:] = 0.9
composition_pm[1,:] = 0.09
composition_pm[2,:] = 0.01
## emission_voc_name
emission_voc_name = ds.createVariable("emission_voc_name",'S1',('nvoc','max_string_length'))
emission_voc_name.long_name = "VOC name"
emission_voc_name.standard_name = "voc_name"
emission_voc_name.units = ""
voc_name = np.array(['NMVOC','MVOC'],dtype= 'S25')
emission_voc_name[:] = nc.stringtochar(voc_name)
if ds.dimensions['nvoc'].size != len(emission_voc_name):
    print('ERROR!! Check nvoc names and size!')
    os.remove(static_pth+static+'_chemistry')
    exit()
## composition_voc (must sum to 1)
composition_voc = ds.createVariable("composition_voc",'f4',('voc_comp','nvoc'))
composition_voc.long_name = "composition of VOC"
composition_voc.standard_name = "composition_VOC"
composition_voc.units =""
composition_voc[0,:] = 0.9
composition_voc[1,:] = 0.1

print('Chemistry driver created!')
ds.close()
