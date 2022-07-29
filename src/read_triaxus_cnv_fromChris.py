import numpy as np
#from seabird.cnv import fCNV
import pandas
#from matplotlib import pyplot as plt
import os
import gsw as TEOS_10
import scipy.signal as scipy_signal
#import scipy.stats as scipy_stats
import xarray

#import bz2
#import gzip
#import linecache
#import re
#import warnings
#import zipfile
#from datetime import datetime
#from io import StringIO
#from pathlib import Path
import cartopy.feature as cfeature


import read_seabird_cnv_fromChris  as seabird_reader 
def compute_CT_SA_z_on_profile(in_situ_temp_profile, prac_salinity_profile,pressure_profile,longitude_profile,latitude_profile):
    
    absolute_salinity_profile  = TEOS_10.SA_from_SP(prac_salinity_profile,pressure_profile,longitude_profile,latitude_profile)
    conservative_temp_profile  = TEOS_10.CT_from_t(absolute_salinity_profile,in_situ_temp_profile,pressure_profile)
    depth_profile              = TEOS_10.z_from_p(pressure_profile,latitude_profile)
    
    return absolute_salinity_profile, conservative_temp_profile, depth_profile


#================================ HEADER - User should set paths and file names =================================================================================#
CNV_FILE_PATH = './../../triaxus/processing/in2022_v06/seasave/in2022_v06_01'
#CNV_FILE_PATH = '/run/user/336728/gvfs/smb-share:server=data.investigator.csiro.au,share=voyages/in2019_v05/triaxus/processing/in2019_v05/seasave/in2019_v05_02'
#CNV_FILE_NAME = 'in2019_v05_02_002-with_ECOtriplet.cnv'
CNV_FILE_NAME = 'in2022_v06_01_002.cnv'

OUTPUT_DATA_PATH = '../DATA/TRIAXUS/Processing/Processed_Data/TOW_1'

#================================================================================================================================================================#

#OUTPUT_FILE_NAME = file_name.split('.cnv')[0] + '_processed_ascending_descending.nc'
OUTPUT_FILE_NAME = CNV_FILE_NAME.split('.cnv')[0] + '_processed_ascending_descending.nc'

triaxus_dataframe = seabird_reader.from_cnv(os.path.join(CNV_FILE_PATH,CNV_FILE_NAME))

triaxus_time         = triaxus_dataframe['timeS'].values
start_time           = pandas.Timestamp(triaxus_dataframe._metadata['time'])

triaxus_time         = start_time + pandas.to_timedelta(triaxus_time, unit='s')

triaxus_latitude     = triaxus_dataframe['latitude'].values
triaxus_longitude    = triaxus_dataframe['longitude'].values
triaxus_pressure     = triaxus_dataframe.index.values


triaxus_temperature_1  = triaxus_dataframe['t090C'].values
triaxus_salinity_1     = triaxus_dataframe['sal00'].values



triaxus_temperature_2  = triaxus_dataframe['t190C'].values
triaxus_salinity_2     = triaxus_dataframe['sal11'].values
#other_variables_to_get = ['sbeox0Mm/L','sbeox1Mm/L','par','CStarTr0','wetCDOM','flECO-AFL','turbWETbb0']
other_variables_to_get = []
#new_var_names          = ['O2_1','O2_2','par','Tras','wetCDOM','flECO-AFL','turbWETbb0']
new_var_names          = ['O2_1','O2_2']



top_peaks,peak_properties    = scipy_signal.find_peaks(-triaxus_pressure,prominence=[100,500],distance=100)
bottom_peaks,peak_properties = scipy_signal.find_peaks(triaxus_pressure,prominence=[100,500],distance=100)
top_peaks = np.sort(top_peaks)
bottom_peaks = np.sort(bottom_peaks)
   
dive_index      = top_peaks
surfacing_index = bottom_peaks[bottom_peaks>top_peaks[0]]



land = cfeature.NaturalEarthFeature(category='physical', name='land', scale='50m',
                                    facecolor='gray')
states = cfeature.NaturalEarthFeature(category='cultural', scale='50m', facecolor='none',
                                      name='admin_1_states_provinces_shp')


descending_time             = []
descending_longitude        = []
descending_latitude         = []
    
descending_pressure           = []
descending_depth              = []
descending_abs_salinity_1     = []
descending_con_temperature_1  = []
descending_sigma0_1           = []

descending_abs_salinity_2     = []
descending_con_temperature_2  = []
descending_sigma0_2           = []


ascending_time      = []
ascending_longitude = []
ascending_latitude  = []
ascending_pressure         = []

ascending_abs_salinity_1     = []
ascending_con_temperature_1  = []
ascending_abs_salinity_2     = []
ascending_con_temperature_2  = []
ascending_depth               = []
ascending_sigma0_1           = []
ascending_sigma0_2           = []


ascending_supplementary_vars = {}
descending_supplementary_vars = {}
for i_var in range(0,len(other_variables_to_get)):
    ascending_supplementary_vars[other_variables_to_get[i_var]]  = []
    descending_supplementary_vars[other_variables_to_get[i_var]] = []

n_dive_cycles = len(dive_index)-1

base_date = pandas.Timestamp(2010,1,1)

for i_dive in range(0,n_dive_cycles):
    
    #===============================================# 
    # Descending Profiles    
    #===============================================#
    
    time_for_profile = triaxus_time[dive_index[i_dive]:surfacing_index[i_dive]]
    seconds_since_base_date = (time_for_profile-base_date).values/np.timedelta64(1, 's')

    descending_time.append(seconds_since_base_date)
    descending_longitude.append(triaxus_longitude[dive_index[i_dive]:surfacing_index[i_dive]])
    descending_latitude.append(triaxus_latitude[dive_index[i_dive]:surfacing_index[i_dive]])
    descending_pressure.append(triaxus_pressure[dive_index[i_dive]:surfacing_index[i_dive]])
    

    #Sensor 1
    descending_in_situ_temperature  = triaxus_temperature_1[dive_index[i_dive]:surfacing_index[i_dive]]
    descending_prac_salinity        = triaxus_salinity_1[dive_index[i_dive]:surfacing_index[i_dive]]

 
    absolute_salinity_profile, conservative_temp_profile, depth_profile = compute_CT_SA_z_on_profile(
                                                                          descending_in_situ_temperature, descending_prac_salinity,descending_pressure[-1],descending_longitude[-1],descending_latitude[-1])

    descending_sigma0_1.append(TEOS_10.sigma0(absolute_salinity_profile,conservative_temp_profile))
    descending_abs_salinity_1.append(absolute_salinity_profile)
    descending_con_temperature_1.append(conservative_temp_profile)
    descending_depth.append(depth_profile)
    
    #Sigma0 with iconstant temperature and salinity
#    const_T_profile = conservative_temp_profile 

    #Sensor 2
    descending_in_situ_temperature  = triaxus_temperature_2[dive_index[i_dive]:surfacing_index[i_dive]]
    descending_prac_salinity        = triaxus_salinity_2[dive_index[i_dive]:surfacing_index[i_dive]]

 
    absolute_salinity_profile, conservative_temp_profile, depth_profile = compute_CT_SA_z_on_profile(
                                                                          descending_in_situ_temperature, descending_prac_salinity,descending_pressure[-1],descending_longitude[-1],descending_latitude[-1])
    descending_sigma0_2.append(TEOS_10.sigma0(absolute_salinity_profile,conservative_temp_profile))

    descending_abs_salinity_2.append(absolute_salinity_profile)
    descending_con_temperature_2.append(conservative_temp_profile)
    
    #Supplementary variables (O2, par, all that jazz...)
    for i_var in range(0,len(other_variables_to_get)):
        descending_supplementary_vars[other_variables_to_get[i_var]].append(
                                     triaxus_dataframe[other_variables_to_get[i_var]].values[dive_index[i_dive]:surfacing_index[i_dive]])

    #===============================================# 
    # Ascending Profiles    
    #===============================================#

    time_for_profile = triaxus_time[surfacing_index[i_dive]:dive_index[i_dive+1]] 
    seconds_since_base_date = (time_for_profile-base_date).values/np.timedelta64(1, 's')

    ascending_time.append(seconds_since_base_date)
    
    ascending_longitude.append(triaxus_longitude[surfacing_index[i_dive]:dive_index[i_dive+1]])
    ascending_latitude.append(triaxus_latitude[surfacing_index[i_dive]:dive_index[i_dive+1]])

    ascending_pressure.append(triaxus_pressure[surfacing_index[i_dive]:dive_index[i_dive+1]])
   
    #Sensor 1 
    ascending_in_situ_temperature = triaxus_temperature_1[surfacing_index[i_dive]:dive_index[i_dive+1]]
    ascending_prac_salinity       = triaxus_salinity_1[surfacing_index[i_dive]:dive_index[i_dive+1]]

    absolute_salinity_profile, conservative_temp_profile, depth_profile = compute_CT_SA_z_on_profile(
                                                                          ascending_in_situ_temperature, ascending_prac_salinity,ascending_pressure[-1],ascending_longitude[-1],ascending_latitude[-1])

    ascending_sigma0_1.append(TEOS_10.sigma0(absolute_salinity_profile,conservative_temp_profile))
    ascending_abs_salinity_1.append(absolute_salinity_profile)
    ascending_con_temperature_1.append(conservative_temp_profile)
    ascending_depth.append(depth_profile)

    #Sensor 2 
    ascending_in_situ_temperature = triaxus_temperature_2[surfacing_index[i_dive]:dive_index[i_dive+1]]
    ascending_prac_salinity       = triaxus_salinity_2[surfacing_index[i_dive]:dive_index[i_dive+1]]

    absolute_salinity_profile, conservative_temp_profile, depth_profile = compute_CT_SA_z_on_profile(
                                                                          ascending_in_situ_temperature, ascending_prac_salinity,ascending_pressure[-1],ascending_longitude[-1],ascending_latitude[-1])

    ascending_sigma0_2.append(TEOS_10.sigma0(absolute_salinity_profile,conservative_temp_profile))
    ascending_abs_salinity_2.append(absolute_salinity_profile)
    ascending_con_temperature_2.append(conservative_temp_profile)


    #Supplementary variables (O2, par, all that jazz...)
    for i_var in range(0,len(other_variables_to_get)):
        ascending_supplementary_vars[other_variables_to_get[i_var]].append(
                                     triaxus_dataframe[other_variables_to_get[i_var]].values[surfacing_index[i_dive]:dive_index[i_dive+1]])

max_points_descending_profiles = 0 
max_points_ascending_profiles = 0 

for i_dive in range(0,n_dive_cycles):
        if max_points_descending_profiles<len(descending_depth[i_dive]):
            max_points_descending_profiles = len(descending_depth[i_dive])
        if max_points_ascending_profiles<len(ascending_depth[i_dive]):
            max_points_ascending_profiles = len(ascending_depth[i_dive])


transect_dataset = xarray.DataArray(np.nan * np.zeros([n_dive_cycles,max_points_descending_profiles],dtype=np.float32),name='descending_depth',
                                        dims=['dive_cycle', 'descending_profile_points'], 
                                        coords={'dive_cycle': np.arange(0,n_dive_cycles,1), 'descending_profile_points': np.arange(0,max_points_descending_profiles,1)} ).to_dataset()
transect_dataset['descending_pressure'] = xarray.DataArray(np.nan * np.zeros([n_dive_cycles,max_points_descending_profiles],dtype=np.float32),
                                        dims=['dive_cycle', 'descending_profile_points'], 
                                        coords={'dive_cycle': np.arange(0,n_dive_cycles,1), 'descending_profile_points': np.arange(0,max_points_descending_profiles,1)} )

transect_dataset['descending_longitude'] = xarray.DataArray(np.nan * np.zeros([n_dive_cycles,max_points_descending_profiles],dtype=np.float32),
                                        dims=['dive_cycle', 'descending_profile_points'], 
                                        coords={'dive_cycle': np.arange(0,n_dive_cycles,1), 'descending_profile_points': np.arange(0,max_points_descending_profiles,1)} )
transect_dataset['descending_latitude'] = xarray.DataArray(np.nan * np.zeros([n_dive_cycles,max_points_descending_profiles],dtype=np.float32),
                                        dims=['dive_cycle', 'descending_profile_points'], 
                                        coords={'dive_cycle': np.arange(0,n_dive_cycles,1), 'descending_profile_points': np.arange(0,max_points_descending_profiles,1)} )
    
transect_dataset['descending_time'] = xarray.DataArray(np.zeros([n_dive_cycles,max_points_descending_profiles],dtype=np.float64),
                                        dims=['dive_cycle', 'descending_profile_points'], 
                                        coords={'dive_cycle': np.arange(0,n_dive_cycles,1), 'descending_profile_points': np.arange(0,max_points_descending_profiles,1)},
                                         attrs={"units":"seconds since 2010-01-01"})
    #transect_dataset['descending_time'].units = "seconds since 1950-01-01"
transect_dataset['ascending_pressure'] = xarray.DataArray(np.nan * np.zeros([n_dive_cycles,max_points_ascending_profiles],dtype=np.float32),
                                        dims=['dive_cycle', 'ascending_profile_points'], 
                                        coords={'dive_cycle': np.arange(0,n_dive_cycles,1), 'ascending_profile_points': np.arange(0,max_points_ascending_profiles,1)} )
    
transect_dataset['ascending_depth'] = xarray.DataArray(np.nan * np.zeros([n_dive_cycles,max_points_ascending_profiles],dtype=np.float32),
                                        dims=['dive_cycle', 'ascending_profile_points'], 
                                        coords={'dive_cycle': np.arange(0,n_dive_cycles,1), 'ascending_profile_points': np.arange(0,max_points_ascending_profiles,1)} )
    
transect_dataset['ascending_latitude'] = xarray.DataArray(np.nan * np.zeros([n_dive_cycles,max_points_ascending_profiles],dtype=np.float32),
                                        dims=['dive_cycle', 'ascending_profile_points'], 
                                        coords={'dive_cycle': np.arange(0,n_dive_cycles,1), 'ascending_profile_points': np.arange(0,max_points_ascending_profiles,1)} )
transect_dataset['ascending_longitude'] = xarray.DataArray(np.nan * np.zeros([n_dive_cycles,max_points_ascending_profiles],dtype=np.float32),
                                        dims=['dive_cycle', 'ascending_profile_points'], 
                                        coords={'dive_cycle': np.arange(0,n_dive_cycles,1), 'ascending_profile_points': np.arange(0,max_points_ascending_profiles,1)} )
transect_dataset['ascending_time'] = xarray.DataArray(np.zeros([n_dive_cycles,max_points_ascending_profiles],dtype=np.float64),
                                        dims=['dive_cycle', 'ascending_profile_points'], 
                                        coords={'dive_cycle': np.arange(0,n_dive_cycles,1), 'ascending_profile_points': np.arange(0,max_points_ascending_profiles,1)},
                                        attrs={"units":"seconds since 2010-01-01"})

    
    
transect_dataset['descending_con_temperature_1'] = xarray.DataArray(np.nan * np.zeros([n_dive_cycles,max_points_descending_profiles],dtype=np.float32),name='decending_depth',
                                        dims=['dive_cycle', 'descending_profile_points'], 
                                        coords={'dive_cycle': np.arange(0,n_dive_cycles,1), 'descending_profile_points': np.arange(0,max_points_descending_profiles,1)} )
    
transect_dataset['ascending_con_temperature_1'] = xarray.DataArray(np.nan * np.zeros([n_dive_cycles,max_points_ascending_profiles],dtype=np.float32),
                                        dims=['dive_cycle', 'ascending_profile_points'], 
                                        coords={'dive_cycle': np.arange(0,n_dive_cycles,1), 'ascending_profile_points': np.arange(0,max_points_ascending_profiles,1)} )
    
transect_dataset['descending_abs_salinity_1'] = xarray.DataArray(np.nan * np.zeros([n_dive_cycles,max_points_descending_profiles],dtype=np.float32),name='decending_depth',
                                        dims=['dive_cycle', 'descending_profile_points'], 
                                        coords={'dive_cycle': np.arange(0,n_dive_cycles,1), 'descending_profile_points': np.arange(0,max_points_descending_profiles,1)} )
    
transect_dataset['ascending_abs_salinity_1'] = xarray.DataArray(np.nan * np.zeros([n_dive_cycles,max_points_ascending_profiles],dtype=np.float32),
                                        dims=['dive_cycle', 'ascending_profile_points'], 
                                        coords={'dive_cycle': np.arange(0,n_dive_cycles,1), 'ascending_profile_points': np.arange(0,max_points_ascending_profiles,1)} )
    
transect_dataset['descending_sigma0_1'] = xarray.DataArray(np.nan * np.zeros([n_dive_cycles,max_points_descending_profiles],dtype=np.float32),
                                        dims=['dive_cycle', 'descending_profile_points'], 
                                        coords={'dive_cycle': np.arange(0,n_dive_cycles,1), 'descending_profile_points': np.arange(0,max_points_descending_profiles,1)} ) 
transect_dataset['ascending_sigma0_1'] = xarray.DataArray(np.nan * np.zeros([n_dive_cycles,max_points_ascending_profiles],dtype=np.float32),
                                        dims=['dive_cycle', 'ascending_profile_points'], 
                                        coords={'dive_cycle': np.arange(0,n_dive_cycles,1), 'ascending_profile_points': np.arange(0,max_points_ascending_profiles,1)} )
    
transect_dataset['descending_con_temperature_2'] = xarray.DataArray(np.nan * np.zeros([n_dive_cycles,max_points_descending_profiles],dtype=np.float32),name='decending_depth',
                                        dims=['dive_cycle', 'descending_profile_points'], 
                                        coords={'dive_cycle': np.arange(0,n_dive_cycles,1), 'descending_profile_points': np.arange(0,max_points_descending_profiles,1)} )
    
transect_dataset['ascending_con_temperature_2'] = xarray.DataArray(np.nan * np.zeros([n_dive_cycles,max_points_ascending_profiles],dtype=np.float32),
                                        dims=['dive_cycle', 'ascending_profile_points'], 
                                        coords={'dive_cycle': np.arange(0,n_dive_cycles,1), 'ascending_profile_points': np.arange(0,max_points_ascending_profiles,1)} )
    
transect_dataset['descending_abs_salinity_2'] = xarray.DataArray(np.nan * np.zeros([n_dive_cycles,max_points_descending_profiles],dtype=np.float32),name='decending_depth',
                                        dims=['dive_cycle', 'descending_profile_points'], 
                                        coords={'dive_cycle': np.arange(0,n_dive_cycles,1), 'descending_profile_points': np.arange(0,max_points_descending_profiles,1)} )
    
transect_dataset['ascending_abs_salinity_2'] = xarray.DataArray(np.nan * np.zeros([n_dive_cycles,max_points_ascending_profiles],dtype=np.float32),
                                        dims=['dive_cycle', 'ascending_profile_points'], 
                                        coords={'dive_cycle': np.arange(0,n_dive_cycles,1), 'ascending_profile_points': np.arange(0,max_points_ascending_profiles,1)} )
    
transect_dataset['descending_sigma0_2'] = xarray.DataArray(np.nan * np.zeros([n_dive_cycles,max_points_descending_profiles],dtype=np.float32),
                                        dims=['dive_cycle', 'descending_profile_points'], 
                                        coords={'dive_cycle': np.arange(0,n_dive_cycles,1), 'descending_profile_points': np.arange(0,max_points_descending_profiles,1)} )
    
    
transect_dataset['ascending_sigma0_2'] = xarray.DataArray(np.nan * np.zeros([n_dive_cycles,max_points_ascending_profiles],dtype=np.float32),
                                        dims=['dive_cycle', 'ascending_profile_points'], 
                                        coords={'dive_cycle': np.arange(0,n_dive_cycles,1), 'ascending_profile_points': np.arange(0,max_points_ascending_profiles,1)} )


for i_var in range(0,len(other_variables_to_get)):
    transect_dataset['descending_' + new_var_names[i_var]] = xarray.DataArray(np.nan * np.zeros([n_dive_cycles,max_points_descending_profiles],dtype=np.float32),
                                                      dims=['dive_cycle', 'descending_profile_points'],
                                                      coords={'dive_cycle': np.arange(0,n_dive_cycles,1), 'descending_profile_points': np.arange(0,max_points_descending_profiles,1)} )

    transect_dataset['ascending_' + new_var_names[i_var]] = xarray.DataArray(np.nan * np.zeros([n_dive_cycles,max_points_ascending_profiles],dtype=np.float32),
                                                      dims=['dive_cycle', 'ascending_profile_points'],
                                                      coords={'dive_cycle': np.arange(0,n_dive_cycles,1), 'ascending_profile_points': np.arange(0,max_points_ascending_profiles,1)} )


print(transect_dataset)
    
for i_dive in range(0,n_dive_cycles):
    n_profile_points = len(descending_depth[i_dive])
    transect_dataset['descending_depth'][i_dive,0:n_profile_points]            = descending_depth[i_dive]
    transect_dataset['descending_pressure'][i_dive,0:n_profile_points]         = descending_pressure[i_dive]
    transect_dataset['descending_longitude'][i_dive,0:n_profile_points]        = descending_longitude[i_dive]
    transect_dataset['descending_latitude'][i_dive,0:n_profile_points]         = descending_latitude[i_dive]
    transect_dataset['descending_time'][i_dive,0:n_profile_points]             = descending_time[i_dive]
        #transect_dataset['descending_time'][i_dive,n_profile_points::]             = np.datetime64("NaT")
        

    transect_dataset['descending_con_temperature_1'][i_dive,0:n_profile_points] = descending_con_temperature_1[i_dive]
    transect_dataset['descending_abs_salinity_1'][i_dive,0:n_profile_points]    = descending_abs_salinity_1[i_dive]
    transect_dataset['descending_sigma0_1'][i_dive,0:n_profile_points]          = descending_sigma0_1[i_dive]

    transect_dataset['descending_con_temperature_2'][i_dive,0:n_profile_points] = descending_con_temperature_2[i_dive]
    transect_dataset['descending_abs_salinity_2'][i_dive,0:n_profile_points]    = descending_abs_salinity_2[i_dive]
    transect_dataset['descending_sigma0_2'][i_dive,0:n_profile_points]          = descending_sigma0_2[i_dive]
    
    for i_var in range(0,len(other_variables_to_get)):
        transect_dataset['descending_' + new_var_names[i_var]][i_dive,0:n_profile_points] = descending_supplementary_vars[other_variables_to_get[i_var]][i_dive]

    n_profile_points = len(ascending_depth[i_dive])
    transect_dataset['ascending_depth'][i_dive,0:n_profile_points]           = ascending_depth[i_dive]
    transect_dataset['ascending_pressure'][i_dive,0:n_profile_points]        = ascending_pressure[i_dive]
    transect_dataset['ascending_longitude'][i_dive,0:n_profile_points]       = ascending_longitude[i_dive]
    transect_dataset['ascending_latitude'][i_dive,0:n_profile_points]        = ascending_latitude[i_dive]
    transect_dataset['ascending_time'][i_dive,0:n_profile_points]            = ascending_time[i_dive]
        #transect_dataset['ascending_time'][i_dive,n_profile_points::]            = np.datetime64("NaT")

    transect_dataset['ascending_con_temperature_1'][i_dive,0:n_profile_points] = ascending_con_temperature_1[i_dive]
    transect_dataset['ascending_abs_salinity_1'][i_dive,0:n_profile_points]    = ascending_abs_salinity_1[i_dive]
    transect_dataset['ascending_sigma0_1'][i_dive,0:n_profile_points]          = ascending_sigma0_1[i_dive]

    
    transect_dataset['ascending_con_temperature_2'][i_dive,0:n_profile_points] = ascending_con_temperature_2[i_dive]
    transect_dataset['ascending_abs_salinity_2'][i_dive,0:n_profile_points]    = ascending_abs_salinity_2[i_dive]
    transect_dataset['ascending_sigma0_2'][i_dive,0:n_profile_points]          = ascending_sigma0_2[i_dive]


    for i_var in range(0,len(other_variables_to_get)):
        transect_dataset['ascending_' + new_var_names[i_var]][i_dive,0:n_profile_points]  = ascending_supplementary_vars[other_variables_to_get[i_var]][i_dive]

processed_data_path = ''

output_file_name = CNV_FILE_NAME.split('.cnv')[0] + '_processed_ascending_descending.nc'

print('Writing data to: ', os.path.join(OUTPUT_DATA_PATH,OUTPUT_FILE_NAME))
transect_dataset.to_netcdf(os.path.join(OUTPUT_DATA_PATH,OUTPUT_FILE_NAME)) 


