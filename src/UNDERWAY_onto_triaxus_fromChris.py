# -*- coding: utf-8 -*-
"""
NOT WORKING!!!!!!!!!!!!!
Created on Thu May 13 19:51:51 2021

@author: Chris
"""

import numpy as np
import xarray
import os
#from matplotlib import pyplot as plt
#import scipy.signal as scipy_signal
#import scipy.stats as scipy_stats
#import geopy
#import cmocean
import pandas
#import xarray

VOYAGE_ID = 'in2022_v06'
TRIAXUS_TOW = '01_002'

#TRIAXUS_DATA_PATH = './'
#OUTPUT_PATH      = './'
#OUTPUT_FILE_NAME = VOYAGE_ID + '_triaxus' + str('%02d' % TRIAXUS_TOW) + '_SADCP' + str(SADCP_FREQUENCY) +  '.nc'
INPUT_DATA_PATH = '../DATA/TRIAXUS/Processing/Processed_Data/TOW_1'
INPUT_FILE_NAME = VOYAGE_ID + '_' +  TRIAXUS_TOW + '_processed_ascending_descending.nc'
#INPUT_FILE_NAME = VOYAGE_ID + '_triaxus_tow_' + TRIAXUS_TOW + '_profiles.nc'

processed_data_path = './'
output_file_name = VOYAGE_ID + '_triaxus_' + TRIAXUS_TOW + '_UWAY.nc'

# Get triaxus
#TRIAXUS_FILE_NAME = VOYAGE_ID + '_triaxus_tow' + str('%02d' % TRIAXUS_TOW) + '_profiles.nc'
processed_triaxus_dataset = xarray.open_dataset(
                                os.path.join(INPUT_DATA_PATH,INPUT_FILE_NAME),
                                decode_times=True,use_cftime=False)
#processed_triaxus_dataset = xarray.open_dataset(
#                                os.path.join(TRIAXUS_DATA_PATH,TRIAXUS_FILE_NAME),
#                                decode_times=True,use_cftime=False)
processed_triaxus_dataset['time'] = pandas.to_datetime(processed_triaxus_dataset['time'].values)
n_dive_cycles = processed_triaxus_dataset['time'].size


# Get underway
#ds_underw_all = xar.open_dataset(r'..\..\underway\in2022_v06uwy.nc')
underway_data_path = '..\..\underway'
underway_file_name = VOYAGE_ID + 'uwy.nc'

underway_dataset = xarray.open_dataset(
                                os.path.join(underway_data_path,underway_file_name),
                                decode_times=True,use_cftime=False)

epoch = pandas.Timestamp(underway_dataset.Epoch[-22:-3])
start_time = epoch + pandas.to_timedelta(int(underway_dataset.Epoch[0:8]),unit='S')
uwy_date_time       = start_time + pandas.to_timedelta(5.0*underway_dataset['sample'].values,unit='S')
underway_dataset    = underway_dataset.reindex(sample=uwy_date_time)
underway_dataset    = underway_dataset.rename({'sample':'time'})

print(underway_dataset['time'][0])

#underway_samples = underway_dataset['sample']


ship_longitude = underway_dataset['longitude']
ship_latitude  = underway_dataset['latitude']



variables_to_get = ['portAirTemp','stbdAirTemp','portTrueWindSpeed','portTrueWindDir','stbdTrueWindSpeed','stbdTrueWindDir',
                    'atmPressure','ultrasonicTrueWindSpeed','ultrasonicTrueWindDir','portRain','stbdRain','portRadiometer','stbdRadiometer',
                    'portPAR','stbdPAR','salinity','waterTemp','fluorescence']

#variables_to_get = ['portAirTemp','stbdAirTemp','portTrueWindSpeed','portTrueWindDir','stbdTrueWindSpeed','stbdTrueWindDir','portDewPoint','stbdDewPoint',
#                    'atmPressure','ultrasonicTrueWindSpeed','ultrasonicTrueWindDir','portRain','stbdRain','portRadiometer','stbdRadiometer',
#                    'portPAR','stbdPAR','salinity','waterTemp','fluorescence']
longitude   = np.zeros(2*n_dive_cycles,dtype=np.float64)
latitude    = np.zeros(2*n_dive_cycles,dtype=np.float64)
time        = np.empty(2*n_dive_cycles,dtype='datetime64[ns]')

#Allocate memory
variables_to_save = []
for i_var in range(0,len(variables_to_get)):
    variables_to_save.append(np.zeros(2*n_dive_cycles,dtype=np.float64))


underway_dataset.close()

for i_dive in range(0,n_dive_cycles):
    
    print("Dive: ", i_dive)
    descending_longitude   = processed_triaxus_dataset['longitude'][i_dive,:].squeeze(drop=True).dropna(dim='descending_profile_points')
    descending_latitude    = processed_triaxus_dataset['latitude'][i_dive,:].squeeze(drop=True).dropna(dim='descending_profile_points')
    descending_time        = pandas.to_datetime(processed_triaxus_dataset['time'][i_dive,:].squeeze(drop=True).dropna(dim='descending_profile_points').values)
    descending_time        = descending_time[time.year>2010]
   

    uwy_for_descent     = underway_dataset.sel(time=slice(descending_time[0],descending_time[-1])).mean(dim='time')

    n_profile_points    = descending_time.size
    longitude[2*i_dive] = descending_longitude[int(n_profile_points/2)].values
    latitude[2*i_dive]  = descending_latitude[int(n_profile_points/2)].values
    time[2*i_dive]      = descending_time[int(n_profile_points/2)]

    for i_var in range(0,len(variables_to_get)):
        variables_to_save[i_var][2*i_dive] = uwy_for_descent[variables_to_get[i_var]].values
            

    ascending_longitude   = processed_triaxus_dataset['ascending_longitude'][i_dive,:].squeeze(drop=True).dropna(dim='ascending_profile_points')
    ascending_latitude    = processed_triaxus_dataset['ascending_latitude'][i_dive,:].squeeze(drop=True).dropna(dim='ascending_profile_points')
    ascending_time       = pandas.to_datetime(processed_triaxus_dataset['ascending_time'][i_dive,:].squeeze(drop=True).dropna(dim='ascending_profile_points').values)
    ascending_time       = ascending_time[ascending_time.year>2010]

    uwy_for_ascent     = underway_dataset.sel(time=slice(ascending_time[0],ascending_time[-1])).mean(dim='time')
   
    n_profile_points        = ascending_time.size
    longitude[(2*i_dive)+1] = ascending_longitude[int(n_profile_points/2)].values
    latitude[(2*i_dive)+1]  = ascending_latitude[int(n_profile_points/2)].values
    time[(2*i_dive)+1]      = ascending_time[int(n_profile_points/2)]

 
    for i_var in range(0,len(variables_to_get)):
        variables_to_save[i_var][2*i_dive+1] = uwy_for_ascent[variables_to_get[i_var]].values
 
transect_output_dataset = xarray.DataArray(longitude,name='longitude',
                                           dims=['time'],coords={'time': time}).to_dataset()
transect_output_dataset['latitude'] = xarray.DataArray(latitude,dims=['time'],coords={'time': time})
for i_var in range(0,len(variables_to_get)):
    transect_output_dataset[variables_to_get[i_var]] = xarray.DataArray(variables_to_save[i_var],dims=['time'],coords={'time': time})


transect_output_dataset.to_netcdf(os.path.join(processed_data_path,output_file_name))

