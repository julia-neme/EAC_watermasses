import numpy as np
import xarray 
import os
#from matplotlib import pyplot as plt
#import scipy.signal as scipy_signal
#import scipy.stats as scipy_stats
#import geopy
#import cmocean
import pandas

#==================HEADER STUFF - User should pay attention===============================================================#
SADCP_FREQUENCY = 75
VOYAGE_ID = 'in2022_v06'
TRIAXUS_TOW = '01_002'


INPUT_DATA_PATH = '../DATA/TRIAXUS/Processing/Processed_Data/TOW_1'
OUTPUT_DATA_PATH = '../DATA/TRIAXUS/Processing/Processed_Data/TOW_1'
SADCP_FILE_PATH = r'..\..\adcp\uhdas\proc\os' + str(SADCP_FREQUENCY) + 'nb/contour'

#OUTPUT_FILE_NAME = VOYAGE_ID + '_triaxus' + str('%02d' % TRIAXUS_TOW) + '_SADCP' + str(SADCP_FREQUENCY) +  '.nc'
OUTPUT_FILE_NAME = VOYAGE_ID + '_triaxus' + TRIAXUS_TOW + '_SADCP' + str(SADCP_FREQUENCY) +  '.nc'

#==========================================================================================================================#

SADCP_file_name = 'os' + str(SADCP_FREQUENCY) + 'nb.nc'
SADCP_dataset = xarray.open_dataset(os.path.join(SADCP_FILE_PATH,SADCP_file_name),decode_times=True,use_cftime=False)
SADCP_dataset['time'] = pandas.to_datetime(SADCP_dataset['time'].values)
SADCP_depth           = SADCP_dataset['depth'][0,:].squeeze()
print('Read SADCP file: ', SADCP_file_name + ' from ' + SADCP_FILE_PATH)


TRIAXUS_DATA_PATH = '../DATA/TRIAXUS/Processing/Processed_Data/TOW_1'
TRIAXUS_FILE_NAME = VOYAGE_ID + '_triaxus_tow_' + TRIAXUS_TOW + '_profiles.nc'

processed_triaxus_dataset = xarray.open_dataset(
                                os.path.join(TRIAXUS_DATA_PATH,TRIAXUS_FILE_NAME),
                                decode_times=True,use_cftime=False)
processed_triaxus_dataset['time'] = pandas.to_datetime(processed_triaxus_dataset['time'].values)


print('Read Triaxus file: ', TRIAXUS_FILE_NAME + ' from ' + TRIAXUS_DATA_PATH)
n_profiles = processed_triaxus_dataset['time'].size


longitude   = np.zeros(n_profiles,dtype=np.float64)
latitude    = np.zeros(n_profiles,dtype=np.float64)
time        = np.empty(n_profiles,dtype='datetime64[ns]')

u_profile   = np.zeros([n_profiles,SADCP_depth.size],dtype=np.float64)
v_profile   = np.zeros([n_profiles,SADCP_depth.size],dtype=np.float64)
pg_profile   = np.zeros([n_profiles,SADCP_depth.size],dtype=np.float64)


for i_profile in range(0,n_profiles):


    triaxus_time_stamp       = processed_triaxus_dataset['time'][i_profile]


    SADCP_for_profile        = SADCP_dataset.interp(time=triaxus_time_stamp,method='nearest')
    u_profile[i_profile,:]   = SADCP_for_profile['u']
    v_profile[i_profile,:]   = SADCP_for_profile['v']

    pg_profile[i_profile,:] = SADCP_for_profile['pg']
    longitude[i_profile]    = SADCP_for_profile['lon']

    latitude[i_profile]    = SADCP_for_profile['lat']
    time[i_profile]        = pandas.to_datetime(SADCP_for_profile['time'].values)



u_profile[np.abs(u_profile)>10] = np.nan
v_profile[np.abs(v_profile)>10] = np.nan

processed_triaxus_dataset.close()
SADCP_dataset.close()


transect_output_dataset = xarray.DataArray(longitude,name='longitude',
                                           dims=['time'],coords={'time': time}).to_dataset()
transect_output_dataset['latitude'] = xarray.DataArray(latitude,dims=['time'],coords={'time': time})
transect_output_dataset['u']     = xarray.DataArray(u_profile,dims=['time','depth'],coords={'time': time,'depth':SADCP_depth.values})
transect_output_dataset['v']     = xarray.DataArray(v_profile,dims=['time','depth'],coords={'time': time,'depth':SADCP_depth.values})
transect_output_dataset['pg']    = xarray.DataArray(pg_profile,dims=['time','depth'],coords={'time': time,'depth':SADCP_depth.values})


print('Writing file: ' + OUTPUT_DATA_PATH, ' to path ', OUTPUT_FILE_NAME)
transect_output_dataset.to_netcdf(os.path.join(OUTPUT_DATA_PATH,OUTPUT_FILE_NAME))



