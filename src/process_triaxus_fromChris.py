import numpy as np
import xarray 
import os
#from matplotlib import pyplot as plt
#import scipy.signal as scipy_signal
import scipy.stats as scipy_stats
import pandas
#import gsw as TEOS_10
import cartopy.feature as cfeature


#=============================== Header: Users should modify input paths and files ================================#
#VOYAGE_ID = 'in2021_v03'
#TRIAXUS_TOW = '03_004'
VOYAGE_ID = 'in2022_v06'
TRIAXUS_TOW = '01_002'


#INPUT_DATA_PATH = './'
##INPUT_FILE_NAME = VOYAGE_ID + '_' +  str('%02d' % TRIAXUS_TOW) + '_002_processed_ascending_descending.nc'
#INPUT_FILE_NAME = VOYAGE_ID + '_' +  TRIAXUS_TOW + '_processed_ascending_descending.nc'
#OUTPUT_DATA_PATH = './'
#OUTPUT_FILE_NAME = VOYAGE_ID + '_triaxus_tow' + TRIAXUS_TOW + '_profiles.nc'


INPUT_DATA_PATH = '../DATA/TRIAXUS/Processing/Processed_Data/TOW_1'
INPUT_FILE_NAME = VOYAGE_ID + '_' +  TRIAXUS_TOW + '_processed_ascending_descending.nc'

CNV_FILE_NAME = 'in2022_v06_01_002.cnv'
OUTPUT_DATA_PATH = '../DATA/TRIAXUS/Processing/Processed_Data/TOW_1'
OUTPUT_FILE_NAME = VOYAGE_ID + '_triaxus_tow_' + TRIAXUS_TOW + '_profiles.nc'



#===================================================================================================================#
#processed_triaxus_dataset = xarray.open_dataset(
#                                os.path.join(INPUT_DATA_PATH,INPUT_FILE_NAME),
#                                decode_times=True,use_cftime=False)
processed_triaxus_dataset = xarray.open_dataset(
                                os.path.join(INPUT_DATA_PATH,INPUT_FILE_NAME),
                                decode_times=True,use_cftime=False)


n_dive_cycles = processed_triaxus_dataset['dive_cycle'].size
land = cfeature.NaturalEarthFeature(category='physical', name='land', scale='50m',
                                    facecolor='gray')
states = cfeature.NaturalEarthFeature(category='cultural', scale='50m', facecolor='none',
                                      name='admin_1_states_provinces_shp')



pressure_bin_size = 2 #Pressure bin size in db
pressure_max = 400
pressure_min = 0

pressure_bins  = np.arange(pressure_min,pressure_max,pressure_bin_size)

total_data_vars = processed_triaxus_dataset.data_vars

variables_to_get = []

not_data_variables = ['profile_points','longitude','latitude','pressure','time']
for i_var in total_data_vars:

    #print(i_var)

    if ('descending' in i_var) and not any(x in i_var for x in not_data_variables):
        variables_to_get.append(i_var.replace('descending_',''))


descending_to_ascending_string = ['descending','ascending']
binned_profile_data = {}

for i_up_or_down in descending_to_ascending_string:
    for i_var in variables_to_get:
        binned_profile_data[i_var] = np.zeros([2*n_dive_cycles,pressure_bins.size-1],dtype=np.float64)


#dyn_hgt_gridded     = np.zeros([depth_bins.size-1,n_dive_cycles],dtype=np.float64)
#spice_gridded       = np.zeros([depth_bins.size-1,n_dive_cycles],dtype=np.float64)

longitude   = np.zeros(2*n_dive_cycles,dtype=np.float64)
latitude    = np.zeros(2*n_dive_cycles,dtype=np.float64)
time        = np.empty(2*n_dive_cycles,dtype='datetime64[ns]')

dive_counter = 0
for i_dive in range(0,n_dive_cycles):

    print('Dive: ', i_dive, ' of ', n_dive_cycles)
    for i_up_or_down in descending_to_ascending_string:
        print('Processing ' + i_up_or_down + ' profiles')
        
        current_scan_number = processed_triaxus_dataset[i_up_or_down + '_profile_points']
        current_longitude   = processed_triaxus_dataset[i_up_or_down + '_longitude'][i_dive,:].squeeze(drop=True) #.dropna(dim='descending_profile_points')
        current_latitude    = processed_triaxus_dataset[i_up_or_down + '_latitude'][i_dive,:].squeeze(drop=True) #.dropna(dim='descending_profile_points')

        current_pressure    = processed_triaxus_dataset[i_up_or_down + '_pressure'][i_dive,:].squeeze(drop=True) #.dropna(dim='descending_profile_points')
        current_time        = pandas.to_datetime(processed_triaxus_dataset[i_up_or_down + '_time'][i_dive,:].squeeze(drop=True).values) #.dropna(dim='descending_profile_points').values)
        current_time        = current_time[current_time.year>2010]

        current_profiles = []
    
        n_profile_points = current_time.size

        longitude[dive_counter] = current_longitude[int(n_profile_points/2)].values
        latitude[dive_counter]  = current_latitude[int(n_profile_points/2)].values
        time[dive_counter]      = current_time[int(n_profile_points/2)]
    
        for i_var in variables_to_get:
            current_variable = processed_triaxus_dataset[i_up_or_down + '_' + i_var][i_dive,:].squeeze(drop=True) #.dropna(dim='descending_profile_points')
   
            binned_profile,_,_ = scipy_stats.binned_statistic(current_pressure, current_variable, statistic='mean', bins=pressure_bins, range=None)
            binned_profile_data[i_var][dive_counter,:] = binned_profile
        
        dive_counter = dive_counter+1


processed_triaxus_dataset.close() 


transect_output_dataset = xarray.DataArray(longitude,name='longitude',
                                           dims=['time'],coords={'time': time}).to_dataset()
transect_output_dataset['latitude'] = xarray.DataArray(latitude,dims=['time'],coords={'time': time})

for i_var in variables_to_get:
    transect_output_dataset[i_var] = xarray.DataArray(binned_profile_data[i_var],dims=['time','pressure'],coords={'time': time,'pressure':pressure_bins[0:pressure_bins.size-1]})    


print('Writing file:', OUTPUT_FILE_NAME + ' to path: ' + OUTPUT_FILE_NAME)

transect_output_dataset.to_netcdf(os.path.join(OUTPUT_DATA_PATH , OUTPUT_FILE_NAME))


