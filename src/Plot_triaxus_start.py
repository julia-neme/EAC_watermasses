# -*- coding: utf-8 -*-
"""
Created on Wed May 12 19:09:00 2021

@author: Amandine
"""
from matplotlib import pyplot as plt
import numpy as np
import datetime
#import pandas as pd # import libraries
import xarray as xar
import cmocean
import matplotlib

ds_underw_all = xar.open_dataset(r'..\..\underway\in2022_v06uwy.nc')

INPUT_DATA_PATH = '../DATA/TRIAXUS/Processing/Processed_Data/TOW_1'
VOYAGE_ID = 'in2022_v06'
TRIAXUS_TOW = '01_002'

plot_folder = '../FIGURES/'
TRIAXUS_file_name = '../DATA/TRIAXUS/Processing/Processed_Data/TOW_1/in2022_v06_triaxus_tow_' + TRIAXUS_TOW + '_profiles.nc'
SADCP_file_name = '../DATA/TRIAXUS/Processing/Processed_Data/TOW_1/in2022_v06_triaxus' + TRIAXUS_TOW + '_SADCP75.nc' # TO CHANGE
SADCP_file_name2 = '../DATA/TRIAXUS/Processing/Processed_Data/TOW_1/in2022_v06_triaxus' + TRIAXUS_TOW + '_SADCP150.nc'
#TRIAXUS_file_name = 'in2021_v03_triaxus_tow03_all_profiles.nc'
#SADCP_file_name = 'in2021_v03_triaxus03_all_SADCP75.nc'
#SADCP_file_name2 = 'in2021_v03_triaxus03_all_SADCP150.nc'

ds_triaxus = xar.open_dataset(TRIAXUS_file_name)
ds_triaxus_SADCP75 = xar.open_dataset(SADCP_file_name)
ds_triaxus_SADCP150 = xar.open_dataset(SADCP_file_name2)

# Parameters;
lon_min = np.floor(np.min(ds_triaxus.longitude)*10)/10
lon_max = np.floor(np.max(ds_triaxus.longitude)*10+1)/10

# Figures
plt.figure()
plt.contourf(ds_triaxus.longitude,ds_triaxus.pressure,ds_triaxus.con_temperature_1.T,20,cmap = cmocean.cm.thermal)
plt.gca().invert_yaxis()
# LABELS / LIMITS
plt.ylim(400,0)     # plt.ylim(-30.6, -29.8)     
plt.xlim(lon_min,lon_max)     # plt.ylim(-30.6, -29.8)     
plt.title('Conservative temperature Triaxus tow' + TRIAXUS_TOW,fontsize=16)       
plt.xlabel('Longitude [$^o$E]',fontsize=13)                  
plt.ylabel('Depth [m]',fontsize=13)                  
cb = plt.colorbar()
cb.set_label('Temperature [$^{o}$C]', size=13)
CS = plt.contour(ds_triaxus_SADCP150.longitude,ds_triaxus_SADCP150.depth,ds_triaxus_SADCP150.v.T,cmap = matplotlib.colors.ListedColormap(['k']), levels=(-1, -0.5, 0),linewidths = 0.3)
plt.clabel(CS, fontsize=10, fmt='%0.1f')
# SAVE
plt.tight_layout()
plt.savefig(plot_folder + 'plot_triaxus_tow' + TRIAXUS_TOW + '_temp.png', bbox_inches='tight', pad_inches=0.5, dpi=300)
plt.show()



plt.figure()
plt.contourf(ds_triaxus.longitude,ds_triaxus.pressure,ds_triaxus.abs_salinity_1.T,20,cmap = cmocean.cm.haline)
plt.gca().invert_yaxis()
# LABELS / LIMITS
plt.ylim(400,0)     # plt.ylim(-30.6, -29.8)     
plt.xlim(lon_min,lon_max)     # plt.ylim(-30.6, -29.8)     
plt.title('Absolute salinity Triaxus tow' + TRIAXUS_TOW,fontsize=16)             
plt.xlabel('Longitude [$^o$E]',fontsize=13)                  
plt.ylabel('Depth [m]',fontsize=13)                  
cb = plt.colorbar()
cb.set_label('Salinity', size=13)
CS = plt.contour(ds_triaxus_SADCP150.longitude,ds_triaxus_SADCP150.depth,ds_triaxus_SADCP150.v.T,cmap = matplotlib.colors.ListedColormap(['k']), levels=(-1, -0.5, 0),linewidths = 0.3)
plt.clabel(CS, fontsize=10, fmt='%0.1f')
#plt.clim(clim0)
# SAVE
plt.tight_layout()
plt.savefig(plot_folder + 'plot_triaxus_tow' + TRIAXUS_TOW + '_sal.png', bbox_inches='tight', pad_inches=0.5, dpi=300)
plt.show()

cmap = cmocean.cm.balance;
plt.figure()
plt.contourf(ds_triaxus.longitude,ds_triaxus.pressure,ds_triaxus.sigma0_1.T,20,cmap = cmocean.cm.dense)
plt.gca().invert_yaxis()
# LABELS / LIMITS
plt.ylim(400,0)     # plt.ylim(-30.6, -29.8)     
plt.xlim(lon_min,lon_max)     # plt.ylim(-30.6, -29.8)     
plt.title('Density Triaxus tow' + TRIAXUS_TOW ,fontsize=16)           
plt.xlabel('Longitude [$^o$E]',fontsize=13)                  
plt.ylabel('Depth [m]',fontsize=13)                  
cb = plt.colorbar()
cb.set_label('Density [kg m$^{-3}$]', size=13)
CS = plt.contour(ds_triaxus_SADCP150.longitude,ds_triaxus_SADCP150.depth,ds_triaxus_SADCP150.v.T,cmap = matplotlib.colors.ListedColormap(['k']), levels=(-1, -0.5, 0),linewidths = 0.3)
plt.clabel(CS, fontsize=10, fmt='%0.1f')

#plt.clim(clim0)
# SAVE
plt.tight_layout()
plt.savefig(plot_folder + 'plot_triaxus_tow' + TRIAXUS_TOW + '_dens.png', bbox_inches='tight', pad_inches=0.5, dpi=300)
plt.show()





cmap = cmocean.cm.balance;
levelsu = np.arange(-1, 1, 0.1)
levelsv = np.arange(-1, 1, 0.1)

plt.figure()
Q = plt.contourf(ds_triaxus_SADCP75.longitude,ds_triaxus_SADCP75.depth,ds_triaxus_SADCP75.u.T,20,cmap = cmap, levels=levelsu)
#plt.contour(ds_triaxus_SADCP75.longitude,ds_triaxus_SADCP75.depth,ds_triaxus_SADCP75.u.T,cmap = cmap, levels=(0),colour = 'k')
plt.gca().invert_yaxis()
# LABELS / LIMITS
plt.ylim(800,0)     # plt.ylim(-30.6, -29.8)     
plt.xlim(lon_min,lon_max)     # plt.ylim(-30.6, -29.8)     
plt.title('Zonal velocity Triaxus tow' + TRIAXUS_TOW ,fontsize=16)          
plt.xlabel('Longitude [$^o$E]',fontsize=13)                  
plt.ylabel('Depth [m]',fontsize=13)                  
cb = plt.colorbar(Q)
cb.set_label('Velocity u [m s$^{-1}$]', size=13)
#plt.clim(clim0)
# SAVE
plt.tight_layout()
plt.savefig(plot_folder + 'plot_triaxus_tow' + TRIAXUS_TOW + '_U.png', bbox_inches='tight', pad_inches=0.5, dpi=300)
plt.show()

plt.figure()
plt.contourf(ds_triaxus_SADCP75.longitude,ds_triaxus_SADCP75.depth,ds_triaxus_SADCP75.v.T,20,cmap = cmap, levels=levelsv)
plt.gca().invert_yaxis()
# LABELS / LIMITS
plt.ylim(800,0)     # plt.ylim(-30.6, -29.8)     
plt.xlim(lon_min,lon_max)     # plt.ylim(-30.6, -29.8)     
plt.title('Meridional velocity Triaxus tow' + TRIAXUS_TOW ,fontsize=16)         
plt.xlabel('Longitude [$^o$E]',fontsize=13)                  
plt.ylabel('Depth [m]',fontsize=13)                  
cb = plt.colorbar()
cb.set_label('Velocity v [m s$^{-1}$]', size=13)
# SAVE
plt.tight_layout()
plt.savefig(plot_folder + 'plot_triaxus_tow' + TRIAXUS_TOW + '_V.png', bbox_inches='tight', pad_inches=0.5, dpi=300)
plt.show()


ds_triaxus_SADCP75.close()
ds_triaxus_SADCP150.close()

ds_triaxus.close()


#
### Select time
##speed = np.sqrt(ds_triaxus_SADCP75.u**2+ds_triaxus_SADCP75.v**2)
##ds_triaxus_SADCP75['speed'] = speed
###aaa =ds_triaxus_SADCP75.sel(time="2021-05-11 00")
##ds_triaxus_SADCP75.speed.sel(time="2021-05-10 21").mean(dim='time').plot()
##ds_triaxus_SADCP75.speed.sel(time="2021-05-11 05").mean(dim='time').plot()
##
##ds_triaxus_SADCP75.u.sel(time="2021-05-10 21").mean(dim='time').plot()
##ds_triaxus_SADCP75.u.sel(time="2021-05-11 05").mean(dim='time').plot()
##ds_triaxus_SADCP75.v.sel(time="2021-05-10 21").mean(dim='time').plot()
##ds_triaxus_SADCP75.v.sel(time="2021-05-11 05").mean(dim='time').plot()
##
##ds_triaxus_SADCP75.longitude.sel(time="2021-05-10 21")
##ds_triaxus_SADCP75.longitude.sel(time="2021-05-11 05")
#
## Or latitude
#ind_lon1 = np.where(ds_triaxus_SADCP75.longitude>151.8)[0][0]
#ind_lon2 = np.where(ds_triaxus_SADCP75.longitude>152.1)[0][0]
#ind_lon3 = np.where(ds_triaxus_SADCP75.longitude>152.4)[0][0]
#
#f, ax = plt.subplots(1,2,figsize=(12,12))
#ds_triaxus_SADCP75.u.isel(time = ind_lon1).plot.line(y="depth", ax=ax[0],label='lon1 151.8')
#ds_triaxus_SADCP75.v.isel(time = ind_lon1).plot.line(y="depth", ax=ax[1],label='lon1 151.8')
#ds_triaxus_SADCP75.u.isel(time = ind_lon2).plot.line(y="depth", ax=ax[0],label='lon1 152.1')
#ds_triaxus_SADCP75.v.isel(time = ind_lon2).plot.line(y="depth", ax=ax[1],label='lon1 152.1')
#ds_triaxus_SADCP75.u.isel(time = ind_lon3).plot.line(y="depth", ax=ax[0],label='lon1 152.4')
#ds_triaxus_SADCP75.v.isel(time = ind_lon3).plot.line(y="depth", ax=ax[1],label='lon1 152.4')
#ax[0].legend()
#ax[1].invert_yaxis()
#ax[0].invert_yaxis()
## SAVE
#plt.tight_layout()
#plt.savefig(plot_folder + 'plot_triaxus_tow02_profiles_UV.png', bbox_inches='tight', pad_inches=0.5, dpi=300)
#plt.show()
#
#f, ax = plt.subplots(1,2,figsize=(12,12))
#ds_triaxus.con_temperature_1.isel(time = ind_lon1).plot.line(y="pressure", ax=ax[0],label='lon1 151.8')
#ds_triaxus.abs_salinity_1.isel(time = ind_lon1).plot.line(y="pressure", ax=ax[1],label='lon1 151.8')
#ds_triaxus.con_temperature_1.isel(time = ind_lon2).plot.line(y="pressure", ax=ax[0],label='lon1 152.1')
#ds_triaxus.abs_salinity_1.isel(time = ind_lon2).plot.line(y="pressure", ax=ax[1],label='lon1 152.1')
#ds_triaxus.con_temperature_1.isel(time = ind_lon3).plot.line(y="pressure", ax=ax[0],label='lon1 152.4')
#ds_triaxus.abs_salinity_1.isel(time = ind_lon3).plot.line(y="pressure", ax=ax[1],label='lon1 152.4')
#ax[0].legend()
#ax[1].invert_yaxis()
#ax[0].invert_yaxis()
## SAVE
#plt.tight_layout()
#plt.savefig(plot_folder + 'plot_triaxus_tow02_profilesTS.png', bbox_inches='tight', pad_inches=0.5, dpi=300)
#plt.show()




