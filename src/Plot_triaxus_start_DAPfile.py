# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 22:18:12 2022

@author: Amandine
"""
from matplotlib import pyplot as plt
import numpy as np
#import datetime
#import pandas as pd # import libraries
import xarray as xar
import cmocean
import matplotlib

ds_underw_all = xar.open_dataset(r'..\..\underway\in2022_v06uwy.nc')

VOYAGE_ID = 'in2022_v06'
# Tow1
TRIAXUS_TOW = '01_002'
TRIAXUS_file_name = '../../triaxus/processing/in2022_v06/cap/cappro/avg/section/tow1Leg2AvgCast.nc'
# Tow 2
#TRIAXUS_TOW = '02_002'
#TRIAXUS_file_name = '../../triaxus/processing/in2022_v06/cap/cappro/avg/section/tow2Leg2AvgCast.nc'
# Tow 2
#TRIAXUS_TOW = '02_002'
#TRIAXUS_file_name = '../../triaxus/processing/in2022_v06/cap/cappro/avg/section/tow2Leg2AvgCast.nc'

plot_folder = '../FIGURES/'
#TRIAXUS_file_name = '../DATA/TRIAXUS/Processing/Processed_Data/TOW_1/in2022_v06_triaxus_tow_' + TRIAXUS_TOW + '_profiles.nc'


ds_triaxus = xar.open_dataset(TRIAXUS_file_name)
# temperature
#salinity
# chlorophyll
# oxygen
# cdom
# nitrate

# Parameters;
lon_min = np.floor(np.min(ds_triaxus.longitude)*10)/10
lon_max = np.floor(np.max(ds_triaxus.longitude)*10+1)/10


# Figures
plt.figure()
#plt.contourf(ds_triaxus.time,ds_triaxus.pressure,ds_triaxus.temperature.T,20,cmap = cmocean.cm.thermal)
plt.pcolor(ds_triaxus.longitude,ds_triaxus.pressure,ds_triaxus.temperature,cmap = cmocean.cm.thermal)
plt.gca().invert_yaxis()
# LABELS / LIMITS
plt.ylim(400,0)     # plt.ylim(-30.6, -29.8)     
plt.xlim(lon_min,lon_max)     # plt.ylim(-30.6, -29.8)     
plt.title('Conservative temperature Triaxus tow' + TRIAXUS_TOW,fontsize=16)       
#plt.xlabel('Longitude [$^o$E]',fontsize=13)                  
plt.ylabel('Depth [m]',fontsize=13)                  
cb = plt.colorbar()
cb.set_label('Temperature [$^{o}$C]', size=13)
#CS = plt.contour(ds_triaxus_SADCP150.longitude,ds_triaxus_SADCP150.depth,ds_triaxus_SADCP150.v.T,cmap = matplotlib.colors.ListedColormap(['k']), levels=(-1, -0.5, 0),linewidths = 0.3)
#plt.clabel(CS, fontsize=10, fmt='%0.1f')
# SAVE
plt.tight_layout()
plt.savefig(plot_folder + 'plot_triaxus_tow' + TRIAXUS_TOW + '_temp.png', bbox_inches='tight', pad_inches=0.5, dpi=300)
plt.show()


plt.figure()
#plt.contourf(ds_triaxus.longitude,ds_triaxus.pressure,ds_triaxus.salinity.T,20,cmap = cmocean.cm.haline)
plt.pcolor(ds_triaxus.longitude,ds_triaxus.pressure,ds_triaxus.salinity,cmap = cmocean.cm.haline)
plt.gca().invert_yaxis()
# LABELS / LIMITS
plt.ylim(400,0)     # plt.ylim(-30.6, -29.8)     
plt.xlim(lon_min,lon_max)     # plt.ylim(-30.6, -29.8)     
plt.title('Absolute salinity Triaxus tow' + TRIAXUS_TOW,fontsize=16)             
plt.xlabel('Longitude [$^o$E]',fontsize=13)                  
plt.ylabel('Depth [m]',fontsize=13)                  
cb = plt.colorbar()
cb.set_label('Salinity', size=13)
#CS = plt.contour(ds_triaxus_SADCP150.longitude,ds_triaxus_SADCP150.depth,ds_triaxus_SADCP150.v.T,cmap = matplotlib.colors.ListedColormap(['k']), levels=(-1, -0.5, 0),linewidths = 0.3)
#plt.clabel(CS, fontsize=10, fmt='%0.1f')
#plt.clim(clim0)
# SAVE
plt.tight_layout()
plt.savefig(plot_folder + 'plot_triaxus_tow' + TRIAXUS_TOW + '_sal.png', bbox_inches='tight', pad_inches=0.5, dpi=300)
plt.show()


plt.figure()
#plt.contourf(ds_triaxus.longitude,ds_triaxus.pressure,ds_triaxus.salinity.T,20,cmap = cmocean.cm.haline)
plt.pcolor(ds_triaxus.longitude,ds_triaxus.pressure,ds_triaxus.oxygen,cmap = cmocean.cm.oxy)
plt.gca().invert_yaxis()
# LABELS / LIMITS
plt.ylim(400,0)     # plt.ylim(-30.6, -29.8)     
plt.xlim(lon_min,lon_max)     # plt.ylim(-30.6, -29.8)     
plt.title('Dissolved oxygen Triaxus tow' + TRIAXUS_TOW,fontsize=16)             
plt.xlabel('Longitude [$^o$E]',fontsize=13)                  
plt.ylabel('Depth [m]',fontsize=13)                  
cb = plt.colorbar()
cb.set_label('Oxygen [uM]', size=13)
#CS = plt.contour(ds_triaxus_SADCP150.longitude,ds_triaxus_SADCP150.depth,ds_triaxus_SADCP150.v.T,cmap = matplotlib.colors.ListedColormap(['k']), levels=(-1, -0.5, 0),linewidths = 0.3)
#plt.clabel(CS, fontsize=10, fmt='%0.1f')
#plt.clim(clim0)
# SAVE
plt.tight_layout()
plt.savefig(plot_folder + 'plot_triaxus_tow' + TRIAXUS_TOW + '_o2.png', bbox_inches='tight', pad_inches=0.5, dpi=300)
plt.show()


plt.figure()
#plt.contourf(ds_triaxus.longitude,ds_triaxus.pressure,ds_triaxus.salinity.T,20,cmap = cmocean.cm.haline)
plt.pcolor(ds_triaxus.longitude,ds_triaxus.pressure,ds_triaxus.chlorophyll,cmap = cmocean.cm.algae)
plt.gca().invert_yaxis()
# LABELS / LIMITS
plt.ylim(400,0)     # plt.ylim(-30.6, -29.8)     
plt.xlim(lon_min,lon_max)     # plt.ylim(-30.6, -29.8)     
plt.title('Chlorophyll Triaxus tow ' + TRIAXUS_TOW,fontsize=16)             
plt.xlabel('Longitude [$^o$E]',fontsize=13)                  
plt.ylabel('Depth [m]',fontsize=13)                  
cb = plt.colorbar()
cb.set_label('Chlorophyll [ug/L]', size=13)
#CS = plt.contour(ds_triaxus_SADCP150.longitude,ds_triaxus_SADCP150.depth,ds_triaxus_SADCP150.v.T,cmap = matplotlib.colors.ListedColormap(['k']), levels=(-1, -0.5, 0),linewidths = 0.3)
#plt.clabel(CS, fontsize=10, fmt='%0.1f')
#plt.clim(clim0)
# SAVE
plt.tight_layout()
plt.savefig(plot_folder + 'plot_triaxus_tow' + TRIAXUS_TOW + '_chll.png', bbox_inches='tight', pad_inches=0.5, dpi=300)
plt.show()


plt.figure()
#plt.contourf(ds_triaxus.longitude,ds_triaxus.pressure,ds_triaxus.salinity.T,20,cmap = cmocean.cm.haline)
plt.pcolor(ds_triaxus.longitude,ds_triaxus.pressure,ds_triaxus.nitrate,cmap = cmocean.cm.matter)
plt.gca().invert_yaxis()
# LABELS / LIMITS
plt.ylim(400,0)     # plt.ylim(-30.6, -29.8)     
plt.xlim(lon_min,lon_max)     # plt.ylim(-30.6, -29.8)     
plt.title('Nitrate Triaxus tow ' + TRIAXUS_TOW,fontsize=16)             
plt.xlabel('Longitude [$^o$E]',fontsize=13)                  
plt.ylabel('Depth [m]',fontsize=13)                  
cb = plt.colorbar()
cb.set_label('Nitrate [umole]', size=13)
#CS = plt.contour(ds_triaxus_SADCP150.longitude,ds_triaxus_SADCP150.depth,ds_triaxus_SADCP150.v.T,cmap = matplotlib.colors.ListedColormap(['k']), levels=(-1, -0.5, 0),linewidths = 0.3)
#plt.clabel(CS, fontsize=10, fmt='%0.1f')
#plt.clim(clim0)
# SAVE
plt.tight_layout()
plt.savefig(plot_folder + 'plot_triaxus_tow' + TRIAXUS_TOW + '_nitrate.png', bbox_inches='tight', pad_inches=0.5, dpi=300)
plt.show()