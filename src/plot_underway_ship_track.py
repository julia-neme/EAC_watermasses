###############################################################################
# Plot all data from underway ADCP with map of positions. 
# -----------------------------------
# 16/07/2022 by Julia Neme and Hannah Dawson (in2022_v06)
###############################################################################

import cmocean as cm
import datetime as dt
import glob
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import xarray as xr
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
plt.rcParams['font.size'] = 14

adcp_75 = xr.open_dataset(r"../../adcp/uhdas/proc/os75nb/contour/os75nb.nc")
adcp_150 = xr.open_dataset(r"../../adcp/uhdas/proc/os150nb/contour/os150nb.nc")
N75 = len(adcp_75['time'])
N150 = len(adcp_150['time'])
# Stop script if both ADCP's have different recording length!
if N75 != N150:
   print('Different recording times of ADCP 75Hz and ADCP 150Hz')

N = np.min([N75, N150])
# Get date labels
label = []
for n in np.arange(0, N-(N%100)+100, 100):
      label.append(str(adcp_75['time'].values[n])[5:10])
axlabel = ['ADCP_150 U-vel', 'ADCP_150 V-vel', 'ADCP_75 U-vel', 'ADCP_75 V-vel']

def figure_map():
   coast = np.loadtxt("BATHY/eaccoast.dat")
   bathy = xr.open_dataset('BATHY/bathy_dbdb2_v30_AUSTRALIA.nc')
   lon_bathy = bathy['lon']
   lat_bathy = bathy['lat']
   depth = bathy['height']

   fig = plt.figure(figsize = (25,8))
   grd = gs.GridSpec(2, 4, hspace = .15, width_ratios = [1,1,1,.05])
   axs = [fig.add_subplot(grd[:,0]),
         fig.add_subplot(grd[0,1]),
         fig.add_subplot(grd[0,2]),
         fig.add_subplot(grd[1,1]),
         fig.add_subplot(grd[1,2]),
         fig.add_subplot(grd[:,3])]
   axs[0].plot(coast[:,0], coast[:,1], 'k')
   c = axs[0].contour(lon_bathy, lat_bathy, depth, 
                     levels = np.arange(-6000,0,500),
                     cmap = cm.cm.gray, zorder = 0)
   axs[0].contourf(lon_bathy, lat_bathy, depth.where(depth>0), 
                  colors = 'peachpuff', linestyles = 'solid')
   axs[0].set_xlim(152, 158)
   axs[0].set_ylim(-30, -25)

   axs[1].set_ylim(-310, 0)
   axs[2].set_ylim(-310, 0)
   axs[3].set_ylim(-610, 0)
   axs[4].set_ylim(-610, 0)
   
   count = 0
   for ax in axs[1:-1]:
      ax.set_xticks(adcp_75.time[np.arange(0, N-(N%100)+100, 100)])
      ax.set_title(axlabel[count])
      count += 1
   for ax in axs[1:3]:
      ax.set_xticklabels([])
   for ax in axs[3:-1]:
      ax.set_xticklabels(label, rotation = 20)
   fig.suptitle('IN2022_V06 Underway ADCP', fontsize = 18)

   return fig, axs

fig, axs = figure_map()
c = axs[0].scatter(adcp_75['lon'].values, adcp_75['lat'].values, 
                  c = np.arange(0, N), cmap = 'gnuplot')
cbar1 = plt.colorbar(c, orientation = 'horizontal', ax = axs[0])
cbar1.set_ticks(np.arange(0, N-(N%100)+100, 100))
cbar1.ax.set_xticklabels(label, rotation = 20)

j = axs[1].contourf(adcp_150['time'], -adcp_150['depth'].isel(time=0), adcp_150['u'].T,
               levels = np.arange(-.6, .65, .05), cmap = 'RdBu_r', 
               extend = 'both')
axs[2].contourf(adcp_150['time'], -adcp_150['depth'].isel(time=0), adcp_150['v'].T,
               levels = np.arange(-.6, .65, .05), cmap = 'RdBu_r', 
               extend = 'both')      

axs[3].contourf(adcp_75['time'], -adcp_75['depth'].isel(time=0), adcp_75['u'].T,
               levels = np.arange(-.6, .65, .05), cmap = 'RdBu_r', 
               extend = 'both')
axs[4].contourf(adcp_75['time'], -adcp_75['depth'].isel(time=0), adcp_75['v'].T,
                  levels = np.arange(-.6, .65, .05), cmap = 'RdBu_r', 
                  extend = 'both')  
cbar2 = plt.colorbar(j, orientation = 'vertical', cax = axs[-1])
cbar2.set_label('Velocity (m/s)')
         
plt.savefig('../FIGURES/IN2022_underway_adcp_alltimes.jpg', bbox_inches = 'tight')
