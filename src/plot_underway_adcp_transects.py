###############################################################################
# Plot all data from underway ADCP with map of positions. 
# -----------------------------------
# 17/07/2022 by Julia Neme and Hannah Dawson (in2022_v06)
###############################################################################

from pyexpat.errors import XML_ERROR_DUPLICATE_ATTRIBUTE
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
plt.rcParams['font.size'] = 13

adcp_75 = xr.open_dataset(r"../../adcp/uhdas/proc/os75nb/contour/os75nb.nc")
adcp_150 = xr.open_dataset(r"../../adcp/uhdas/proc/os150nb/contour/os150nb.nc")
N = len(adcp_150['time'])

coast = np.loadtxt("BATHY/eaccoast.dat")
bathy = xr.open_dataset('BATHY/bathy_dbdb2_v30_AUSTRALIA.nc')
lon_bathy = bathy['lon']
lat_bathy = bathy['lat']
depth = bathy['height']

def figure_map(transect, ti, tf, xlabels):
   
   fig = plt.figure(figsize = (25,10))
   grd = gs.GridSpec(2, 4, hspace = .3, width_ratios = [1,1,1,.05])
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
                  colors = 'gainsboro', linestyles = 'solid')
   axs[0].scatter(153, -27.5, marker = 'o', c = 'dimgrey', s = 80, zorder = 2)
   axs[0].text(152.1, -27.8, 'Brisbane', c = 'dimgrey', fontsize = 15)
    
   axs[0].set_xlim(152, 158)
   axs[0].set_ylim(-30, -25)
   axs[1].set_ylim(-310, 0)
   axs[2].set_ylim(-310, 0)
   axs[3].set_ylim(-610, 0)
   axs[4].set_ylim(-610, 0)
   
   count = 0
   for ax in axs[1:-1]:
      ax.set_xticks(adcp_150['time'][np.arange(0, N-(N%20)+20, 20)])
      ax.set_title(axlabel[count])
      ax.set_xticklabels(xlabels, rotation = 45, fontsize=11)
      count += 1
   fig.suptitle('IN2022_V06 Underway ADCP Transect '+str(transect)+'\n Start time: '+ 
   str(adcp_75['time'].values[ti])[5:16] + ', End time: '+
   str(adcp_75['time'].values[tf])[5:16], fontsize = 18)

   return fig, axs

tis = [0, 240, 332, 440, 558, 715, 840, 1099]
tfs = [168, 275, 375, 555, 594, 765, 875, 1130]
tr = [1, -1, 1, -1, 1, -1, 1, 1]

for transect in [1]:#, 2, 3, 4, 5]:
    ti = tis[transect-1]
    tf = tfs[transect-1]
    N = tf - ti
    tbool = tr[transect-1]
    label = []
    xlabels = []
    if transect == 5:
        lvls = np.arange(-.8, .85, .05)
    else:
        lvls = np.arange(-.6, .65, .05)
    for n in np.arange(0, N-(N%20)+20, 20):
        label.append(str(adcp_75['time'].values[n])[5:16])
        xlabels.append(adcp_75['lon'].values[n])
    xlabels = np.round(xlabels, decimals = 1)
    axlabel = ['ADCP_150 U-vel', 'ADCP_150 V-vel', 'ADCP_75 U-vel', 'ADCP_75 V-vel', 'ADCP_150 Speed','ADCP_75 Speed']

    fig, axs = figure_map(transect, ti, tf, xlabels)
    c = axs[0].scatter(adcp_75['lon'][ti:tf].values, adcp_75['lat'][ti:tf].values,
            c = np.arange(0, N), cmap = 'rainbow')
    cbar1 = plt.colorbar(c, orientation = 'horizontal', ax = axs[0])
    cbar1.set_ticks(np.arange(0, N-(N%20)+20, 20))
    cbar1.ax.set_xticklabels(label, rotation = 90)
    axs[1].contourf(adcp_150['time'][ti:tf], -adcp_150['depth'].isel(time=0), 
                    adcp_150['u'][ti:tf].T, levels = lvls, 
                    cmap = 'RdBu_r', extend = 'both')
    axs[2].contourf(adcp_150['time'][ti:tf], -adcp_150['depth'].isel(time=0), 
                    adcp_150['v'][ti:tf].T, levels = lvls, 
                    cmap = 'RdBu_r', extend = 'both')      
    axs[3].contourf(adcp_75['time'][ti:tf], -adcp_75['depth'].isel(time=0), 
                    adcp_75['u'][ti:tf].T, levels = lvls, 
                    cmap = 'RdBu_r', extend = 'both')
    c = axs[4].contourf(adcp_75['time'][ti:tf], -adcp_75['depth'].isel(time=0), 
                        adcp_75['v'][ti:tf].T, levels = lvls, 
                        cmap = 'RdBu_r', extend = 'both')  
    cbar2 = plt.colorbar(c, orientation = 'vertical', cax = axs[-1])
    cbar2.set_label('Velocity (m/s)')   
    # Add time references
    for ax in axs[3:5]:
        ax.scatter(adcp_75['time'][ti:tf], 
                -np.zeros(N)*adcp_75['depth'].isel(time = 0)[0].values,
                c = np.arange(0, N), s = 180, cmap = 'rainbow')
        if tbool == -1:
            ax.invert_xaxis()
    for ax in axs[1:3]:
        ax.scatter(adcp_150['time'][ti:tf], 
                -np.zeros(N)*adcp_150['depth'].isel(time = 0)[0].values,
                c = np.arange(0, N), s = 180, cmap = 'rainbow')
        if tbool == -1:
            ax.invert_xaxis()
    plt.savefig('../FIGURES/IN2022_underway_adcp_transect'+str(transect)+'.jpg', 
                bbox_inches = 'tight')
    plt.show()


ti=1210
for tf in range(1300, 1350, 10):

    plt.figure(figsize = (15,15))
    plt.plot(coast[:,0], coast[:,1], 'k')
    c = plt.contour(lon_bathy, lat_bathy, depth, 
                    levels = np.arange(-6000,0,500),
                    cmap = cm.cm.gray, zorder = 0)
    plt.contourf(lon_bathy, lat_bathy, depth.where(depth>0), 
                colors = 'gainsboro', linestyles = 'solid')
    plt.scatter(153, -27.5, marker = 'o', c = 'dimgrey', s = 80, zorder = 2)
    plt.text(152.5, -27.7, 'Brisbane', c = 'dimgrey', fontsize = 15)
        
    plt.scatter(adcp_150['lon'][ti:tf].values, adcp_150['lat'][ti:tf].values, 
                c = np.arange(0, tf-ti), cmap = 'rainbow')
    plt.xlim(153, 156)
    plt.ylim(-28, -26)
    plt.show()
plt.close()







# testing with speed


def figure_map(transect, ti, tf, xlabels):
   
   fig = plt.figure(figsize = (30,9))
   grd = gs.GridSpec(2, 7, hspace = .3, wspace=.22, width_ratios = [1,1,1,1,.05,.01,.05])
   axs = [fig.add_subplot(grd[:,0]),
         fig.add_subplot(grd[0,1]),
         fig.add_subplot(grd[0,2]),
         fig.add_subplot(grd[0,3]),
         fig.add_subplot(grd[1,1]),
         fig.add_subplot(grd[1,2]),
         fig.add_subplot(grd[1,3]),
         fig.add_subplot(grd[:,4]),
         fig.add_subplot(grd[:,6])]
   axs[0].plot(coast[:,0], coast[:,1], 'k')
   c = axs[0].contour(lon_bathy, lat_bathy, depth, 
                     levels = np.arange(-6000,0,500),
                     cmap = cm.cm.gray, zorder = 0)
   axs[0].contourf(lon_bathy, lat_bathy, depth.where(depth>0), 
                  colors = 'gainsboro', linestyles = 'solid')
   axs[0].scatter(153, -27.5, marker = 'o', c = 'dimgrey', s = 80, zorder = 2)
   axs[0].text(152.1, -27.8, 'Brisbane', c = 'dimgrey', fontsize = 15)
    
   axs[0].set_xlim(152, 158)
   axs[0].set_ylim(-30, -25)
   axs[1].set_ylim(-310, 0)
   axs[2].set_ylim(-310, 0)
   axs[3].set_ylim(-310, 0)
   axs[4].set_ylim(-610, 0)
   axs[5].set_ylim(-610, 0)
   axs[6].set_ylim(-610, 0)
   
   count = 0
   for ax in axs[1:-2]:
      ax.set_xticks(adcp_150['time'][np.arange(0, N-(N%20)+20, 20)])
      ax.set_title(axlabel[count])
      ax.set_xticklabels(xlabels, rotation = 45, fontsize=11)
      #ax.set_yticklabels(fontsize=12)
      if count >=3:
          ax.set_xlabel('Longitude ($^\circ$E)')
      count += 1
   fig.suptitle('IN2022_V06 Underway ADCP Transect '+str(transect)+'\n Start time: '+ 
   str(adcp_75['time'].values[ti])[5:16] + ', End time: '+
   str(adcp_75['time'].values[tf])[5:16], fontsize = 18)

   return fig, axs




transect = 1
ti = tis[transect-1]
tf = tfs[transect-1]
N = tf - ti
tbool = tr[transect-1]
label = []
xlabels = []
if transect == 5:
    lvls = np.arange(-.8, .85, .05)
else:
    lvls = np.arange(-.6, .65, .05)
spdlvls = np.arange(0, 1, .1)
for n in np.arange(0, N-(N%20)+20, 20):
    label.append(str(adcp_75['time'].values[n])[5:16])
    xlabels.append(adcp_75['lon'].values[n])
xlabels = np.round(xlabels, decimals = 1)
#axlabel = ['ADCP_150 U-vel', 'ADCP_150 V-vel', 'ADCP_75 U-vel', 'ADCP_75 V-vel']
axlabel = ['ADCP_150 U-vel', 'ADCP_150 V-vel', 'ADCP_150 Speed', 'ADCP_75 U-vel', 'ADCP_75 V-vel' ,'ADCP_75 Speed']
yy = str(adcp_75['time'].values[0])[0:4]
dd = str(adcp_75['time'].values[0])[5:7]
mm = str(adcp_75['time'].values[0])[8:10]


fig, axs = figure_map(transect, ti, tf, xlabels)
c = axs[0].scatter(adcp_75['lon'][ti:tf].values, adcp_75['lat'][ti:tf].values,
        c = np.arange(0, N), cmap = 'rainbow')
cbar1 = plt.colorbar(c, orientation = 'horizontal', ax = axs[0])
cbar1.set_ticks(np.arange(0, N-(N%20)+20, 20))
cbar1.ax.set_xticklabels(label, rotation = 90)
adcp_150['speed'] = np.sqrt(adcp_150.u**2 + adcp_150.v**2)
adcp_75['speed'] = np.sqrt(adcp_75.u**2 + adcp_75.v**2)

axs[1].contourf(adcp_150['time'][ti:tf], -adcp_150['depth'].isel(time=0), 
                adcp_150['u'][ti:tf].T, levels = lvls, 
                cmap = 'RdBu_r', extend = 'both')
axs[2].contourf(adcp_150['time'][ti:tf], -adcp_150['depth'].isel(time=0), 
                adcp_150['v'][ti:tf].T, levels = lvls, 
                cmap = 'RdBu_r', extend = 'both')
axs[3].contourf(adcp_150['time'][ti:tf], -adcp_150['depth'].isel(time=0), 
                adcp_150['speed'][ti:tf].T, levels = spdlvls, 
                cmap = cm.cm.speed, extend = 'max')                   
axs[4].contourf(adcp_75['time'][ti:tf], -adcp_75['depth'].isel(time=0), 
                adcp_75['u'][ti:tf].T, levels = lvls, 
                cmap = 'RdBu_r', extend = 'both')
c = axs[5].contourf(adcp_75['time'][ti:tf], -adcp_75['depth'].isel(time=0), 
                    adcp_75['v'][ti:tf].T, levels = lvls, 
                    cmap = 'RdBu_r', extend = 'both')  
sp = axs[6].contourf(adcp_75['time'][ti:tf], -adcp_75['depth'].isel(time=0), 
                adcp_75['speed'][ti:tf].T, levels = spdlvls, 
                cmap = cm.cm.speed, extend = 'max')    
cbar2 = plt.colorbar(c, orientation = 'vertical', cax = axs[-2])
cbar2.set_label('Velocity (m/s)')  
cbar3 = plt.colorbar(sp, orientation = 'vertical', cax = axs[-1])
cbar3.set_label('Speed (m/s)')   
# Add time references
for ax in axs[4:7]:
    ax.scatter(adcp_75['time'][ti:tf], 
            -np.zeros(N)*adcp_75['depth'].isel(time = 0)[0].values,
            c = np.arange(0, N), s = 180, cmap = 'rainbow')
    if tbool == -1:
        ax.invert_xaxis()
for ax in axs[1:4]:
    ax.scatter(adcp_150['time'][ti:tf], 
            -np.zeros(N)*adcp_150['depth'].isel(time = 0)[0].values,
            c = np.arange(0, N), s = 180, cmap = 'rainbow')
    if tbool == -1:
        ax.invert_xaxis()
plt.savefig('../FIGURES/IN2022_underway_adcp_transect'+str(transect)+'_trial.jpg', 
            bbox_inches = 'tight')
plt.show()

