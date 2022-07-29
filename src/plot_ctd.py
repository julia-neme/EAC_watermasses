###############################################################################
# This script makes a single plot for all CTD casts saved as .nc files listed
# in the path_to_files directory (using function figure). It also creates a 
# figure with CTD station location and T and S profiles (using figure_map 
# function). Run from command line as python3 plot_ctd.py
# -----------------------------------
# 13/07/2022 by Julia Neme (in2022_v06)
###############################################################################

import cmocean as cm
import glob
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import xarray as xr
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
plt.rcParams['font.size'] = 14

path_to_files = glob.glob(r"../../ctd/processing/in2022_v06/cap/cappro/avg/*.nc")

def figure(path):

    data = xr.open_dataset(path)
    lon = np.round(data['longitude'].squeeze().values, 2)
    lat = np.round(data['latitude'].squeeze().values, 2)
    tim = str(data['time'].values[0])[:-10]
    ttl = path[-15:-9]

    fig = plt.figure(figsize = (18,10))
    ax1 = fig.add_subplot(131)
    ax1s = ax1.twiny()
    ax2 = fig.add_subplot(132, sharey = ax1)
    ax3 = fig.add_subplot(133, sharey = ax1)

    fig.suptitle('CAST: '+ttl+'  TIME: '+tim+' \n LON: '+str(lon)+'; LAT: '
                 +str(lat), fontsize = 16)
    ax1.set_ylabel('Pressure (dbar)'); 
    ax1.yaxis.set_minor_locator(MultipleLocator(100))
    ax1.yaxis.set_minor_locator(MultipleLocator(50))
    ax1.set_xlim(0, 26)
    ax1s.set_xlim(34.4, 35.8)
    ax2.set_xlim(140, 270)
    #ax3.set_xlim(0, 3)

    p = data['pressure'].squeeze()
    # Temperature
    if 'temperature' in data.data_vars:
        t = data['temperature'].squeeze()
        ax1.plot(t, p, color = 'C0')
        ax1.set_xlabel('Temperature ($^{\circ}$C)', color = 'C0')
        ax1.tick_params(axis = 'x', color = 'C0', labelcolor = 'C0')
        ax1.set_ylim(0, data['pressure'][-1])
        ax1.axes.invert_yaxis()
        ax1.grid(linestyle = '--', color = 'lightgrey')
    if 'salinity' in data.data_vars:
        s = data['salinity'].squeeze()
        ax1s.plot(s, p, color = 'orangered')
        ax1s.set_xlabel('Salinity (PSU)', color = 'orangered')
        ax1s.tick_params(axis = 'x', color = 'orangered', labelcolor = 'orangered')
        ax1s.set_ylim(0, data['pressure'][-1])
        ax1s.axes.invert_yaxis()
    if 'oxygen' in data.data_vars:
        o = data['oxygen'].squeeze()
        ax2.plot(o, p, color = 'k')
        ax2.set_xlabel('Oxygen (uM)', color = 'k')
        ax2.tick_params(axis = 'x', color = 'k', labelcolor = 'k')
        ax2.set_ylim(0, data['pressure'][-1])
        ax2.axes.invert_yaxis()
        ax2.grid(linestyle = '--', color = 'lightgrey')
    if 'fluorometerFlag' in data.data_vars:
        mask = data['fluorometerFlag'].squeeze() == 0
        f = data['fluorometer'].squeeze()
        ax3.plot(f[mask], p[mask], color = 'C2')
        ax3.set_xlabel('Fluorescence (ug/l)', color = 'C2')
        ax3.tick_params(axis = 'x', color = 'C2', labelcolor = 'C2')
        ax3.set_ylim(0, data['pressure'][-1])
        ax3.axes.invert_yaxis()
        ax3.grid(linestyle = '--', color = 'lightgrey')
    plt.savefig('../FIGURES/ctd_cast'+ttl+'.jpg', bbox_inches = 'tight')
    plt.close()
    return lon, lat, ttl, t, s, p

def figure_map():
    coast = np.loadtxt("BATHY/eaccoast.dat")
    bathy = xr.open_dataset('BATHY/bathy_dbdb2_v30_AUSTRALIA.nc')
    lon_bathy = bathy['lon']
    lat_bathy = bathy['lat']
    depth = bathy['height']

    fig = plt.figure(figsize = (20,11))
    grd = gs.GridSpec(2, 3, height_ratios = [1,.3])
    axs = [fig.add_subplot(grd[0,0]),
           fig.add_subplot(grd[:,1]),
           fig.add_subplot(grd[:,2])]
    axs[0].plot(coast[:,0], coast[:,1], 'k')
    c = axs[0].contour(lon_bathy, lat_bathy, depth, 
                       levels = np.arange(-6000,0,500),
                       cmap = cm.cm.gray, zorder = 0)
    axs[0].contourf(lon_bathy, lat_bathy, depth.where(depth>0), 
                    colors = 'peachpuff', linestyles = 'solid')

    axs[0].set_xlim(152, 156)
    axs[0].set_ylim(-30, -25)
    axs[1].set_xlim(0, 26)
    axs[1].set_xlabel('Temperature ($^{o}$C)')
    axs[2].set_xlim(34.4, 35.8)
    axs[2].set_xlabel('Salinity (PSU)')
    axs[1].set_ylim(0,4000)
    axs[1].set_ylabel('Pressure (dbar)')
    axs[2].set_ylim(0,4000)
    axs[1].grid(linestyle = '--', color = 'lightgrey')
    axs[2].grid(linestyle = '--', color = 'lightgrey')
    fig.suptitle('INV2022_V06 CTD casts', fontsize = 16)

    return fig, axs

locs = np.empty([len(path_to_files), 2])
labl = []
cmap = plt.get_cmap('gnuplot')
ls = np.linspace(0,1,len(path_to_files))
n = 0
for cast in path_to_files: 
    locs[n,0], locs[n,1], l, t, s, p = figure(cast)
    labl.append(l)
    n += 1

fig, axs = figure_map(); 
n = 0
for cast in path_to_files:
    locs[n,0], locs[n,1], l, t, s, p = figure(cast);
    labl.append(l)
    axs[0].scatter(locs[n,0], locs[n,1], color = cmap(ls[n]), 
                   label = labl[n]);
    axs[1].plot(t, p, color = cmap(ls[n]));
    axs[2].plot(s, p, color = cmap(ls[n]));
    n += 1
axs[1].axes.invert_yaxis()
axs[2].axes.invert_yaxis()
axs[0].legend(bbox_to_anchor = (.3,-.4,.4,.4), loc = 'lower center', ncol = 3, 
              fontsize = 8)
plt.savefig('../FIGURES/INV2022_V06_CTD_all.jpg', bbox_inches = 'tight')
