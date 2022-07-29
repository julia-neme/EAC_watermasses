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

t1l2 = xr.open_dataset(r'../../triaxus/processing/in2022_v06/cap/cappro/avg/section/tow1Leg2AvgCast.nc')
t2l2 = xr.open_dataset(r'../../triaxus/processing/in2022_v06/cap/cappro/avg/section/tow2Leg2AvgCast.nc')
t3l2 = xr.open_dataset(r'../../triaxus/processing/in2022_v06/cap/cappro/avg/section/tow3Leg2AvgCast.nc')
t4l2 = xr.open_dataset(r'../../triaxus/processing/in2022_v06/cap/cappro/avg/section/tow4Leg2AvgCast.nc')
t4l4 = xr.open_dataset(r'../../triaxus/processing/in2022_v06/cap/cappro/avg/section/tow4Leg4AvgCast.nc')
t5l2 = xr.open_dataset(r'../../triaxus/processing/in2022_v06/cap/cappro/avg/section/tow5Leg2AvgCast.nc')



def figure():
    fig = plt.figure(figsize = (20,16))
    grd = gs.GridSpec(3, 2, hspace = .15)
    axs = [fig.add_subplot(grd[0,0]),
           fig.add_subplot(grd[0,1]),
           fig.add_subplot(grd[1,0]),
           fig.add_subplot(grd[1,1]),
           fig.add_subplot(grd[2,0])]
    for ax in [axs[1], axs[3]]:
        ax.set_yticklabels([])
    for ax in [axs[0], axs[2], axs[4]]:
        ax.set_ylabel('Pressure (dbar)')
    for ax in axs[-2:]:
        ax.set_xlabel('Longitude ($^{\circ}$E)')
    return fig, axs

# Full fields

var = ['temperature', 'salinity', 'oxygen', 'chlorophyll', 'nitrate']
cmaps = [cm.cm.thermal, cm.cm.haline, cm.cm.oxy, 'YlGnBu_r', cm.cm.matter, ]

fig, axs = figure()
n = 0
for ax in axs:
    if var[n] != 'temperature':
        units = t1l2[var[n]].units
    else:
        units = '$^{\circ}$C'
    c = axs[n].pcolormesh(t1l2['longitude'], -t1l2['pressure'],
                          t1l2[var[n]], cmap = cmaps[n])
    plt.colorbar(c, ax = axs[n]).set_label(var[n]+' ('+units+')')
    n += 1
fig.suptitle('IN2022_V06 - TRIAXUS tow 1', x = 0.5, y = 0.95, fontsize = 20)
plt.savefig('../FIGURES/triaxus_tow1_in2022_v06.jpg', bbox_inches = 'tight')

fig, axs = figure()
n = 0
for ax in axs:
    if var[n] != 'temperature':
        units = t2l2[var[n]].units
    else:
        units = '$^{\circ}$C'
    c = axs[n].pcolormesh(t2l2['longitude'], -t2l2['pressure'],
                          t2l2[var[n]], cmap = cmaps[n])
    plt.colorbar(c, ax = axs[n]).set_label(var[n]+' ('+units+')')
    n += 1
fig.suptitle('IN2022_V06 - TRIAXUS tow 2', x = 0.5, y = 0.95, fontsize = 20)
plt.savefig('../FIGURES/triaxus_tow2_in2022_v06.jpg', bbox_inches = 'tight')
plt.show()

fig, axs = figure()
n = 0
for ax in axs:
    if var[n] != 'temperature':
        units = t3l2[var[n]].units
    else:
        units = '$^{\circ}$C'
    c = axs[n].pcolormesh(t3l2['longitude'], -t3l2['pressure'],
                          t3l2[var[n]], cmap = cmaps[n])
    plt.colorbar(c, ax = axs[n]).set_label(var[n]+' ('+units+')')
    n += 1
fig.suptitle('IN2022_V06 - TRIAXUS tow 3', x = 0.5, y = 0.95, fontsize = 20)
plt.savefig('../FIGURES/triaxus_tow3_in2022_v06.jpg', bbox_inches = 'tight')
plt.show()

fig, axs = figure()
n = 0
for ax in axs:
    if var[n] != 'temperature':
        units = t4l2[var[n]].units
    else:
        units = '$^{\circ}$C'
    c = axs[n].pcolormesh(t4l2['longitude'], -t4l2['pressure'],
                          t4l2[var[n]], cmap = cmaps[n])
    plt.colorbar(c, ax = axs[n]).set_label(var[n]+' ('+units+')')
    n += 1
fig.suptitle('IN2022_V06 - TRIAXUS tow 4 leg 2', x = 0.5, y = 0.95, fontsize = 20)
plt.savefig('../FIGURES/triaxus_tow4l2_in2022_v06.jpg', bbox_inches = 'tight')
plt.show()

fig, axs = figure()
n = 0
for ax in axs:
    if var[n] != 'temperature':
        units = t4l4[var[n]].units
    else:
        units = '$^{\circ}$C'
    c = axs[n].pcolormesh(t4l4['latitude'], -t4l4['pressure'],
                          t4l4[var[n]], cmap = cmaps[n])
    plt.colorbar(c, ax = axs[n]).set_label(var[n]+' ('+units+')')
    n += 1
    ax.set_xlabel('Latitude ($^{\circ}$N)')
fig.suptitle('IN2022_V06 - TRIAXUS tow 4 leg 4', x = 0.5, y = 0.95, fontsize = 20)
plt.savefig('../FIGURES/triaxus_tow4l3_in2022_v06.jpg', bbox_inches = 'tight')
plt.show()

fig, axs = figure()
n = 0
for ax in axs:
    if var[n] != 'temperature':
        units = t5l2[var[n]].units
    else:
        units = '$^{\circ}$C'
    c = axs[n].pcolormesh(t5l2['longitude'], -t5l2['pressure'],
                          t5l2[var[n]], cmap = cmaps[n])
    plt.colorbar(c, ax = axs[n]).set_label(var[n]+' ('+units+')')
    n += 1
fig.suptitle('IN2022_V06 - TRIAXUS tow 5 leg 2', x = 0.5, y = 0.95, fontsize = 20)
plt.savefig('../FIGURES/triaxus_tow5l2_in2022_v06.jpg', bbox_inches = 'tight')
plt.show()

# Zonal diff

fig, axs = figure()
n = 0
for ax in axs:
    if var[n] != 'temperature':
        units = t1l2[var[n]].units
    else:
        units = '$^{\circ}$C'
    c = axs[n].contourf(t1l2['longitude'][:-1], -t1l2['pressure'],
                        -t1l2[var[n]].diff('time'), levels = 20, 
                        cmap = 'RdBu_r')
    plt.colorbar(c, ax = axs[n]).set_label(var[n]+' ('+units+')')
    n += 1
fig.suptitle('IN2022_V06 - TRIAXUS tow 1 - zonal diff', x = 0.5, y = 0.95, 
             fontsize = 20)
plt.savefig('../FIGURES/triaxus_tow1_in2022_v06_zonal_diff.jpg', bbox_inches = 'tight')

fig, axs = figure()
n = 0
for ax in axs:
    if var[n] != 'temperature':
        units = t2l2[var[n]].units
    else:
        units = '$^{\circ}$C'
    c = axs[n].contourf(t2l2['longitude'][:-1], -t2l2['pressure'],
                        -t2l2[var[n]].diff('time'), levels = 20, 
                        cmap = 'RdBu_r')
    plt.colorbar(c, ax = axs[n]).set_label(var[n]+' ('+units+')')
    n += 1
fig.suptitle('IN2022_V06 - TRIAXUS tow 2 - zonal diff', x = 0.5, y = 0.95, 
             fontsize = 20)
plt.savefig('../FIGURES/triaxus_tow2_in2022_v06_zonal_diff.jpg', bbox_inches = 'tight')
