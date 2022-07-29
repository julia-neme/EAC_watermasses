#!/usr/bin/env python
# coding: utf-8

# # Plot underway
# Created by Amandine Schaeffer 8 May 2021

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd 
import xarray as xar
import cmocean
import cartopy.crs as ccrs
import datetime
import glob


# Paths
path_bathy = './BATHY/'
path_fig = '../FIGURES/'

# Bathymetry and coastine
coast=np.loadtxt(path_bathy + "eaccoast.dat")
bathy = xar.open_dataset('BATHY/bathy_dbdb2_v30_AUSTRALIA.nc')
lon_bathy = bathy['lon']
lat_bathy = bathy['lat']
depth = bathy['height']

# set figure parameters and other variable arrays
xlim0=153; xlim1=155.5 ; ylim0=-27.8;  ylim1=-26.4 # Middle moorings
xticks = [152,152.5,153,153.5,154,154.5,155,155.5,156]
#xticklabels = [152,'',153,'',154,'',155,'',156]
yticks = [-30,-29.5,-29,-28.5,-28,-27.5,-27,-26.5,-26,-25.5,-25]
brislon, brislat = 153, -27.5
mooringlons = [153.90, 154.002, 154.134, 154.286, 154.639, 155.30230]
mooringlats = [-27.325, -27.316, -27.281, -27.244, -27.202, -27.1026]

# Set date to plot
DATE_START = '2022-07-14'; DATE_END = '2022-07-22';  
DATE_END = str(datetime.datetime.now())[0:10];
#str(datetime.datetime.now())[0:10]

# Get data Underway
ds_underw_all = xar.open_dataset(r'..\..\underway\in2022_v06uwy.nc')
# Create variable time 
epoch = pd.Timestamp(ds_underw_all.Epoch[-22:-3])
start_time = epoch + pd.to_timedelta(int(ds_underw_all.Epoch[0:8]),unit='S')
uwy_date_time       = start_time + pd.to_timedelta(5.0*ds_underw_all['sample'].values,unit='S')
#ds_underw_all['time'] = uwy_date_time
ds_underw_all    = ds_underw_all.reindex(sample=uwy_date_time)
ds_underw_all    = ds_underw_all.rename({'sample':'time'})
# Subset 
ds_underw = ds_underw_all.sel(time=slice(DATE_START, DATE_END))
TIME = ds_underw['time'].values
# Loop to fill in grid with densities
import gsw
DENS = np.zeros(len(ds_underw.rawTsgSensorTemp.values))
for i in range(len(ds_underw.rawTsgSensorTemp.values)):
    DENS[i]=gsw.rho(ds_underw.salinity.values[i],ds_underw.rawTsgSensorTemp.values[i],0)

ds_underw_all.close()
ds_underw.close()



# Get data ADCP
## Deep ADCP
#if ds_adcp_all:    
#    del ds_adcp_all
#    del ds_adcp75
#    del ds_adcp150
ds_adcp_all = xar.open_dataset(r'..\..\adcp\uhdas\proc\os75nb\contour\os75nb.nc')
ds_adcp75 = ds_adcp_all.sel(time=slice(DATE_START, DATE_END))
ds_adcp_all.close()
# Check: DEPTH_adcp75[:,10] # 200m depth

## Shallow ADCP
ds_adcp_all = xar.open_dataset(r'..\..\adcp\uhdas\proc\os150nb\contour\os150nb.nc')
ds_adcp150 = ds_adcp_all.sel(time=slice(DATE_START, DATE_END))
ds_adcp_all.close()
# Check: DEPTH_adcp150[:,4] # 50m depth


## Plots
# Subplots temp and sal and O2 and FLuo
vminT=18; vmaxT=23
vminS=35.3; vmaxS=35.6
vminO=205; vmaxO=220
vminF=0; vmaxF=3

def figure_map():
    # plot coastline
    ax.plot(coast[:,0],coast[:,1], c='dimgrey')
    # shade land light grey
    ax.contourf(lon_bathy, lat_bathy, depth.where(depth>0),colors = 'gainsboro', linestyles = 'solid',transform=ccrs.PlateCarree())
    # label bathymetry contours
    ax.contour(lon_bathy, lat_bathy[800:1000], depth[800:1000,:], levels = [-2000, -1000, -200], colors = 'lightgrey', linestyles='-')
    # plot location of Brisbane
    ax.scatter(brislon, brislat, marker='o', c='dimgrey', s=80, zorder=2,transform=ccrs.PlateCarree())
    ax.text(153, -27.7, 'Brisbane', c='dimgrey', fontsize=15)
    # scatter plot mooring locations
    ax.scatter(mooringlons, mooringlats, marker='o', c='k', s=80, zorder=2,transform=ccrs.PlateCarree())
    # set xlimits and labels
    ax.set_xlim([xlim0, xlim1])
    ax.set_ylim([ylim0, ylim1])
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xticklabels(xticks, fontsize=16)
    ax.set_yticklabels(yticks, fontsize=16)
    ax.set_xlabel('Longitude ($^o$E)',fontsize=18)                  
    ax.set_ylabel('Latitude ($^o$N)',fontsize=18) 
    ax.grid()
    return ax


fig, axes = plt.subplots(2,2, figsize=(25, 20), subplot_kw=dict(projection=ccrs.PlateCarree()))
scale = '10m'

#Subplot 1
ax = axes[0,0]
ax = figure_map()
ax.set_title("Underway temperature" + DATE_START + '_' + DATE_END, fontsize=24)
Q1 = ax.scatter(ds_underw.longitude,ds_underw.latitude, c=ds_underw.rawTsgSensorTemp, s=80, cmap=cmocean.cm.thermal, vmin=vminT, vmax=vmaxT, label='underrway',transform=ccrs.PlateCarree())
#ax.quiver(ds_adcp75.lon,ds_adcp75.lat,ds_adcp75.u[:,1],ds_adcp75.v[:,1],color='k', units='width',scale = 15, label=str(int(np.nanmean(ds_adcp75.depth[:,1]))) + ' m',transform=ccrs.PlateCarree())
ax.quiver(ds_adcp150.lon,ds_adcp150.lat,ds_adcp150.u[:,1],ds_adcp150.v[:,1], color='b',units='width',scale = 10, label=str(int(np.nanmean(ds_adcp150.depth[:,1]))) + ' m',transform=ccrs.PlateCarree())
ax.quiver(ds_adcp150.lon,ds_adcp150.lat,ds_adcp150.u[:,4],ds_adcp150.v[:,4], color='k',units='width',scale = 10, label=str(int(np.nanmean(ds_adcp150.depth[:,4]))) + ' m',transform=ccrs.PlateCarree())
#ax.quiver(LON_adcp150,LAT_adcp150,U_adcp150[:,4],V_adcp150[:,4], color='k',units='width',scale = 10, label=str(int(np.nanmean(DEPTH_adcp150[:,3]))) + ' m',transform=ccrs.PlateCarree())
#ax.quiver(LON_adcp150,LAT_adcp150,U_adcp150[:,4],V_adcp150[:,4], color='r',units='width',scale = 10,width=0.004, headwidth =3, label='50 m ',transform=ccrs.PlateCarree())
ax.set_extent([xlim0, xlim1, ylim0, ylim1])
plt.colorbar(Q1, ax=ax, orientation="horizontal", pad=0.05).set_label("$^o$C",fontsize=22)
ax.legend(fontsize=22)

#Subplot 2
ax = axes[0,1]
ax = figure_map()
ax.set_title("Underway salinity", fontsize=24)
Q = ax.scatter(ds_underw.longitude,ds_underw.latitude, c=ds_underw.salinity, s=50, cmap=cmocean.cm.haline, vmin=vminS, vmax=vmaxS, label='underrway',transform=ccrs.PlateCarree())
ax.set_extent([xlim0, xlim1, ylim0, ylim1])
plt.colorbar(Q, ax=ax, orientation="horizontal", pad=0.05).set_label("PSU",fontsize=22)
#
#Subplot 3
ax = axes[1,0]
ax = figure_map()
ax.set_title("Underway DO", fontsize=24)
Q = ax.scatter(ds_underw.longitude,ds_underw.latitude, c=ds_underw.do, s=50, cmap=cmocean.cm.oxy, vmin=vminO, vmax=vmaxO, label='underrway',transform=ccrs.PlateCarree())
ax.set_extent([xlim0, xlim1, ylim0, ylim1])
plt.colorbar(Q, ax=ax, orientation="horizontal", pad=0.05).set_label("uM",fontsize=22)

#Subplot 4
ax = axes[1,1]
ax = figure_map()
ax.set_title("Underway fluorescence", fontsize=24)
Q1 = ax.scatter(ds_underw.longitude,ds_underw.latitude, c=ds_underw.fluorescenceConcentration, s=80, cmap='rainbow', vmin=vminF, vmax=vmaxF, label='underrway',transform=ccrs.PlateCarree())
#cmap=cmocean.cm.algae
ax.set_extent([xlim0, xlim1, ylim0, ylim1])
plt.colorbar(Q1, ax=ax, orientation="horizontal", pad=0.05).set_label("ug/L", fontsize=22)


# SAVE
plt.tight_layout()
plt.savefig(path_fig + 'plot_underway_from_day_'+ DATE_START + '_' + DATE_END + '.png', bbox_inches='tight', pad_inches=0.5, dpi=300)
plt.show()

ds_adcp75.close()
ds_adcp150.close()



######### PLOT all track
## Plots
##
#def figure_map():
#    # plot coastline
#    ax.plot(coast[:,0],coast[:,1], c='dimgrey')
#    # shade land light grey
#    ax.contourf(lon_bathy, lat_bathy, depth.where(depth>0),colors = 'gainsboro', linestyles = 'solid',transform=ccrs.PlateCarree())
#    # label bathymetry contours
#    ax.contour(lon_bathy, lat_bathy[800:1000], depth[800:1000,:], levels = [-2000, -1000, -200], colors = 'lightgrey', linestyles='-')
#    # plot location of Brisbane
#    ax.scatter(brislon, brislat, marker='o', c='dimgrey', s=80, zorder=2,transform=ccrs.PlateCarree())
#    ax.text(123, -27.7, 'Brisbane', c='dimgrey', fontsize=22)
#    # scatter plot mooring locations
#    Q0 = ax.scatter(mooringlons, mooringlats, marker='o', c='k', s=80, zorder=2,transform=ccrs.PlateCarree())
#    # set xlimits and labels
#    ax.set_xlim([xlim0, xlim1])
#    ax.set_ylim([ylim0, ylim1])
#    ax.set_xticks(xticks)
#    ax.set_yticks(yticks)
#    ax.set_xticklabels(xticks, fontsize=16)
#    ax.set_yticklabels(yticks, fontsize=16)
#    ax.set_xlabel('Longitude ($^o$E)',fontsize=18)                  
#    ax.set_ylabel('Latitude ($^o$N)',fontsize=18) 
#    ax.grid()
#    return ax


# Function to read Rapid cast
import os
import seawater.eos80 as sw # pip install seawater
def fun_read(file):
    with open(file,'r') as rapidctdfile:
        line = rapidctdfile.readline()
        print(line)
        while "[DATA]" not in line:
            line = rapidctdfile.readline()
            if "Latitude" in line:                
                    print(line)
                    latitude = float(line.split("=",1)[1])
            if "Longitude" in line:                
                    print(line)
                    longitude = float(line.split("=",1)[1])
        #read the data headers
        header = rapidctdfile.readline().split()
        #read the data units
        units = rapidctdfile.readline().split()
        print(header)
        print(units)
        dataline = rapidctdfile.readline()
        pressure = []
        temp = []
        conductivity = []
        rapid_file_salinity = []
        
        #expected data format
        #2019/09/15 03:46:16.358	2.869	2.885	21.322	49.694	35.340	1525.495	1024.677	731486       
        rapidctd_dataset = []
        while dataline:
            datastrings = dataline.split()
            
            #print(datastrings)
            numbers = []
            for i in range(2,len(datastrings)):
                #print(dataline[i])
                numbers.append(float(datastrings[i]))
            pressure.append(numbers[1])
            temp.append(numbers[2])
            conductivity.append(numbers[3])
            rapid_file_salinity.append(numbers[4])
           # dataset.append(numbers)
            dataline = rapidctdfile.readline()            
       # derive salinity     
        ref_cond = 42.914
        rt = [x/ref_cond for x in conductivity]    
        salinity = sw.salt(rt,temp,pressure)    
        # plt.plot(salinity,pressure)
        return pressure, temp, salinity, datastrings, header, latitude, longitude


fig, ax = plt.subplots(1,1, figsize=(15, 20), subplot_kw=dict(projection=ccrs.PlateCarree()))
scale = '10m'
# Underway
ax = figure_map()
Q0 = ax.scatter(mooringlons, mooringlats, marker='o', c='k', s=80, zorder=2,transform=ccrs.PlateCarree())
ax.set_title("IN2022_v06", fontsize=24)
Q1 = ax.scatter(ds_underw.longitude,ds_underw.latitude, c=ds_underw.rawTsgSensorTemp, s=80, cmap=cmocean.cm.thermal, vmin=18, vmax=23, label='underrway',transform=ccrs.PlateCarree())

# Add CTDs
#fig, ax = plt.subplots(1,1, figsize=(15, 20), subplot_kw=dict(projection=ccrs.PlateCarree()))
#ax = figure_map()
list_FILES = glob.glob(r"../../ctd/processing/in2022_v06/cap/cappro/avg/*.nc")
N_FILES = len(list_FILES)
for f in range(N_FILES):
    file = list_FILES[f]
    data = xar.open_dataset(file)
    lon = np.round(data['longitude'].squeeze().values, 2)
    lat = np.round(data['latitude'].squeeze().values, 2)
    tim = str(data['time'].values[0])[:-10]
    Q2 = ax.scatter(lon,lat,c ='cyan',s=70,label='CTD')
    
# Add Rapidcast
list_FILES = sorted(glob.glob(r'..\DATA\RAPIDCAST\VL_*vp2'))
N_FILES = len(list_FILES)
for f in range(N_FILES):
    file = list_FILES[f]
    [pressure, temp, salinity, datastrings, header, latitude, longitude] = fun_read(file)
    Q3 = ax.scatter(longitude,latitude,c ='blue',s=70,label='Rapid Cast',linewidth=3)
    
    
# Add SBP line
Q4, = ax.plot([154.050859, 154.068007],[-27.556572, -28.371804],'grey',transform=ccrs.PlateCarree(),label='SBP / XBT line');

# Add Triaxus
t1l2 = xar.open_dataset(r'../../triaxus/processing/in2022_v06/cap/cappro/avg/section/tow1Leg2AvgCast.nc')
t2l2 = xar.open_dataset(r'../../triaxus/processing/in2022_v06/cap/cappro/avg/section/tow2Leg2AvgCast.nc')
t3l2 = xar.open_dataset(r'../../triaxus/processing/in2022_v06/cap/cappro/avg/section/tow3Leg2AvgCast.nc')
t4l2 = xar.open_dataset(r'../../triaxus/processing/in2022_v06/cap/cappro/avg/section/tow4Leg2AvgCast.nc')
t4l4 = xar.open_dataset(r'../../triaxus/processing/in2022_v06/cap/cappro/avg/section/tow4Leg4AvgCast.nc')
t5l2 = xar.open_dataset(r'../../triaxus/processing/in2022_v06/cap/cappro/avg/section/tow5Leg2AvgCast.nc')
Q6, = ax.plot(t1l2.longitude, t1l2.latitude,'-k',transform=ccrs.PlateCarree(),label='Triaxus lines',linewidth=3);
Q5 = ax.plot(t2l2.longitude, t2l2.latitude,'-k',transform=ccrs.PlateCarree(),linewidth=3);
Q5 = ax.plot(t3l2.longitude, t3l2.latitude,'-k',transform=ccrs.PlateCarree(),linewidth=3);
Q5 = ax.plot(t4l2.longitude, t4l2.latitude,'-k',transform=ccrs.PlateCarree(),linewidth=3);
Q5 = ax.plot(t4l4.longitude, t4l4.latitude,'-k',transform=ccrs.PlateCarree(),linewidth=3);
Q5 = ax.plot(t5l2.longitude, t5l2.latitude,'-k',transform=ccrs.PlateCarree(),linewidth=3);

plt.legend([Q0,Q2, Q3,Q4, Q6],['EAC moorings', 'CTDs','RapidCasts', 'SBP / XBT line','Triaxus lines'],fontsize=22)
plt.colorbar(Q1, ax=ax, orientation="horizontal", pad=0.05).set_label("Underway temperature [$^o$C]",fontsize=22)
ax.set_extent([152.2, 156, -29, -26])
plt.savefig(path_fig + 'plot_ALLvoyage_from_day_'+ DATE_START + '_' + DATE_END + '.png', bbox_inches='tight', pad_inches=0.5, dpi=300)
plt.show()    

    
    
    
    
    
    
    
    
    
    
    
    