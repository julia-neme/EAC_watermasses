#################################################################
# PLOT SATELLITE DATA (SST, CHL, or SSH)                        #
# Created by Hannah Dawson 15 July 2022                         #
# run from command line with python3 plot_IMOS_satellite.py var #
#################################################################


from matplotlib import pyplot as plt
import numpy as np
import xarray as xr
from glob import glob
import scipy.io
import cmocean as cmo
#import sys
import re
import os
# Read in variable passed to python script from command line
# must be one of 'chl', 'sst', or 'ssh'. 
var = 'SST_mat'
print(var)
#var = 'chl'

# set paths (assumes you're running from SCRIPTS folder)
path_bathy = './BATHY/'
path_fig = '../FIGURES/'
path_data = '../DATA/SATELLITE/'



# Read in bathymetry and coast file
coast=np.loadtxt(path_bathy + "eaccoast.dat")
bathy = xr.open_dataset('BATHY/bathy_dbdb2_v30_AUSTRALIA.nc')
lon_bathy = bathy['lon']
lat_bathy = bathy['lat']
depth = bathy['height']

# set figure parameters and other variable arrays
xlim1=156; xlim0=152.5; ylim0=-29;ylim1=-25
xticks = [152,152.5,153,153.5,154,154.5,155,155.5,156]
xticklabels = [152,'',153,'',154,'',155,'',156]
yticks = [-30,-29,-28,-27,-26,-25,-24]
fig_x = 10
fig_y = fig_x * np.cos(abs(ylim1-ylim0)/180*np.pi)
brislon, brislat = 153, -27.5
mooringlons = [153.90, 154.002, 154.134, 154.286, 154.639, 155.30230]
mooringlats = [-27.325, -27.316, -27.281, -27.244, -27.202, -27.1026]

# define figure plot
def figure_sat():
    fig = plt.figure(figsize=(fig_x,fig_y), dpi=80, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(1,1,1)
    # plot coastline
    ax.plot(coast[:,0],coast[:,1], c='dimgrey')
    # shade land light grey
    ax.contourf(lon_bathy, lat_bathy, depth.where(depth>0),colors = 'gainsboro', linestyles = 'solid')
    # label bathymetry contours
    cs = ax.contour(lon_bathy, lat_bathy[800:1000], depth[800:1000,:], levels = [-2000, -1000, -200], colors = 'lightgrey', linestyles='-')
    # plot location of Brisbane
    ax.scatter(brislon, brislat, marker='o', c='dimgrey', s=80, zorder=2)
    ax.text(152.5, -27.7, 'Brisbane', c='dimgrey', fontsize=15)
    # scatter plot mooring locations
    ax.scatter(mooringlons, mooringlats, marker='o', c='k', s=80, zorder=2)
    # set xlimits and labels
    ax.set_xlim([xlim0, xlim1])
    ax.set_ylim([ylim0, ylim1])
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xticklabels(xticklabels, fontsize=14)
    ax.set_yticklabels(yticks, fontsize=14)
    ax.set_xlabel('Longitude ($^o$E)',fontsize=16)                  
    ax.set_ylabel('Latitude ($^o$N)',fontsize=16)
    
    return ax



# Define files containing variable
#files = sorted(glob(path_data+f'/{var.upper()}/*.nc'))
files = sorted(glob(path_data + '/SST/Brisbane2_4hr/*.mat'))
files_SSH = sorted(glob(path_data + '/SSH/IMOS_OceanCurrent_HV_*.nc'))


### CHANGE TO GET ALL
for f, file in enumerate(files[-3:]):
    # load variable data and extract time for figure title
    
    if var == 'ssh':
        ds = xr.open_dataset(file)
        date_str = str(ds.TIME.values)[2:12]
    elif var == 'SST_mat':
        mat = scipy.io.loadmat(file)
        sst = mat['sst1']
        lons = mat['lns'][1,:]
        lats = mat['lts'][:,0]
        date_str = str(file[-14:-4])
        yy, mm, dd, hh = int(date_str[0:4]), int(date_str[4:6]), int(date_str[6:8]), int(date_str[8:])
    else:
        ds = xr.open_dataset(file)
        date_str = str(ds.time.values)[2:12]
        yy, mm, dd = int(date_str[0:4]), int(date_str[5:7]), int(date_str[8:])
    print(date_str)

    # Read in corresponding ssh file
    if (var =='sst' or var == 'chl' or var =='SST_mat'):
        #sshfile = sorted(glob(path_data+f'SSH/IMOS_OceanCurrent_HV_{yy:04d}{mm:02d}{dd:02d}*.nc'))[0]
        pattern = re.compile('IMOS_OceanCurrent_HV_'+ date_str[0:8]) #+ str(dd))
        dir = path_data + '\SSH/'
        same_date_ssh = 0
        for filepath in os.listdir(dir):
            if pattern.match(filepath):
                print('ok') 
                same_date_ssh = 1
                sshfile = sorted(glob(path_data+f'SSH/IMOS_OceanCurrent_HV_{yy:04d}{mm:02d}{dd:02d}*.nc')) # get the exact day
                print(sshfile[0])
                sshds = xr.open_dataset(sshfile[0])
        if same_date_ssh == 0:
            sshfile = sorted(glob(path_data+f'SSH/IMOS_OceanCurrent_HV_{yy:04d}{mm:02d}*.nc')) # get the last one (not ideal)
            sshds = xr.open_dataset(sshfile[-1])
        date_ssh = str(sshds.TIME.values)[2:12]
        
    # call figure plot
    ax = figure_sat()
    
    if var == 'sst':
        # plot variable
        p1 = ax.pcolormesh(ds.lon, ds.lat, ds.sea_surface_temperature.squeeze()-273.15, cmap=cmo.cm.thermal, vmin=20, vmax=25)
        ax.quiver(sshds.LONGITUDE, sshds.LATITUDE,sshds.UCUR.squeeze(),sshds.VCUR.squeeze(), units='width',scale = 10, width=0.004, headwidth =3)
        # plot colorbar
        cb = plt.colorbar(p1, extend='both', shrink=0.8)
        cb.set_label('SST ($^o$C)', size=14)
        plt.title(f'SST:  {date_str}',fontsize=18) 
        
    elif var == 'SST_mat':
    #    sst = mat['sst1']
    #lons = mat['lns'][1,:]
    #lats = mat['lts'][:,0]
#        p1 = ax.pcolormesh(lons, lats, sst.squeeze(), cmap=cmo.cm.thermal, vmin=19, vmax=23)
        p1 = ax.contourf(lons, lats, sst.squeeze(), cmap=cmo.cm.thermal,
                levels = np.arange(19.5, 23, .05), extend = 'both')
        ax.quiver(sshds.LONGITUDE, sshds.LATITUDE,sshds.UCUR.squeeze(),sshds.VCUR.squeeze(), units='width',scale = 10, width=0.004, headwidth =3)        
        cb = plt.colorbar(p1, extend='both', shrink=0.8)
        cb.set_label('SST ($^o$C)', size=14)
        plt.title(f'SST:  {date_str}, SSH: {date_ssh}',fontsize=18)         

    elif var == 'chl':
        #plot variable
        p1 = ax.pcolormesh(ds.longitude, ds.latitude, ds.chl_oc3.squeeze(), cmap=cmo.cm.haline, vmin=0, vmax=100)
        ax.quiver(sshds.LONGITUDE, sshds.LATITUDE,sshds.UCUR.squeeze(),sshds.VCUR.squeeze(), units='width',scale = 10, width=0.004, headwidth =3)
        # plot colorbar
        cb = plt.colorbar(p1, extend='max', shrink=0.8)
        cb.set_label('CHL (mg/m^{3})', size=14)
        plt.title(f'CHL:  {date_str}',fontsize=18) 
    else:
        # plot variable
        p1 = ax.pcolormesh(ds.LONGITUDE, ds.LATITUDE, ds.GSLA.squeeze(), cmap=cmo.cm.delta, vmin=-0.4, vmax=0.4, zorder=1)
        ax.quiver(ds.LONGITUDE, ds.LATITUDE,ds.UCUR.squeeze(),ds.VCUR.squeeze(), units='width',scale = 10, width=0.004, headwidth =3)
        # plot colorbar
        cb = plt.colorbar(p1, extend='both', shrink=0.8)
        cb.set_label('SLA (m)', size=14)
        plt.title(f'SLA:  {date_str}',fontsize=18)

    outfile = path_fig + f'satellite_{var}_{date_str}.png'
    print(outfile)
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.show()
    plt.close()