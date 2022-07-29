# PLOT sst and ssh vectors
# Created by Amandine Schaeffer 8 May 2021

from matplotlib import pyplot as plt
import numpy as np
import datetime
#import pandas as pd # import libraries
import netCDF4 # import libraries
from datetime import date
from datetime import timedelta
import glob
import pandas as pd 
import xarray as xar
from netCDF4 import Dataset

import gsw
#OR from geopy import distance

# Housekeeping

base_path_bathy = './BATHY/'
plot_folder = '../FIGURES/'

# Path files Opendap

# SSH GEO
file_SSH = 'http://thredds.aodn.org.au/thredds/dodsC/IMOS/OceanCurrent/GSLA/NRT00/2021/IMOS_OceanCurrent_HV_20210509T180000Z_GSLA_FV02_NRT00_C-20210512T231907Z.nc.gz'
#file_SSH = 'http://thredds.aodn.org.au/thredds/dodsC/IMOS/OceanCurrent/GSLA/NRT00/2021/IMOS_OceanCurrent_HV_20210506T180000Z_GSLA_FV02_NRT00_C-20210509T231806Z.nc.gz'
#file_SSH = 'http://thredds.aodn.org.au/thredds/dodsC/IMOS/OceanCurrent/GSLA/NRT00/2021/IMOS_OceanCurrent_HV_20210504T180000Z_GSLA_FV02_NRT00_C-20210507T231818Z.nc.gz'
#file_SSH = 'http://thredds.aodn.org.au/thredds/dodsC/IMOS/OceanCurrent/GSLA/NRT00/2021/IMOS_OceanCurrent_HV_20210503T180000Z_GSLA_FV02_NRT00_C-20210506T231824Z.nc.gz'
#file_SSH = 'http://thredds.aodn.org.au/thredds/dodsC/IMOS/OceanCurrent/GSLA/NRT00/2021/IMOS_OceanCurrent_HV_20210502T180000Z_GSLA_FV02_NRT00_C-20210505T231805Z.nc.gz'
#file_SSH = '/home/z3340777/hdrive/My_documents/AUSTRALIE2/SATELLITE/DATA/SSH/IMOS_Yearfiles/MEAN_1992_2019_IMOS_OceanCurrent_HV_DM01.nc' # your file name with the eventual path

# SST
file_sst = 'http://thredds.aodn.org.au/thredds/dodsC/IMOS/SRS/SST/ghrsst/L3S-1d/ngt/2021/20210507152000-ABOM-L3S_GHRSST-SSTskin-AVHRR_D-1d_night.nc'
#file_sst = 'http://thredds.aodn.org.au/thredds/dodsC/IMOS/SRS/SST/ghrsst/L3S-1d/dn/2021/20210507092000-ABOM-L3S_GHRSST-SSTfnd-AVHRR_D-1d_dn.nc'
#file_sst = 'http://thredds.aodn.org.au/thredds/dodsC/IMOS/SRS/SST/ghrsst/L3S-1d/dn/2021/20210506092000-ABOM-L3S_GHRSST-SSTfnd-AVHRR_D-1d_dn.nc'
#file_sst = 'http://thredds.aodn.org.au/thredds/dodsC/IMOS/SRS/SST/ghrsst/L3S-1d/ngt/2021/20210506152000-ABOM-L3S_GHRSST-SSTskin-AVHRR_D-1d_night.nc'
#file_sst = 'http://thredds.aodn.org.au/thredds/dodsC/IMOS/SRS/SST/ghrsst/L3S-1d/ngt/2021/20210505152000-ABOM-L3S_GHRSST-SSTskin-AVHRR_D-1d_night.nc'
#file_sst = 'http://thredds.aodn.org.au/thredds/dodsC/IMOS/SRS/SST/ghrsst/L3S-1d/ngt/2021/20210505152000-ABOM-L3S_GHRSST-SSTskin-AVHRR_D-1d_night.nc'
#file_sst = 'http://thredds.aodn.org.au/thredds/dodsC/IMOS/SRS/SST/ghrsst/L3S-1d/dn/2021/20210505092000-ABOM-L3S_GHRSST-SSTfnd-AVHRR_D-1d_dn.nc'
#file_sst = 'http://thredds.aodn.org.au/thredds/dodsC/IMOS/SRS/SST/ghrsst/L3S-1d/dn/2021/20210504092000-ABOM-L3S_GHRSST-SSTfnd-AVHRR_D-1d_dn.nc'
#file_sst = 'http://thredds.aodn.org.au/thredds/dodsC/IMOS/SRS/SST/ghrsst/L4/RAMSSA/2021/20210506120000-ABOM-L4_GHRSST-SSTfnd-RAMSSA_09km-AUS-v02.0-fv01.0.nc'
#file_sst = 'http://thredds.aodn.org.au/thredds/dodsC/IMOS/SRS/SST/ghrsst/L3S-1d/ngt/2021/20210504152000-ABOM-L3S_GHRSST-SSTskin-AVHRR_D-1d_night.nc'
#file_sst = 'http://thredds.aodn.org.au/thredds/dodsC/IMOS/SRS/SST/ghrsst/L3S-1d/ngt/2021/20210503152000-ABOM-L3S_GHRSST-SSTskin-AVHRR_D-1d_night.nc'
#file_sst = 'https://thredds.aodn.org.au/thredds/dodsC/IMOS/SRS/SST/ghrsst/L3S-1d/ngt/2021/20210501152000-ABOM-L3S_GHRSST-SSTskin-AVHRR_D-1d_night.nc'
#file = 'imos-data/CSIRO/Climatology/SSTAARS/2017/AODN-product/SSTAARS_daily_fit_028.nc'
#file_sst = 'http://thredds.aodn.org.au/thredds/dodsC/IMOS/SRS/SST/ghrsst/L3S-1d/dn/2021/20210501092000-ABOM-L3S_GHRSST-SSTfnd-AVHRR_D-1d_dn.nc'

# CHL MODIS
file_chl = 'http://thredds.aodn.org.au/thredds/dodsC/IMOS/SRS/OC/gridded/aqua/P1D/2021/05/A.P1D.20210506T053000Z.aust.chl_oc3.nc'
#file_chl = 'http://thredds.aodn.org.au/thredds/dodsC/IMOS/SRS/OC/gridded/aqua/P1D/2021/05/A.P1D.20210505T053000Z.aust.chl_oc3.nc'
#file_chl = 'http://thredds.aodn.org.au/thredds/dodsC/IMOS/SRS/OC/gridded/aqua/P1D/2021/05/A.P1D.20210504T053000Z.aust.chl_oc3.nc'
#file_chl = 'http://thredds.aodn.org.au/thredds/dodsC/IMOS/SRS/OC/gridded/aqua/P1D/2021/05/A.P1D.20210502T053000Z.aust.chl_oc3.nc'
#file_chl = 'http://thredds.aodn.org.au/thredds/dodsC/IMOS/SRS/OC/gridded/aqua/P1D/2021/05/A.P1D.20210504T053000Z.aust.chl_oc3.nc'


# Bathymetry

# Load coastine
coast=np.loadtxt(base_path_bathy + "eaccoast.dat")
File_bathy = Dataset(base_path_bathy + 'bathy_dbdb2_v30_AUSTRALIA.nc')
latBathy = File_bathy.variables['lat'][:]        # units and all
lonBathy = File_bathy.variables['lon'][:]        # units and all
heightBathy = File_bathy.variables['height'][:]        # units and all

# Get SSH from IMOS

ds_SSH = xar.open_dataset(file_SSH)
ds_SSH

LAT_SSH = ds_SSH['LATITUDE'].values
LON_SSH = ds_SSH['LONGITUDE'].values
#read velocity
UCUR = ds_SSH['UCUR'][-1,:,:].values
VCUR = ds_SSH['VCUR'][-1,:,:].values
SPEED = np.sqrt(UCUR**2+VCUR**2)

GSLA = ds_SSH['GSLA'][-1,:,:].values
ds_SSH.TIME.values[-1]

# PLOT Vectors

#plt.pcolor(LON,LAT,VCUR)
xlim1=155 ; xlim0=150; ylim0=-38;ylim1=-29
#xlim1=158 ; xlim0=148; ylim0=-38;ylim1=-30
clim0 = (0, 1.5)
import cmocean
cmap = cmocean.cm.speed;

fig_x = 10
fig_y = fig_x * np.cos(abs(ylim1-ylim0)/180*np.pi)
figure1 = plt.figure(figsize=(fig_x,fig_y), dpi=80, facecolor='w', edgecolor='k')
ax = plt.subplot(1,1,1)
Q_bathy = plt.contour(lonBathy, latBathy, heightBathy, levels = [-2000, -1000, -200], colors = 'grey', linestyles='solid')

Q1 = plt.pcolor(LON_SSH,LAT_SSH,SPEED)
plt.clim(clim0)
Q = plt.quiver(LON_SSH,LAT_SSH,UCUR,VCUR, units='width',scale = 10, width=0.004, headwidth =3)
Q2 = plt.contour(LON_SSH,LAT_SSH,GSLA, levels = 20, colors = 'black', linestyles='dashed')
    #skip = (slice(None,None,2),slice(None,None,2))
plt.plot(coast[:,0],coast[:,1])
#     #plt.plot(lon[ind_lon],lat[ind_lat],'kX',markersize=12)
#     #plt.title(dates[0].strftime('%d/%m/%Y') + ' - ' + dates[-1].strftime('%d/%m/%Y'),fontsize=15)                  

# LABELS / LIMITS
plt.xlabel('Longitude [$^o$E]',fontsize=16)                  
plt.ylabel('Latitude [$^o$N]',fontsize=16)                  
plt.ylim(ylim0,ylim1)     # plt.ylim(-30.6, -29.8)     
plt.xlim(xlim0,xlim1) # plt.xlim(152.8, 154)
cb = plt.colorbar(Q1)
cb.set_label('Speed [m s$^-1$]', size=14)
plt.title('Geostrophy SSH (' + str(ds_SSH.TIME.values[-1])[0:16] +')',fontsize=18)    

#plt.savefig( plot_folder +'plot_SVPs_ALL_20210506.png', bbox_inches='tight', pad_inches=0.5, dpi=150)
plt.savefig(plot_folder + 'plot_SSH_'+ str(ds_SSH.TIME.values[-1])[0:13] + '.png')
plt.show()

#ymin,ymax = ax.get_ylim()
##aspect_ratio = np.cos((ymax + ymin)/2/180*np.pi)  # set(gca,'dataaspectratio',[1 cos((Y(2)+Y(1))/2/180*pi) 1]) 
##adjustFigAspect(figure1,aspect=aspect_ratio)
#     cb.set_label('Speed [m s$^{-1}]$',fontsize=18)
#     ax.quiverkey(Q, 0.1, 0.1, 1, '1' + 'm s$^{-1}$', labelpos='E', coordinates='axes')
#     plt.plot([153.397], [-30.268], 'ko') # CH100
# #    plt.plot(VCUR_min_lon,lat,'xk')
#     plt.rcParams.update({'font.size': 12})
#     plt.savefig('../PLOTS2/'+ str(file_name[25:29]) + '_yearly_' + str(ds_yearly_OK.TIME[ind].values)[0:4] + '.png')
#     plt.show()


# In[ ]:


print(plot_folder + 'plot_SSH_'+ str(ds_SSH.TIME.values[-1])[0:13] + '.png')


# Get SST

# In[ ]:


ds = xar.open_dataset(file_sst+'#fillmismatch')
#ds = xar.open_dataset(file_sst)
ds

# SLICE so it is smalle
# BIG ds2 = ds.sel(lat=slice(-25, -40), lon=slice(148, 160))
ds_SST = ds.sel(lat=slice(-29, -40), lon=slice(149, 158))
#ds_SST = ds  # FOR RAMSSA

LAT = ds_SST['lat'].values
LON = ds_SST['lon'].values
SST_K = ds_SST['sea_surface_temperature'][-1,:,:].values
#SST_K = ds_SST['analysed_sst'][0,:,:].values   # FOR RAMSSA

SST = SST_K - 273.15  # Celcius
TIME = ds_SST['time'].values
ds_SST.time.values[-1]


# In[ ]:


# PLOT
# LIMITS for plot
xlim1=155 ; xlim0=150; ylim0=-38;ylim1=-29
#xlim1=158 ; xlim0=148; ylim0=-38;ylim1=-30
clim0 = (15, 24)
import cmocean
cmap = cmocean.cm.thermal;

# FIGURE
fig_x = 12
fig_y = fig_x * np.cos(abs(ylim1-ylim0)/180*np.pi)
figure1 = plt.figure(figsize=(fig_x,fig_y), dpi=80, facecolor='w', edgecolor='k')
ax = plt.subplot(1,1,1)

# BATHY contours
Q_bathy = plt.contour(lonBathy, latBathy, heightBathy, levels = [-5000, -2000, -1000, -200], colors = 'grey', linestyles='solid')

# PCOLOR SST
#Q = plt.pcolor(ds_SST.lon[4000:4400], ds_SST.lat[2000:3000],SST[2000:3000,4000:4400])
#Q = plt.pcolor(ds_SST.lon[slice(None,None,2)], ds_SST.lat[slice(None,None,2)],SST[slice(None,None,2),slice(None,None,2)])
Q = plt.pcolor(ds_SST.lon, ds_SST.lat,SST, cmap=cmap)
plt.clim(clim0)
Q2 = plt.quiver(LON_SSH,LAT_SSH,UCUR,VCUR, units='width',scale = 10, width=0.004, headwidth =3)

plt.plot(coast[:,0],coast[:,1],'k')

# SYMBOLS
plt.plot(151.33,-34.3666,'bo',markersize=10)

# LABELS / LIMITS
plt.xlabel('Longitude [$^o$E]',fontsize=16)                  
plt.ylabel('Latitude [$^o$N]',fontsize=16)                  
plt.ylim(ylim0,ylim1)     # plt.ylim(-30.6, -29.8)     
plt.xlim(xlim0,xlim1) # plt.xlim(152.8, 154)
cb = plt.colorbar()
cb.set_label('SST [$^o$C]', size=14)
plt.clim(clim0)
plt.title('SST (' + str(ds_SST.time.values[-1])[0:16] + ') / SSH (' + str(ds_SSH.TIME.values[-1])[0:16] +')',fontsize=18)    

# RATIO
ymin,ymax = ax.get_ylim()


# SAVE
plt.savefig(plot_folder + 'plot_SST_'+ str(ds_SST.time.values[-1])[0:13] + '.png', bbox_inches='tight', pad_inches=0.5, dpi=150)
plt.show()


# Get CHL

# In[ ]:



ds = xar.open_dataset(file_chl+'#fillmismatch')
#ds = xar.open_dataset(file_sst)

# SLICE so it is smalle
# BIG ds2 = ds.sel(lat=slice(-25, -40), lon=slice(148, 160))
ds_CHL = ds.sel(latitude=slice(-29, -40), longitude=slice(149, 158))


LAT_chl = ds_CHL['latitude'].values
LON_chl = ds_CHL['longitude'].values
CHL = ds_CHL['chl_oc3'][-1,:,:].values
TIME = ds_CHL['time'].values


# In[ ]:


# PLOT
# LIMTITS for plot
xlim1=155 ; xlim0=150; ylim0=-38;ylim1=-29
#xlim1=158 ; xlim0=148; ylim0=-38;ylim1=-30
clim0 = (0, 4)
import cmocean
cmap = cmocean.cm.algae ;

# FIGURE
fig_x = 10
fig_y = fig_x * np.cos(abs(ylim1-ylim0)/180*np.pi)
figure1 = plt.figure(figsize=(fig_x,fig_y), dpi=80, facecolor='w', edgecolor='k')
ax = plt.subplot(1,1,1)

# BATHY contours
Q_bathy = plt.contour(lonBathy, latBathy, heightBathy, levels = [-5000, -2000, -1000, -200], colors = 'grey', linestyles='solid')

# PCOLOR SST
#Q = plt.pcolor(ds_CHL.lon[4000:4400], ds_CHL.lat[2000:3000],SST[2000:3000,4000:4400])
skip = (slice(None,None,5),slice(None,None,5))
Q = plt.pcolor(LON_chl[slice(None,None,5)], LAT_chl[slice(None,None,5)],CHL[skip])
plt.clim(clim0)
Q2 = plt.quiver(LON_SSH,LAT_SSH,UCUR,VCUR, units='width',scale = 10, width=0.004, headwidth =3)

plt.plot(coast[:,0],coast[:,1],'k')

# SYMBOLS
plt.plot(151.33,-34.3666,'bo',markersize=10)

# LABELS / LIMITS
plt.xlabel('Longitude [$^o$E]',fontsize=16)                  
plt.ylabel('Latitude [$^o$N]',fontsize=16)                  
plt.ylim(ylim0,ylim1)     # plt.ylim(-30.6, -29.8)     
plt.xlim(xlim0,xlim1) # plt.xlim(152.8, 154)
cb = plt.colorbar()
cb.set_label('CHL oc3 [mg m$^{-3}$]', size=14)
plt.clim(clim0)
plt.title('CHL (' + str(ds_CHL.time.values[-1])[0:16] + ') / SSH (' + str(ds_SSH.TIME.values[-1])[0:16] +')',fontsize=18)    

# RATIO
ymin,ymax = ax.get_ylim()


# SAVE
plt.savefig(plot_folder + 'plot_CHL_'+ str(ds_CHL.time.values[-1])[0:13] + '.png', bbox_inches='tight', pad_inches=0.5, dpi=150)
plt.show()


# In[ ]:


#NOTE TO CALCULATE DISTANCE
# Install TEOS
#!pip install git+https://github.com/TEOS-10/GSW-Python.gitw
#Help https://github.com/TEOS-10/python-gsw/blob/master/gsw/gibbs/earth.py
#!conda install -c conda-forge gsw
#pip install gsw
# import gsw
#lon = [159, 220]
#lat = [-35, 35]
#gsw.distance(lon, lat)

#!pip install geopy
from geopy import distance
#help(geopy.distance)

coords_1 = (52.2296756, 21.0122287)
coords_2 = (52.406374, 16.9251681)
dist = distance.distance(coords_1, coords_2).nm
dist


# In[ ]:




