###############################################################################
# Making plots of ctd profile with SST SSH map of location. Plots will help 
# identify ctd's from the EAC, shelf, and eddies to make TS diagrams
# run as:
# python3 'voyage_year' 'voyage_code'
# Also, remember to create a directory with voyage name in EAS_class
# -----------------------------------
# by Julia Neme (in2022_v06)
###############################################################################

import cmocean as cm
import glob
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import xarray as xr
import warnings
warnings.filterwarnings("ignore")
from siphon.catalog import TDSCatalog
plt.rcParams['font.size'] = 14

coast = np.loadtxt("BATHY/eaccoast.dat")
bathy = xr.open_dataset('BATHY/bathy_dbdb2_v30_AUSTRALIA.nc')
lon_bathy = bathy['lon']
lat_bathy = bathy['lat']
depth = bathy['height']
mooringlons = [153.90, 154.002, 154.134, 154.286, 154.639, 155.30230]
mooringlats = [-27.325, -27.316, -27.281, -27.244, -27.202, -27.1026]

def figure():
   fig = plt.figure(figsize = (20, 8))
   grd = gs.GridSpec(2, 3, hspace = .3, width_ratios = [1,1,.8], height_ratios = [.07, 1])
   axs = [fig.add_subplot(grd[:,0]),
          fig.add_subplot(grd[:,1]),
          fig.add_subplot(grd[1,2])]
   for ax in axs[:-1]:
        ax.plot(coast[:,0], coast[:,1], 'k')
        ax.contourf(lon_bathy, lat_bathy, depth.where(depth>0), 
                    colors = 'gainsboro', linestyles = 'solid')
        ax.scatter(mooringlons, mooringlats, c = 'k', s = 10, zorder = 2)
        ax.scatter(153, -27.5, marker = 'o', c = 'dimgrey', s = 40, zorder = 2)
        ax.text(152.1, -27.8, 'Brisbane', c = 'dimgrey', fontsize = 8)
        ax.set_xlim(152, 158)
        ax.set_ylim(-32, -25)
   axs[2].grid(linestyle = '--', color = 'lightgrey')
   axs[2].tick_params(axis = 'x', color = 'C0', labelcolor = 'C0')
   axs[2].set_xlabel('Temperature ($^{\circ}$C)', color = 'C0')
   axss = axs[2].twiny()
   axss.tick_params(axis = 'x', color = 'orangered', labelcolor = 'orangered')
   axss.grid(linestyle = '--', color = 'lightgrey')
   axss.set_xlabel('Salinity (PSU)', color = 'orangered')
   axs.append(axss)
   axso = axs[2].twiny()
   axso.spines['top'].set_position(("axes", 1.15))
   axso.tick_params(axis = 'x', color = 'k', labelcolor = 'k')
   axso.grid(linestyle = '--', color = 'lightgrey')
   axso.set_xlabel('Oxygen (uM)', color = 'k')
   axs.append(axso)
   return fig, axs

year = sys.argv[1] 

# Gets list of all ssh files in the catalog
path_to_catalog_ssh = 'https://thredds.aodn.org.au/thredds/catalog/IMOS/OceanCurrent/GSLA/' \
                      '/NRT00/'+year+'/catalog.html'
ssh_cat = TDSCatalog(path_to_catalog_ssh)
# However, the catalog link isn't the same path for download (opendap). These 
# are of the form: 
# https://thredds.aodn.org.au/thredds/dodsC/IMOS/OceanCurrent/GSLA/NRT00/2022
# so what we want is to replace '/catalog/ by dodsC in the path.
path_to_opendap_ssh = path_to_catalog_ssh[:36]+'dodsC'+path_to_catalog_ssh[43:-12]

# Get YY-MM-DD of all files in the catalog
ssh_dates = []
for i in range(0, len(ssh_cat.datasets), 1):
    year = str(ssh_cat.datasets[i])[21:25]
    month = str(ssh_cat.datasets[i])[25:27]
    day = str(ssh_cat.datasets[i])[27:29]
    ssh_dates.append(year+'-'+month+'-'+day)
ssh_dates = np.array(ssh_dates)

# Repeat for SST
path_to_catalog_sst = 'https://thredds.aodn.org.au/thredds/catalog/IMOS/SRS/SST/ghrsst/L3S-1d/ngt/' \
                      +year+'/catalog.html'
sst_cat = TDSCatalog(path_to_catalog_sst)
# However, the catalog link isn't the same path for download (opendap). These 
# are of the form: 
# https://thredds.aodn.org.au/thredds/dodsC/IMOS/OceanCurrent/GSLA/NRT00/2022
# so what we want is to replace '/catalog/ by dodsC in the path.
path_to_opendap_sst = path_to_catalog_sst[:36]+'dodsC'+path_to_catalog_sst[43:-12]

# Get YY-MM-DD of all files in the catalog
sst_dates = []
for i in range(0, len(sst_cat.datasets), 1):
    year = str(sst_cat.datasets[i])[:4]
    month = str(sst_cat.datasets[i])[4:6]
    day = str(sst_cat.datasets[i])[6:8]
    sst_dates.append(year+'-'+month+'-'+day)
sst_dates = np.array(sst_dates)


# Get YY-MM-DD of profile
voyage = sys.argv[2]

if int(year) > 2019:
    path_to_profiles = glob.glob(r'../../../'+voyage+'/ctd/processing/'+voyage+'/cap/cappro/avg/*.nc')
else:
    path_to_profiles = glob.glob(r'../Past_Voyage_CTD/'+voyage+'/*/*.nc')
    
for i in range(0, len(path_to_profiles), 1):
    profile = xr.open_dataset(path_to_profiles[i])
    date = str(profile['time'].values[0])[:10]

    print('\n \nDeployment N: '+profile.attrs['Deployment']+'\n')
    print('CTD date: '+date)
    
    fig, axs = figure()
    # Find where they are the same
    idx = np.where(ssh_dates == date)[0]
    if len(idx) != 0:
        ssh = xr.open_dataset(path_to_opendap_ssh+str(ssh_cat.datasets[idx[0]])[:-3]+'.gz')
        ssh_r = ssh.sel(LONGITUDE = slice(151, 159), LATITUDE = slice(-33, -24))
        c = axs[0].pcolormesh(ssh_r['LONGITUDE'], ssh_r['LATITUDE'], ssh_r['GSLA'].squeeze(),
                            cmap = cm.cm.delta)
        plt.colorbar(c, ax = axs[0], orientation = 'horizontal').set_label('SSH (m)')
        axs[0].quiver(ssh_r['LONGITUDE'], ssh_r['LATITUDE'], ssh_r['UCUR'].squeeze(), 
                    ssh_r['VCUR'].squeeze(), units='width',scale = 10, width=0.004, 
                    headwidth =3)
        print('\nSSH date: '+str(ssh_r['TIME'].values[0])[:10])
        del ssh, ssh_r
    else:
        print('No SSH field for this day')
    idx = np.where(sst_dates == date)[0]
    if len(idx) != 0:
        sst = xr.open_dataset(path_to_opendap_sst+str(sst_cat.datasets[idx[0]]))
        sst_r = sst.sel(lon = slice(151, 159), lat = slice(-24, -33))
        c = axs[1].pcolormesh(sst_r['lon'], sst_r['lat'], sst_r['sea_surface_temperature'].squeeze()-273.15,
                            cmap = cm.cm.thermal)
        plt.colorbar(c, ax = axs[1], orientation = 'horizontal').set_label('SST ($^{\circ}}$C)')
        print('\nSST date: '+str(sst_r['time'].values[0])[:10])
        del sst, sst_r
    else:
        print('No SST field for this day')

    for ax in axs[:-3]:
        ax.scatter(profile['longitude'].values[0], profile['latitude'].values[0],
                c = 'm', s = 60)
    axs[2].plot(profile['temperature'].squeeze(), profile['pressure'].squeeze())
    axs[3].plot(profile['salinity'].squeeze(), profile['pressure'].squeeze(), color = 'orangered')
    axs[4].plot(profile['oxygen'].squeeze(), profile['pressure'].squeeze(), color = 'k')
    axs[2].set_ylim(0, profile['pressure'][-1])
    axs[2].axes.invert_yaxis()
    fig.suptitle('Voyage: '+voyage+ ' Deployment N: '+profile.attrs['Deployment']+
                '\n Date: '+date+' Lat: '+str(np.round(profile['latitude'].values[0], 2))+
                ' Lon: '+str(np.round(profile['longitude'].values[0], 2)),
                x = 0.4)
    plt.savefig('../FIGURES/EAC_class/'+voyage+'/'+path_to_profiles[i][-15:-9]+'.jpg',
                bbox_inches = 'tight')
    plt.close()

    print('\n'+path_to_profiles[i][-15:-9]+' done!')