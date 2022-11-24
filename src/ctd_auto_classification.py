############################################################################
# Auto classify ctd profiles based on criteria, and plot for reference:    #
#  CCE: <= -0.2 m                                                          #
#  WCE: >= 0.2 m                                                           #
#  EAC: < -0.3 m/s                                                         #
#  SHW: < 200 m bathymetry                                                 #
# run as:                                                                  #
# python3 src/ctd_auto_classification.py 'voyage_year' 'voyage_code'       #
# -----------------------------------                                      #
# Author: Julia Neme & Hannah Dawson (in2022_v06)                          #
# Date:  3 AUgust 2022                                                     #
############################################################################

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

coast = np.loadtxt("data/eaccoast.dat")
bathy = xr.open_dataset('data/bathy_dbdb2_v30_AUSTRALIA.nc')
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

def classify_ctd(profile, ssh_r):

    # interpolate values to ctd profile location
    lon, lat = profile.longitude.mean().values, profile.latitude.mean().values
    ssh_interp = ssh_r.GSLA.interp(LONGITUDE=lon, LATITUDE=lat)[0].values
    v_interp = ssh_r.VCUR.interp(LONGITUDE=lon, LATITUDE=lat)[0].values
    bathy_interp = -bathy.height.interp(lon=lon, lat=lat).values

    # define ctd classification
    if ssh_interp <= -0.2:
        class_tag = 'CCE'
    elif ssh_interp >= 0.2:
        class_tag = 'WCE'
    elif v_interp < -0.3:
        class_tag = 'EAC'
    elif bathy_interp <= 200:
        class_tag = 'SHW'
    else:
        class_tag = 'NAN'

    return class_tag

# Get year of voyage
year = sys.argv[1] 
# Get YY-MM-DD of profile
voyage = sys.argv[2]
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

if int(year) == 2022:
    path_to_profiles = glob.glob(r'data/CTD_data/'+voyage+'/processed/ctd/*.nc')
else:
    path_to_profiles = glob.glob(r'data/CTD_data/'+voyage+'/*.nc')

# Define colormap
vcmap = cm.cm.balance
vcmap.set_under('cyan')
vcmap.set_over('m')

class_tags = []
profile_tags = []
tag_colors = dict({'CCE':'b','WCE':'r', 'EAC':'limegreen', 'SHW':'sienna', 'NAN':'yellow' })

for i in range(0, len(path_to_profiles), 1):
    profile = xr.open_dataset(path_to_profiles[i])
    date = str(profile['time'].values[0])[:10]

    #print('\n \nDeployment N: '+profile.attrs['Deployment']+'\n')
    print('CTD date: '+date)

    fig, axs = figure()
    # Find where they are the same
    idx = np.where(ssh_dates == date)[0]
    if len(idx) != 0:
        ssh = xr.open_dataset(path_to_opendap_ssh+str(ssh_cat.datasets[idx[0]])[:-3]+'.gz')
        ssh_r = ssh.sel(LONGITUDE = slice(151, 159), LATITUDE = slice(-33, -24))

        class_tags.append(classify_ctd(profile, ssh_r))
        if int(year) in [2012,2013,2015]:
            profile_tags.append(path_to_profiles[i][-19:-3])
        else:
            profile_tags.append(path_to_profiles[i][-22:-3])

        c = axs[0].pcolormesh(ssh_r['LONGITUDE'], ssh_r['LATITUDE'], ssh_r['GSLA'].squeeze(),
                            cmap = cm.cm.delta)
        plt.colorbar(c, ax = axs[0], orientation = 'horizontal').set_label('SSH (m)')
        axs[0].quiver(ssh_r['LONGITUDE'], ssh_r['LATITUDE'], ssh_r['UCUR'].squeeze(), 
                    ssh_r['VCUR'].squeeze(), units='width',scale = 10, width=0.004, 
                    headwidth =3)
        axs[0].contour(bathy['lon'], bathy['lat'], -bathy['height'], levels =[200], colors = ['k'], linewidths = [2],)
        c = axs[1].pcolormesh(ssh_r['LONGITUDE'], ssh_r['LATITUDE'], ssh_r['VCUR'].squeeze(),
                            cmap = vcmap, vmin=-0.3, vmax=0.3)
        axs[1].contour(bathy['lon'], bathy['lat'], -bathy['height'], levels =[200], colors = ['k'], linewidths = [2],)
        plt.colorbar(c, ax = axs[1], orientation = 'horizontal', extend='both').set_label('Surface meridional velocity ($ms^{-1}$)')
        print('\nSSH date: '+str(ssh_r['TIME'].values[0])[:10])
        del ssh, ssh_r
        
        for ax in axs[:-3]:
            ax.scatter(profile['longitude'].values[0], profile['latitude'].values[0],
                    c = tag_colors[class_tags[i]], s = 60);
        axs[2].plot(profile['temperature'].squeeze(), profile['pressure'].squeeze());
        axs[3].plot(profile['salinity'].squeeze(), profile['pressure'].squeeze(), color = 'orangered');
        axs[4].plot(profile['oxygen'].squeeze(), profile['pressure'].squeeze(), color = 'k');
        axs[2].set_ylim(0, profile['pressure'][-1]);
        axs[2].axes.invert_yaxis();
        if int(year) == 2015:
            #fig.suptitle('Voyage: '+voyage+ ' Deployment N: '+profile.attrs['Deployment']+
            #            '\n Date: '+date+' Lat: '+str(np.round(profile['latitude'].values[0], 2))+
            #            ' Lon: '+str(np.round(profile['longitude'].values[0], 2)),
            #            x = 0.4);
            plt.savefig('results/figures/EAC_class/'+path_to_profiles[i][-19:-6]+'.jpg',
                        bbox_inches = 'tight')
        elif int(year) in [2012, 2013]:
            plt.savefig('results/figures/EAC_class/'+path_to_profiles[i][-19:-6]+'.jpg',
                        bbox_inches = 'tight')
        else:
            plt.savefig('results/figures/EAC_class/'+path_to_profiles[i][-22:-9]+'.jpg',
                        bbox_inches = 'tight')
            
        plt.close()

        print('\n'+path_to_profiles[i][-15:-9]+' done!')
    else:
        print('No SSH field for this day, couldnt classify')

tagarr = xr.DataArray(class_tags, dims = 'profile', coords = {'profile':profile_tags})
tagarr.to_netcdf('data/CTD_data/'+voyage+'_classification.nc')