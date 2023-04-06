import cmocean
import glob
import gsw
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.lines import Line2D
import seaborn as sns
import os
plt.rcParams['font.size'] = 12

tags = xr.open_mfdataset(glob.glob('data/CTD_data/*.nc'), combine = 'nested')
tags = tags.rename({'__xarray_dataarray_variable__': 'Class'})
#clr = dict({'EAC':'orange', 'WCE':'orangered', 'CCE':'b', 
#            'SHW':'sienna', 'NAN':'k'})
colors_set = sns.color_palette("coolwarm",11) # plt.style.use('seaborn-dark-palette')
clr = dict({'ss2012_v01':colors_set[0],'ss2013_v05':colors_set[1],'in2015_v02':colors_set[2], 'in2016_v04':colors_set[3], 
            'in2016_v06':colors_set[4], 
            'in2018_v03':colors_set[7], 'in2019_v05':colors_set[8],'in2021_v03':colors_set[9],
            'in2022_v06':colors_set[10]})
months = dict({'ss2012_v01':'April','ss2013_v05':'Aug','in2015_v02':'May', 'in2016_v04':'Sep', 'in2016_v06':'Oct', 
            'in2018_v03':'April', 'in2019_v05':'Sep','in2021_v03':'May','in2022_v06':'Jul'})

def figure_sa():
    fig = plt.figure(figsize = (14,12), facecolor=(1, 1, 1))
    grd = gs.GridSpec(2, 5, hspace = .22, wspace=.25)
    axs = [fig.add_subplot(grd[0,0]),
           fig.add_subplot(grd[0,1]),
           fig.add_subplot(grd[0,2]),
           fig.add_subplot(grd[0,3]),
           fig.add_subplot(grd[0,4]),
           fig.add_subplot(grd[1,0]),
           fig.add_subplot(grd[1,1]),
           fig.add_subplot(grd[1,2]),
           fig.add_subplot(grd[1,3])]
    for ax in axs:
        ax.set_facecolor("white")
        ax.set_xlim(34.45, 36.1)
        ax.set_ylim(0, 4500)
        #ax.xaxis.tick_top()
        ax.set_xlabel('Absolute salinity (g/kg)');
        #ax.xaxis.set_label_position('top');
        ax.set_ylabel('Pressure');
        ax.invert_yaxis()
        ax.grid()
    return fig, axs


def figure_ct():
    fig = plt.figure(figsize = (14,12))
    grd = gs.GridSpec(2, 3, hspace = .22, wspace=.25)
    axs = [fig.add_subplot(grd[0,0]),
           fig.add_subplot(grd[0,1]),
           fig.add_subplot(grd[0,2]),
           fig.add_subplot(grd[1,0]),
           fig.add_subplot(grd[1,1])]
    for ax in axs:
        ax.set_xlim(0, 27)
        ax.set_ylim(0, 4500)
        ax.xaxis.tick_top()
        ax.set_xlabel('Conservative temp. ($^{\circ}$C)');
        ax.xaxis.set_label_position('top');
        ax.set_ylabel('Pressure');
        ax.invert_yaxis()
        ax.grid()
    return fig, axs

def figure_oxy():
    fig = plt.figure(figsize = (14,12))
    grd = gs.GridSpec(2, 3, hspace = .22, wspace=.25)
    axs = [fig.add_subplot(grd[0,0]),
           fig.add_subplot(grd[0,1]),
           fig.add_subplot(grd[0,2]),
           fig.add_subplot(grd[1,0]),
           fig.add_subplot(grd[1,1])]
    for ax in axs:
        ax.set_xlim(150, 250)
        ax.set_ylim(0, 4500)
        ax.xaxis.tick_top()
        ax.set_xlabel('Oxygen (micromol l$^{-1}$)');
        ax.xaxis.set_label_position('top');
        ax.set_ylabel('Pressure');
        ax.invert_yaxis()
        ax.grid()
    return fig, axs


# legend_elements = [Line2D([0], [0], color = clr[t], lw = 5, label = t)
#                    for t in ['WCE', 'CCE', 'EAC', 'SHW', 'NAN']]

#####  SALINITY #####
## Plot salinity coloured by classification
fig, axs = figure_sa()
# Assign each axs to a different class
#ass_axs = dict({'EAC':0, 'WCE':1, 'CCE':2, 'SHW':3, 'NAN':4})
#n = dict({'EAC':0, 'WCE':0, 'CCE':0, 'SHW':0, 'NAN':0})
ass_axs = dict({'ss2012_v01':0,'ss2013_v05':1,'in2015_v02':2, 'in2016_v04':3, 'in2016_v06':4, 
            'in2018_v03':5, 'in2019_v05':6,'in2021_v03':7,
            'in2022_v06':8})
n = dict({'ss2012_v01':0,'ss2013_v05':0,'in2015_v02':0, 'in2016_v04':0, 'in2016_v06':0, 
            'in2018_v03':0, 'in2019_v05':0,'in2021_v03':0,
            'in2022_v06':0})
for i in range(0, len(tags['profile']), 1):
#   data = xr.open_dataset('data/CTD_data/'+tags['profile'][i].item()[:10]+ \
                             # '/'+tags['profile'][i].item()[:]+'.nc');
    data = xr.open_dataset('data/CTD_data/'+tags['profile'][i].item()[:10].upper()+ \
                                '/'+tags['profile'][i].item()+'.nc')
    
    sa = gsw.SA_from_SP(data['salinity'], data['pressure'], 
                        data['longitude'], data['latitude'])
    ct = gsw.CT_from_t(sa, data['temperature'], data['pressure'])

    profile_class = tags.isel(profile=i)['Class'].values.item();
    profile_voyage = tags.isel(profile=i)['profile'].values.item()[0:10]
    profile_year = tags.isel(profile=i)['profile'].values.item()[2:6]
    profile_month = pd.to_datetime(data.time[0].values).month

    axs[ass_axs[profile_voyage]].scatter(sa[:,0], data.pressure, color = clr[profile_voyage],
                                        s = 8, zorder=2);
    n[profile_voyage] = n[profile_voyage] + 1
    for ax in axs:
        if ax != axs[ass_axs[profile_voyage]]:
            ax.scatter(sa[:,0], data.pressure, color = 'grey', s = 8, zorder=1);
    
for title in ['ss2012_v01','ss2013_v05','in2015_v02', 'in2016_v04', 'in2016_v06', 
            'in2018_v03', 'in2019_v05','in2021_v03','in2022_v06']:
    axs[ass_axs[title]].text(0.55, 0.05, 'N = '+str(n[title]), 
                              fontsize=12,
                              transform=axs[ass_axs[title]].transAxes)
    
    axs[ass_axs[title]].set_title(title + months[title], fontsize=16);
#fig.legend(handles = legend_elements, bbox_to_anchor = (.75, .48), fontsize=14);
plt.savefig('results/sa_depth_profiles_by_voyage.jpg', bbox_inches = 'tight')


################################################################################################
###### MEAN per cruise
list_cruise = next(os.walk('data/CTD_data'))[1]

bins = np.arange(0,5000,2)
cruise_prof_ave_sa = np.empty((len(list_cruise), len(bins),)) * np.nan
cruise_prof_ave_ct = np.empty((len(list_cruise), len(bins),)) * np.nan
cruise_prof_ave_oxy = np.empty((len(list_cruise), len(bins),)) * np.nan

for l in range(0, len(list_cruise), 1):
    list_cruise_prof = os.listdir('data/CTD_data/'+ list_cruise[l] + '/' )
    list_cruise_prof_add = glob.glob('data/CTD_data/'+ list_cruise[l] + '/*.nc')
    sa_all = []
    ct_all = []
    oxy_all = []
    pressure_all = []
    for i in range(0, len(list_cruise_prof_add), 1):
        data = xr.open_dataset(list_cruise_prof_add[i])
        print(data)

        sa = gsw.SA_from_SP(data['salinity'], data['pressure'], 
                                data['longitude'], data['latitude'])
        ct = gsw.CT_from_t(sa, data['temperature'], data['pressure'])
        sa_all = np.append(sa_all, sa[:,0])
        ct_all = np.append(ct_all, ct[:,0])
        oxy_all = np.append(oxy_all, data.oxygen)
        pressure_all = np.append(pressure_all, data.pressure)
           
    A = np.vstack((np.digitize(pressure_all, bins), sa_all)).T
    res_lst = [np.nanmean(A[A[:, 0] == i, 1]) for i in range(len(bins))]
    cruise_prof_ave_sa[l,:] = res_lst
    
    A = np.vstack((np.digitize(pressure_all, bins), ct_all)).T
    res_lst = [np.nanmean(A[A[:, 0] == i, 1]) for i in range(len(bins))]
    cruise_prof_ave_ct[l,:] = res_lst
    
    A = np.vstack((np.digitize(pressure_all, bins), oxy_all)).T
    res_lst = [np.nanmean(A[A[:, 0] == i, 1]) for i in range(len(bins))]
    cruise_prof_ave_oxy[l,:] = res_lst


plt.plot(cruise_prof_ave_oxy.T)
plt.plot(cruise_prof_ave_ct.T)
    # plt.scatter(sa_all,pressure_all)
    # plt.scatter(res_lst,bins)
cruise_months = dict({'ss2012_v01':'ss2012_v01, April','ss2013_v05':'ss2013_v05, Aug','in2015_v02':'in2015_v02, May', 
                      'in2016_v04':'in2016_v04, Sep', 'in2016_v06':'in2016_v06, Oct', 
            'in2018_v03':'in2018_v03, April', 'in2019_v05':'in2019_v05, Sep','in2021_v03':'in2021_v03, May','in2022_v06':'in2022_v06, Jul'})    

########## FIGURES
fig, ax = plt.subplots(figsize = (5,8))
for i in range(0, len(tags['profile']), 1):
#   data = xr.open_dataset('data/CTD_data/'+tags['profile'][i].item()[:10]+ \
                             # '/'+tags['profile'][i].item()[:]+'.nc');
    data = xr.open_dataset('data/CTD_data/'+tags['profile'][i].item()[:10].upper()+ \
                                '/'+tags['profile'][i].item()+'.nc')
    
    sa = gsw.SA_from_SP(data['salinity'], data['pressure'], 
                        data['longitude'], data['latitude'])
    ct = gsw.CT_from_t(sa, data['temperature'], data['pressure'])

    profile_class = tags.isel(profile=i)['Class'].values.item();
    profile_voyage = tags.isel(profile=i)['profile'].values.item()[0:10]
    profile_year = tags.isel(profile=i)['profile'].values.item()[2:6]
    profile_month = pd.to_datetime(data.time[0].values).month

    ax.scatter(sa[:,0], data.pressure, color = 'grey', s = 8, zorder=1);
    
h2 = ax.plot(cruise_prof_ave_sa.T,bins)

ax.invert_yaxis()
ax.set_xlim(34.45, 36.1)
ax.set_ylim(0, 4500)
ax.set_xlabel('Absolute salinity (g/kg)');
ax.set_ylabel('Pressure');
ax.invert_yaxis()
ax.grid()
ax.legend(h2,cruise_months.values())
#fig.legend(handles = legend_elements, bbox_to_anchor = (.75, .48), fontsize=14);
plt.savefig('results/sa_depth_profiles_average.jpg', bbox_inches = 'tight')
plt.show()
    
   #### TEMP
fig, ax = plt.subplots(figsize = (5,8))
for i in range(0, len(tags['profile']), 1):
#   data = xr.open_dataset('data/CTD_data/'+tags['profile'][i].item()[:10]+ \
                             # '/'+tags['profile'][i].item()[:]+'.nc');
    data = xr.open_dataset('data/CTD_data/'+tags['profile'][i].item()[:10].upper()+ \
                                '/'+tags['profile'][i].item()+'.nc')
    
    sa = gsw.SA_from_SP(data['salinity'], data['pressure'], 
                        data['longitude'], data['latitude'])
    ct = gsw.CT_from_t(sa, data['temperature'], data['pressure'])

    profile_class = tags.isel(profile=i)['Class'].values.item();
    profile_voyage = tags.isel(profile=i)['profile'].values.item()[0:10]
    profile_year = tags.isel(profile=i)['profile'].values.item()[2:6]
    profile_month = pd.to_datetime(data.time[0].values).month

    ax.scatter(ct[:,0], data.pressure, color = 'grey', s = 8, zorder=1);
    
h2 = ax.plot(cruise_prof_ave_ct.T,bins)

ax.invert_yaxis()
ax.set_xlim(0, 27)
ax.set_ylim(0, 4500)
ax.set_xlabel('Absolute salinity (g/kg)');
ax.set_ylabel('Pressure');
ax.invert_yaxis()
ax.grid()
ax.legend(h2,cruise_months.values())
#fig.legend(handles = legend_elements, bbox_to_anchor = (.75, .48), fontsize=14);
plt.savefig('results/ct_depth_profiles_average.jpg', bbox_inches = 'tight')
plt.show()

 ### OXYGEN
 fig, ax = plt.subplots(figsize = (5,8))
 for i in range(0, len(tags['profile']), 1):
 #   data = xr.open_dataset('data/CTD_data/'+tags['profile'][i].item()[:10]+ \
                              # '/'+tags['profile'][i].item()[:]+'.nc');
     data = xr.open_dataset('data/CTD_data/'+tags['profile'][i].item()[:10].upper()+ \
                                 '/'+tags['profile'][i].item()+'.nc')
     
     sa = gsw.SA_from_SP(data['salinity'], data['pressure'], 
                         data['longitude'], data['latitude'])
     ct = gsw.CT_from_t(sa, data['temperature'], data['pressure'])

     profile_class = tags.isel(profile=i)['Class'].values.item();
     profile_voyage = tags.isel(profile=i)['profile'].values.item()[0:10]
     profile_year = tags.isel(profile=i)['profile'].values.item()[2:6]
     profile_month = pd.to_datetime(data.time[0].values).month

     ax.scatter(data.oxygen, data.pressure, color = 'grey', s = 8, zorder=1);
     
 h2 = ax.plot(cruise_prof_ave_oxy.T,bins)

 ax.invert_yaxis()
 #ax.set_xlim(0, 27)
 ax.set_ylim(0, 4500)
 ax.set_xlabel('Absolute salinity (g/kg)');
 ax.set_ylabel('Pressure');
 ax.invert_yaxis()
 ax.grid()
 ax.legend(h2,cruise_months.values())
 #fig.legend(handles = legend_elements, bbox_to_anchor = (.75, .48), fontsize=14);
 plt.savefig('results/oxy_depth_profiles_average.jpg', bbox_inches = 'tight')
 plt.show()   
    
# #####  OXYGEN  #####
# ## Plot salinity coloured by classification
# fig, axs = figure_oxy()
# # Assign each axs to a different class
# ass_axs = dict({'EAC':0, 'WCE':1, 'CCE':2, 'SHW':3, 'NAN':4})
# n = dict({'EAC':0, 'WCE':0, 'CCE':0, 'SHW':0, 'NAN':0})
# for i in range(0, len(tags['profile']), 1):
# #   data = xr.open_dataset('data/CTD_data/'+tags['profile'][i].item()[:10]+ \
#                              # '/'+tags['profile'][i].item()[:]+'.nc');
#     data = xr.open_dataset('data/CTD_data/'+tags['profile'][i].item()[:10].upper()+ \
#                                 '/'+tags['profile'][i].item()+'.nc')

#     sa = gsw.SA_from_SP(data['salinity'], data['pressure'], 
#                         data['longitude'], data['latitude'])
#     ct = gsw.CT_from_t(sa, data['temperature'], data['pressure'])
#     ox = data['oxygen'];

#     profile_class = tags.isel(profile=i)['Class'].values.item();

#     axs[ass_axs[profile_class]].scatter(ox, data.pressure, color = clr[profile_class],
#                                         s = 8, zorder=2);
#     n[profile_class] = n[profile_class] + 1
#     for ax in axs:
#         if ax != axs[ass_axs[profile_class]]:
#             ax.scatter(ox, data.pressure, color = 'grey', s = 8, zorder=1);
    
# for title in ['EAC', 'WCE', 'CCE', 'SHW', 'NAN']:
#     axs[ass_axs[title]].text(0.55, 0.05, 'N profiles = '+str(n[title]), 
#                               fontsize=12,
#                               transform=axs[ass_axs[title]].transAxes)
# fig.legend(handles = legend_elements, bbox_to_anchor = (.75, .48), fontsize=14);
# plt.savefig('results/ox_depth_profiles_by_classification.jpg', bbox_inches = 'tight')



# ## Now plot salinity coloured by oxygen
# fig, axs = figure_sa()
# # Assign each axs to a different class
# ass_axs = dict({'EAC':0, 'WCE':1, 'CCE':2, 'SHW':3, 'NAN':4})
# n = dict({'EAC':0, 'WCE':0, 'CCE':0, 'SHW':0, 'NAN':0})
# for i in range(0, len(tags['profile']), 1):
# #   data = xr.open_dataset('data/CTD_data/'+tags['profile'][i].item()[:10]+ \
#                              # '/'+tags['profile'][i].item()[:]+'.nc');
#     data = xr.open_dataset('data/CTD_data/'+tags['profile'][i].item()[:10].upper()+ \
#                                 '/'+tags['profile'][i].item()+'.nc')

#     sa = gsw.SA_from_SP(data['salinity'], data['pressure'], 
#                         data['longitude'], data['latitude'])
#     ct = gsw.CT_from_t(sa, data['temperature'], data['pressure'])
#     ox = data['oxygen'];

#     profile_class = tags.isel(profile=i)['Class'].values.item();

#     c = axs[ass_axs[profile_class]].scatter(sa[:,0], data.pressure, c = ox, cmap = cmocean.cm.solar,
#                                         vmin = 150, vmax = 250,
#                                         s = 8, zorder=2);
#     n[profile_class] = n[profile_class] + 1
#     for ax in axs:
#         if ax != axs[ass_axs[profile_class]]:
#             ax.scatter(sa[:,0], data.pressure, color = 'grey', s = 8, zorder=1);
    
# for title in ['EAC', 'WCE', 'CCE', 'SHW', 'NAN']:
#     axs[ass_axs[title]].set_title(title)
#     axs[ass_axs[title]].text(0.55, 0.05, 'N profiles = '+str(n[title]), 
#                               fontsize=12,
#                               transform=axs[ass_axs[title]].transAxes)
# # add colorbar
# plt.colorbar(c, cax = fig.add_axes([0.68, 0.44, 0.2, 0.03]), orientation = 'horizontal').set_label('oxy (micromol l$^{-1}$)')
# plt.savefig('results/sa_depth_profiles_by_classification_woxy.jpg', bbox_inches = 'tight')


# ##### TEMPERATURE #####
# ## Plot temp coloured by classification
# fig, axs = figure_ct()
# # Assign each axs to a different class
# ass_axs = dict({'EAC':0, 'WCE':1, 'CCE':2, 'SHW':3, 'NAN':4})
# n = dict({'EAC':0, 'WCE':0, 'CCE':0, 'SHW':0, 'NAN':0})
# for i in range(0, len(tags['profile']), 1):
# #   data = xr.open_dataset('data/CTD_data/'+tags['profile'][i].item()[:10]+ \
#                              # '/'+tags['profile'][i].item()[:]+'.nc');
#     data = xr.open_dataset('data/CTD_data/'+tags['profile'][i].item()[:10].upper()+ \
#                                 '/'+tags['profile'][i].item()+'.nc')

#     sa = gsw.SA_from_SP(data['salinity'], data['pressure'], 
#                         data['longitude'], data['latitude'])
#     ct = gsw.CT_from_t(sa, data['temperature'], data['pressure'])

#     profile_class = tags.isel(profile=i)['Class'].values.item();

#     axs[ass_axs[profile_class]].scatter(ct[:,0], data.pressure, color = clr[profile_class],
#                                         s = 8, zorder=2);
#     n[profile_class] = n[profile_class] + 1
#     for ax in axs:
#         if ax != axs[ass_axs[profile_class]]:
#             ax.scatter(ct[:,0], data.pressure, color = 'grey', s = 8, zorder=1);
    
# for title in ['EAC', 'WCE', 'CCE', 'SHW', 'NAN']:
#     axs[ass_axs[title]].text(0.55, 0.05, 'N profiles = '+str(n[title]), 
#                               fontsize=12,
#                               transform=axs[ass_axs[title]].transAxes)
# fig.legend(handles = legend_elements, bbox_to_anchor = (.75, .48), fontsize=14);
# plt.savefig('results/ct_depth_profiles_by_classification.jpg', bbox_inches = 'tight')

# ## Now plot temperature coloured by oxygen
# fig, axs = figure_ct()
# # Assign each axs to a different class
# ass_axs = dict({'EAC':0, 'WCE':1, 'CCE':2, 'SHW':3, 'NAN':4})
# n = dict({'EAC':0, 'WCE':0, 'CCE':0, 'SHW':0, 'NAN':0})
# for i in range(0, len(tags['profile']), 1):
# #for i in range(0, 50, 1):
# #   data = xr.open_dataset('data/CTD_data/'+tags['profile'][i].item()[:10]+ \
#                              # '/'+tags['profile'][i].item()[:]+'.nc');
#    data = xr.open_dataset('data/CTD_data/'+tags['profile'][i].item()[:10].upper()+ \
#                                 '/'+tags['profile'][i].item()+'.nc')

#    sa = gsw.SA_from_SP(data['salinity'], data['pressure'], 
#                         data['longitude'], data['latitude'])
#    ct = gsw.CT_from_t(sa, data['temperature'], data['pressure'])
#    ox = data['oxygen'];

#    profile_class = tags.isel(profile=i)['Class'].values.item();

#    c = axs[ass_axs[profile_class]].scatter(ct[:,0], data.pressure, c = ox, cmap = cmocean.cm.solar,
#                                         vmin = 150, vmax = 250,
#                                         s = 8, zorder=2);
#    n[profile_class] = n[profile_class] + 1
#    for ax in axs:
#        if ax != axs[ass_axs[profile_class]]:
#            ax.scatter(ct[:,0], data.pressure, color = 'grey', s = 8, zorder=1);
    
# for title in ['EAC', 'WCE', 'CCE', 'SHW', 'NAN']:
#     axs[ass_axs[title]].set_title(title)
#     axs[ass_axs[title]].text(0.55, 0.05, 'N profiles = '+str(n[title]), 
#                              fontsize=12,
#                              transform=axs[ass_axs[title]].transAxes)
# # add colorbar
# plt.colorbar(c, cax = fig.add_axes([0.68, 0.44, 0.2, 0.03]), orientation = 'horizontal').set_label('oxy (micromol l$^{-1}$)')
# plt.savefig('results/ct_depth_profiles_by_classification_woxy_colormap.jpg', bbox_inches = 'tight')
