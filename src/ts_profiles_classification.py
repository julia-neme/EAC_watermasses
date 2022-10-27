import glob
import gsw
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib.lines import Line2D
plt.rcParams['font.size'] = 12

tags = xr.open_mfdataset(glob.glob('data/CTD_data/*.nc'), combine = 'nested')
tags = tags.rename({'__xarray_dataarray_variable__': 'Class'})
clr = dict({'EAC':'orange', 'WCE':'orangered', 'CCE':'b','SHW':'sienna', 'NAN':'k'})

def figure_sa():
    fig = plt.figure(figsize = (14,12))
    grd = gs.GridSpec(2, 3, hspace = .22, wspace=.25)
    axs = [fig.add_subplot(grd[0,0]),
           fig.add_subplot(grd[0,1]),
           fig.add_subplot(grd[0,2]),
           fig.add_subplot(grd[1,0]),
           fig.add_subplot(grd[1,1])]
    for ax in axs:
        ax.set_xlim(34.45, 36.1)
        ax.set_ylim(0, 4500)
        ax.xaxis.tick_top()
        ax.set_xlabel('Absolute salinity (g/kg)');
        ax.xaxis.set_label_position('top');
        ax.set_ylabel('Pressure');
        ax.invert_yaxis()
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
    return fig, axs

legend_elements = [Line2D([0], [0], color = clr[t], lw = 5, label = t)
                   for t in ['WCE', 'CCE', 'EAC', 'SHW', 'NAN']]

#####  SALINITY #####
## Plot salinity coloured by classification
fig, axs = figure_sa()
# Assign each axs to a different class
ass_axs = dict({'EAC':0, 'WCE':1, 'CCE':2, 'SHW':3, 'NAN':4})
n = dict({'EAC':0, 'WCE':0, 'CCE':0, 'SHW':0, 'NAN':0})
for i in range(0, len(tags['profile']), 1):
    data = xr.open_dataset('data/CTD_data/'+tags['profile'][i].item()[:10]+ \
                             '/'+tags['profile'][i].item()[:]+'.nc');

    sa = gsw.SA_from_SP(data['salinity'], data['pressure'], 
                        data['longitude'], data['latitude'])
    ct = gsw.CT_from_t(sa, data['temperature'], data['pressure'])

    profile_class = tags.isel(profile=i)['Class'].values.item();

    axs[ass_axs[profile_class]].scatter(sa[:,0], data.pressure, color = clr[profile_class],
                                        s = 8, zorder=2);
    n[profile_class] = n[profile_class] + 1
    for ax in axs:
        if ax != axs[ass_axs[profile_class]]:
            ax.scatter(sa[:,0], data.pressure, color = 'grey', s = 8, zorder=1);
    
for title in ['EAC', 'WCE', 'CCE', 'SHW', 'NAN']:
    axs[ass_axs[title]].text(0.55, 0.05, 'N profiles = '+str(n[title]), 
                             fontsize=12,
                             transform=axs[ass_axs[title]].transAxes)
fig.legend(handles = legend_elements, bbox_to_anchor = (.75, .48), fontsize=14);
plt.savefig('results/sa_depth_profiles_by_classification.jpg', bbox_inches = 'tight')

## Now plot salinity coloured by oxygen
fig, axs = figure_sa()
# Assign each axs to a different class
ass_axs = dict({'EAC':0, 'WCE':1, 'CCE':2, 'SHW':3, 'NAN':4})
n = dict({'EAC':0, 'WCE':0, 'CCE':0, 'SHW':0, 'NAN':0})
for i in range(0, len(tags['profile']), 1):
    data = xr.open_dataset('data/CTD_data/'+tags['profile'][i].item()[:10]+ \
                             '/'+tags['profile'][i].item()[:]+'.nc');

    sa = gsw.SA_from_SP(data['salinity'], data['pressure'], 
                        data['longitude'], data['latitude'])
    ct = gsw.CT_from_t(sa, data['temperature'], data['pressure'])
    ox = data['oxygen'];

    profile_class = tags.isel(profile=i)['Class'].values.item();

    c = axs[ass_axs[profile_class]].scatter(sa[:,0], data.pressure, c = ox, cmap = 'gist_ncar',
                                        vmin = 160, vmax = 210,
                                        s = 8, zorder=2);
    n[profile_class] = n[profile_class] + 1
    for ax in axs:
        if ax != axs[ass_axs[profile_class]]:
            ax.scatter(sa[:,0], data.pressure, color = 'grey', s = 8, zorder=1);
    
for title in ['EAC', 'WCE', 'CCE', 'SHW', 'NAN']:
    axs[ass_axs[title]].set_title(title)
    axs[ass_axs[title]].text(0.55, 0.05, 'N profiles = '+str(n[title]), 
                             fontsize=12,
                             transform=axs[ass_axs[title]].transAxes)
# add colorbar
plt.colorbar(c, cax = fig.add_axes([0.68, 0.44, 0.2, 0.03]), orientation = 'horizontal').set_label('oxy')
plt.savefig('results/sa_depth_profiles_by_classification_woxy.jpg', bbox_inches = 'tight')


##### TEMPERATURE #####
## Plot temp coloured by classification
fig, axs = figure_ct()
# Assign each axs to a different class
ass_axs = dict({'EAC':0, 'WCE':1, 'CCE':2, 'SHW':3, 'NAN':4})
n = dict({'EAC':0, 'WCE':0, 'CCE':0, 'SHW':0, 'NAN':0})
for i in range(0, len(tags['profile']), 1):
    data = xr.open_dataset('data/CTD_data/'+tags['profile'][i].item()[:10]+ \
                             '/'+tags['profile'][i].item()[:]+'.nc');

    sa = gsw.SA_from_SP(data['salinity'], data['pressure'], 
                        data['longitude'], data['latitude'])
    ct = gsw.CT_from_t(sa, data['temperature'], data['pressure'])

    profile_class = tags.isel(profile=i)['Class'].values.item();

    axs[ass_axs[profile_class]].scatter(ct[:,0], data.pressure, color = clr[profile_class],
                                        s = 8, zorder=2);
    n[profile_class] = n[profile_class] + 1
    for ax in axs:
        if ax != axs[ass_axs[profile_class]]:
            ax.scatter(ct[:,0], data.pressure, color = 'grey', s = 8, zorder=1);
    
for title in ['EAC', 'WCE', 'CCE', 'SHW', 'NAN']:
    axs[ass_axs[title]].text(0.55, 0.05, 'N profiles = '+str(n[title]), 
                             fontsize=12,
                             transform=axs[ass_axs[title]].transAxes)
fig.legend(handles = legend_elements, bbox_to_anchor = (.75, .48), fontsize=14);
plt.savefig('results/ct_depth_profiles_by_classification.jpg', bbox_inches = 'tight')

## Now plot temperature coloured by oxygen
fig, axs = figure_ct()
# Assign each axs to a different class
ass_axs = dict({'EAC':0, 'WCE':1, 'CCE':2, 'SHW':3, 'NAN':4})
n = dict({'EAC':0, 'WCE':0, 'CCE':0, 'SHW':0, 'NAN':0})
for i in range(0, len(tags['profile']), 1):
    data = xr.open_dataset('data/CTD_data/'+tags['profile'][i].item()[:10]+ \
                             '/'+tags['profile'][i].item()[:]+'.nc');

    sa = gsw.SA_from_SP(data['salinity'], data['pressure'], 
                        data['longitude'], data['latitude'])
    ct = gsw.CT_from_t(sa, data['temperature'], data['pressure'])
    ox = data['oxygen'];

    profile_class = tags.isel(profile=i)['Class'].values.item();

    c = axs[ass_axs[profile_class]].scatter(ct[:,0], data.pressure, c = ox, cmap = 'gist_ncar',
                                        vmin = 160, vmax = 210,
                                        s = 8, zorder=2);
    n[profile_class] = n[profile_class] + 1
    for ax in axs:
        if ax != axs[ass_axs[profile_class]]:
            ax.scatter(ct[:,0], data.pressure, color = 'grey', s = 8, zorder=1);
    
for title in ['EAC', 'WCE', 'CCE', 'SHW', 'NAN']:
    axs[ass_axs[title]].set_title(title)
    axs[ass_axs[title]].text(0.55, 0.05, 'N profiles = '+str(n[title]), 
                             fontsize=12,
                             transform=axs[ass_axs[title]].transAxes)
# add colorbar
plt.colorbar(c, cax = fig.add_axes([0.68, 0.44, 0.2, 0.03]), orientation = 'horizontal').set_label('oxy')
plt.savefig('results/ct_depth_profiles_by_classification_woxy.jpg', bbox_inches = 'tight')
