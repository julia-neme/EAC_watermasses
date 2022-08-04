import glob
import gsw
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib.lines import Line2D
plt.rcParams['font.size'] = 14

tags = xr.open_mfdataset(glob.glob('data/CTD_data/*.nc'), combine = 'nested')
tags = tags.rename({'__xarray_dataarray_variable__': 'Class'})
clr = dict({'EAC':'orange', 'WCE':'orangered', 'CCE':'b', 
            'SHW':'sienna', 'NAN':'k'})

def figure():
    x = np.arange(34.2, 36.11, .01)
    y = np.arange(0, 27.1, .1)
    X, Y = np.meshgrid(x, y)
    Z = gsw.density.sigma0(X, Y)

    fig = plt.figure(figsize = (18,10))
    grd = gs.GridSpec(2, 3, hspace = .3)
    axs = [fig.add_subplot(grd[0,0]),
           fig.add_subplot(grd[0,1]),
           fig.add_subplot(grd[0,2]),
           fig.add_subplot(grd[1,0]),
           fig.add_subplot(grd[1,1])]
    for ax in axs:
        c = ax.contour(X, Y, Z, levels = 20, colors = ['grey'], zorder = 0)
        ax.clabel(c, fmt = '%0.1f')
        ax.set_xlim(34.45, 36.1)
        ax.set_ylim(0, 27)
        ax.set_xlabel('Absolute salinity');
        ax.set_ylabel('Conservative temp. ($^{\circ}$C)');
    return fig, axs

legend_elements = [Line2D([0], [0], color = clr[t], lw = 5, label = t)
                   for t in ['WCE', 'CCE', 'EAC', 'SHW', 'NAN']]

fig, axs = figure()
# Assign each axs to a different class
ass_axs = dict({'EAC':0, 'WCE':1, 'CCE':2, 'SHW':3, 'NAN':4})
n = dict({'EAC':0, 'WCE':0, 'CCE':0, 'SHW':0, 'NAN':0})
for i in range(0, len(tags['profile']), 1):
    data = xr.open_dataset('data/CTD_data/'+tags['profile'][i].item()[:10]+ \
                             '/'+tags['profile'][i].item()+'.nc')

    sa = gsw.SA_from_SP(data['salinity'], data['pressure'], 
                        data['longitude'], data['latitude'])
    ct = gsw.CT_from_t(sa, data['temperature'], data['pressure'])

    profile_class = tags.isel(profile=i)['Class'].values.item()

    axs[ass_axs[profile_class]].scatter(sa, ct, color = clr[profile_class],
                                        s = 8, zorder=2);
    n[profile_class] = n[profile_class] + 1
    for ax in axs:
        if ax != axs[ass_axs[profile_class]]:
            ax.scatter(sa, ct, color = 'grey', s = 8, zorder=1);
    
for title in ['EAC', 'WCE', 'CCE', 'SHW', 'NAN']:
    axs[ass_axs[title]].set_title('N profiles = '+str(n[title]))
fig.legend(handles = legend_elements, bbox_to_anchor = (.83, .3));
plt.savefig('results/ts_auto_classification.jpg', bbox_inches = 'tight')


fig, axs = figure()
# Assign each axs to a different class
ass_axs = dict({'EAC':0, 'WCE':1, 'CCE':2, 'SHW':3, 'NAN':4})
n = dict({'EAC':0, 'WCE':0, 'CCE':0, 'SHW':0, 'NAN':0})
for i in range(0, len(tags['profile']), 1):
    data = xr.open_dataset('data/CTD_data/'+tags['profile'][i].item()[:10]+ \
                             '/'+tags['profile'][i].item()+'.nc')

    sa = gsw.SA_from_SP(data['salinity'], data['pressure'], 
                        data['longitude'].mean(), data['latitude'].mean())
    ct = gsw.CT_from_t(sa, data['temperature'], data['pressure'])
    ox = data['oxygen']
    profile_class = tags.isel(profile=i)['Class'].values.item()

    c = axs[ass_axs[profile_class]].scatter(sa, ct, c = ox, cmap = 'gist_ncar',
                                            vmin = 160, vmax = 210,
                                            s = 8, zorder=2);
    n[profile_class] = n[profile_class] + 1
    for ax in axs:
        if ax != axs[ass_axs[profile_class]]:
            ax.scatter(sa, ct, color = 'grey', s = 8, zorder=1);
plt.colorbar(c, cax = fig.add_axes([0.68, 0.3, 0.2, 0.05]), orientation = 'horizontal').set_label('oxy')
for title in ['EAC', 'WCE', 'CCE', 'SHW', 'NAN']:
    axs[ass_axs[title]].set_title(title+' - N profiles = '+str(n[title]))
plt.savefig('results/ts_auto_classification_woxy.jpg', bbox_inches = 'tight')
