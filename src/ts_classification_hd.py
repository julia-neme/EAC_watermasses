import cmocean as cm
import glob
import gsw
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import sys
import xarray as xr
from matplotlib.lines import Line2D
plt.rcParams['font.size'] = 14

tags = pd.read_csv('../FIGURES/EAC_class/CTD_profile_classification_EAC.csv').to_xarray()
clr = np.zeros(len(tags['index']))
clr = xr.where(tags['Class'] == 'EAC', 'orange', clr)
clr = xr.where(tags['Class'] == 'SH', 'sienna', clr)
clr = xr.where(tags['Class'] == 'CCE', 'royalblue', clr)
clr = xr.where(tags['Class'] == 'WCE', 'orangered', clr)

def figure():
    x = np.arange(34.2, 36.11, .01)
    y = np.arange(0, 27.1, .1)
    X, Y = np.meshgrid(x, y)
    Z = gsw.density.sigma0(X, Y)

    fig = plt.figure(figsize = (11, 10))
    grd = gs.GridSpec(2, 3, hspace = .2, wspace=.2, width_ratios = [1,1,.05])
    axs = [fig.add_subplot(grd[0,0]),
         fig.add_subplot(grd[0,1]),
         fig.add_subplot(grd[1,0]),
         fig.add_subplot(grd[1,1]),
         fig.add_subplot(grd[:,2]),]
    axcount = 0
    for ax in axs[:-1]:
         c = ax.contour(X, Y, Z, levels = 20, colors = ['grey'], ls = '--', zorder = 0)
         ax.clabel(c, fmt = '%0.1f')
         ax.set_xlim(34.45, 36.1)
         ax.set_ylim(0, 27)
         if axcount == 0 or axcount == 2:
             ax.set_ylabel('Conservative temp. ($^{\circ}$C)', fontsize = 14)
         if axcount == 2 or axcount == 3:
             ax.set_xlabel('Absolute salinity', fontsize = 14)
         axcount += 1

    return fig, axs


path = []
for v in ['in2022_v06', 'in2021_v03']:
    path.append(glob.glob('../../../'+v+'/ctd/processing/'+v+ \
                          '/cap/cappro/avg/*.nc'))
for v in ['in2019_v05', 'in2018_v03', 'in2016_v06']:
    path.append(glob.glob('../Past_Voyage_CTD/'+v+'/*/*.nc'))


fig, axs = figure()
count = 0
for t in ['EAC', 'CCE', 'WCE', 'SH']:
    legend_color = clr.where(tags['Class'] == t).dropna('index')[0].item()
    legend_elements = [Line2D([0], [0], color = legend_color, lw = 5, label = t)]

    for i in range(0, len(tags['index']), 1):
        year = str(tags['Year'][i].values)
        voyage = str(tags['Voyage'][i].values)
        cast = str(tags['CTD'][i].values)
        if int(year) > 2019:
            path = glob.glob('../../../in'+year+'_'+voyage+'/ctd/processing/in'+year+'_'+voyage+ \
                '/cap/cappro/avg/*'+cast+'*.nc')
        else:
            path = glob.glob('../Past_Voyage_CTD/in'+year+'_'+voyage+'/*/*'+cast+'*.nc')    
        data = xr.open_dataset(path[0])

        sa = gsw.SA_from_SP(data['salinity'], data['pressure'], data['longitude'], data['latitude'])
        ct = gsw.CT_from_t(sa, data['temperature'], data['pressure'])
        if tags['Class'][i] == t:
            c = axs[count].scatter(sa, ct, c = data['oxygen'], cmap = cm.cm.solar, s = 8,
                    vmin = 170, vmax=230,
                    edgecolor = 'none', zorder = 2);
        else:
            axs[count].scatter(sa, ct, color = 'darkgrey', s = 8,
                    edgecolor = 'none', zorder = 1);
    axs[count].set_title(t, fontsize=18);
    count += 1
cbar = plt.colorbar(c, orientation = 'vertical', cax = axs[-1], extend='both')
cbar.set_label('Oxygen (uM/l)') 
plt.savefig('../FIGURES/EAC_class/ts_oxygen.jpg', bbox_inches = 'tight', dpi=200)
plt.show()












# ---------------------------------------------------------------------------------------
legend_elements = [Line2D([0], [0], color = 'royalblue', lw = 5, label = 'in2022_v06'),
                   Line2D([0], [0], color = 'violet', lw = 5, label = 'in2021_v03'),
                   Line2D([0], [0], color = 'forestgreen', lw = 5, label = 'in2019_v05'),
                   Line2D([0], [0], color = 'navy', lw = 5, label = 'in2018_v03'),
                   Line2D([0], [0], color = 'coral', lw = 5, label = 'in2016_v06')]

fig, ax = figure()
for i, clr in zip(range(0, len(path), 1), ['royalblue', 'violet', 'forestgreen', 'navy', 'coral']):
    
    paths = path[i]

    for j in np.arange(0,len(paths),1):
        data = xr.open_dataset(paths[j])

        sa = gsw.SA_from_SP(data['salinity'], data['pressure'], data['longitude'], data['latitude'])
        ct = gsw.CT_from_t(sa, data['temperature'], data['pressure'])
        ax.scatter(sa, ct, color = clr, s = 8,
                    edgecolor = 'none', zorder = 2);
fig.legend(handles = legend_elements, fontsize = 14, bbox_to_anchor = (1.05, .8));
plt.savefig('../FIGURES/EAC_class/ts_per_cruise.jpg', bbox_inches = 'tight')
plt.show()


for t in ['EAC', 'CCE', 'WCE', 'SH']:
    legend_color = clr.where(tags['Class'] == t).dropna('index')[0].item()
    legend_elements = [Line2D([0], [0], color = legend_color, lw = 5, label = t)]

    fig, ax = figure()
    for i in range(0, len(tags['index']), 1):
        year = str(tags['Year'][i].values)
        voyage = str(tags['Voyage'][i].values)
        cast = str(tags['CTD'][i].values)
        if int(year) > 2019:
            path = glob.glob('../../../in'+year+'_'+voyage+'/ctd/processing/in'+year+'_'+voyage+ \
                '/cap/cappro/avg/*'+cast+'*.nc')
        else:
            path = glob.glob('../Past_Voyage_CTD/in'+year+'_'+voyage+'/*/*'+cast+'*.nc')    
        data = xr.open_dataset(path[0])

        sa = gsw.SA_from_SP(data['salinity'], data['pressure'], data['longitude'], data['latitude'])
        ct = gsw.CT_from_t(sa, data['temperature'], data['pressure'])
        if tags['Class'][i] == t:
            ax.scatter(sa, ct, color = clr[i].item(), s = 8,
                    edgecolor = 'none', zorder = 2);
        else:
            ax.scatter(sa, ct, color = 'darkgrey', s = 8,
                    edgecolor = 'none', zorder = 1);
    fig.legend(handles = legend_elements, fontsize = 14, bbox_to_anchor = (1.05, .8));
    #plt.savefig('../FIGURES/EAC_class/ts_'+t+'_highlight.jpg', bbox_inches = 'tight')
    plt.show()


