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
    ax = fig.add_subplot()
    c = ax.contour(X, Y, Z, levels = 20, colors = ['grey'], ls = '--', zorder = 0)
    ax.clabel(c, fmt = '%0.1f')
    ax.set_xlim(34.45, 36.1)
    ax.set_ylim(0, 27)
    ax.set_xlabel('Absolute salinity', fontsize = 14)
    ax.set_ylabel('Conservative temp. ($^{\circ}$C)', fontsize = 14)
    return fig, ax

def figure_6():
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
        c = ax.contour(X, Y, Z, levels = 20, colors = ['grey'], ls = '--', zorder = 0)
        ax.clabel(c, fmt = '%0.1f')
        ax.set_xlim(34.45, 36.1)
        ax.set_ylim(0, 27)
        ax.set_xlabel('Absolute salinity');
        ax.set_ylabel('Conservative temp. ($^{\circ}$C)');
    return fig, axs

path = []
for v in ['in2022_v06', 'in2021_v03']:
    path.append(glob.glob('../../../'+v+'/ctd/processing/'+v+ \
                          '/cap/cappro/avg/*.nc'))
for v in ['in2019_v05', 'in2018_v03', 'in2016_v06']:
    path.append(glob.glob('../Past_Voyage_CTD/'+v+'/*/*.nc'))

legend_elements = [Line2D([0], [0], color = 'royalblue', lw = 5, label = 'in2022_v06'),
                   Line2D([0], [0], color = 'violet', lw = 5, label = 'in2021_v03'),
                   Line2D([0], [0], color = 'forestgreen', lw = 5, label = 'in2019_v05'),
                   Line2D([0], [0], color = 'navy', lw = 5, label = 'in2018_v03'),
                   Line2D([0], [0], color = 'coral', lw = 5, label = 'in2016_v06')]

fig, axs = figure_6()
past = 0
for i, clr in zip(range(0, len(path), 1), ['royalblue', 'violet', 'forestgreen', 'navy', 'coral']):
    paths = path[i]
    for j in np.arange(0, len(paths), 1):
        data = xr.open_dataset(paths[j])
        sa = gsw.SA_from_SP(data['salinity'], data['pressure'], data['longitude'], data['latitude'])
        ct = gsw.CT_from_t(sa, data['temperature'], data['pressure'])
        for ax in axs[i:]:
            if ax == axs[i]:
                axs[i].scatter(sa, ct, color = clr, s = 8,
                            edgecolor = 'none', zorder = 2);
                past += 1
            else:
                ax.scatter(sa, ct, color = 'grey', s = 8,
                           edgecolor = 'none', zorder = 2);

fig.legend(handles = legend_elements, bbox_to_anchor = (.85, .3));
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
    plt.savefig('../FIGURES/EAC_class/ts_'+t+'_highlight.jpg', bbox_inches = 'tight')
    plt.show()
