import numpy as np
import pandas as pd
import xarray as xr
import sacpy as scp

import re, os, sys, glob
sys.path.append('../library')
from cbasic import *
from cdata  import *
from cplot  import *

import cmocean
import cartopy.crs as ccrs
import matplotlib.gridspec
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

####################################################################################################################################
projection_map = ccrs.PlateCarree(central_longitude=180); transform_map = ccrs.PlateCarree()
colors = [
    (1.0, 0.5, 0.0),  # orange
    (1.0, 1.0, 0.0),  # yellow
    (0.5, 0.8, 0.4),  # greenish
    (0.0, 0.5, 1.0),  # blue
    (1.0, 1.0, 1.0)   # white
    ]
cmap_custom = mcolors.LinearSegmentedColormap.from_list('custom_cmap', colors, N=256)

plt.rcParams['figure.dpi']  = 200
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['xtick.labelsize'] = 11
plt.rcParams['ytick.labelsize'] = 11

def data_map(ax, data, domain, cmap, levs, ticks, transform_map):
    clim = data_domain(data, domain); cf = ccrs_contourf(ax, clim, levs, cmap, transform_map)
    cb = plt.colorbar(cf, shrink=0.7, pad=0.07, orientation='vertical', ticks=ticks)
    ax.set_aspect('auto')
    return cb

def reg_map(ax, data, domain, cmap, levs, transform_map, siglev=0.1):
    slope = data_domain(data.slope, domain); cf = ccrs_contourf(ax, slope, levs, cmap, transform_map)
    cb = plt.colorbar(cf, shrink=0.7, pad=0.07, orientation='vertical', ticks=levs[::5])
    p_val = data_domain(data.p_value, domain); ccrs_sig(p_val, transform_map, siglev)
    ax.set_aspect('auto')
    return cb

def ccrs_sig(sig, transform, siglev, hatch='///', color='k'):
    sig = plt.contourf(sig.lon, sig.lat, sig, levels=[0., siglev, 1.0], colors='none', transform=transform)
    for i, contour in enumerate(sig.collections):
        if i == 0: contour.set_hatch(hatch); contour.set_edgecolor(color); contour.set_linewidth(0.0)

def comp_map(ax, comp, pval, domain, cmap, levs, transform_map, siglev=0.1):
    comp = data_domain(comp, domain); cf = ccrs_contourf(ax, comp, levs, cmap, transform_map)
    cb = plt.colorbar(cf, shrink=0.7, pad=0.07, orientation='vertical', ticks=levs[::5])
    p_val = data_domain(pval, domain); ccrs_sig(p_val, transform_map, siglev)
    ax.set_aspect('auto')
    return cb
  
####################################################################################################################################
ds = xr.open_dataset("../process/Fig06_STA_clim.nc")
clim_STA = ds['clim']

ds = xr.open_dataset("../process/Fig06_KE_reg_TPR_P1.nc")
P1_TPR_slope = ds['TPR_slope']; P1_TPR_p_value = ds['TPR_p_value']
P1_FRQ_slope = ds['FRQ_slope']; P1_FRQ_p_value = ds['FRQ_p_value']
ds = xr.open_dataset("../process/Fig06_KE_reg_TPR_P2.nc")
P2_TPR_slope = ds['TPR_slope']; P2_TPR_p_value = ds['TPR_p_value']
P2_FRQ_slope = ds['FRQ_slope']; P2_FRQ_p_value = ds['FRQ_p_value']

####################################################################################################################################
fig = plt.figure(figsize=(13,7)); gs = gridspec.GridSpec(nrows=2, ncols=2, height_ratios=[1,1], width_ratios=[1,1])

syear = '1977'; eyear = '1999'; levs = np.arange(-0.8,0.81,0.08)
domain = [18, 62, 110, 250]; cmap = cmap_white_center(plt.cm.BrBG)
ax = plt.subplot(gs[0,0], projection=projection_map); ax.set_title('a', weight='bold', fontsize=20, x=-0.07, y=1.07)
ax.set_title('[1977-1999] PrA/AR$_{freq}$A', loc='left', fontsize=13); #ax.set_title('ERA5: 1977-1999', loc='right', fontsize=13)
FRQ_reg = data_domain(P1_FRQ_slope, domain); TPR_reg = xr.Dataset(); TPR_reg['slope'] = P1_TPR_slope; TPR_reg['p_value'] = P1_TPR_p_value; clim = data_domain(clim_STA, domain)
cb = reg_map(ax, TPR_reg, domain, cmap, levs, transform_map); cb.ax.set_title(r'(mm day$^{-1}$)', pad=13, fontsize=11)
levels = [-0.9, -0.6, -0.3, 0.3, 0.6, 0.9]; colors = ['blue', 'blue', 'blue', 'red', 'red', 'red']; styles = ['dashed', 'dashed', 'dashed', 'solid', 'solid', 'solid']
cs = ax.contour(FRQ_reg.lon, FRQ_reg.lat, FRQ_reg, levels=levels, colors=colors, linewidths=1, linestyles=styles, transform=transform_map)#; ax.clabel(cs, fmt='%0.1f', fontsize=9, colors='k')
ax.clabel(cs, inline=True, fontsize=11, fmt='%.1f', colors='k')
STL  = clim.lat[clim.argmax(axis=0)].sel(lon=slice(140,238)).rolling(lon=4, center=True).mean()
ax.plot(STL.lon, STL, c='k', linewidth=2.4, transform=transform_map)
ccrs_grid(ax, np.arange(120,240.1,30), np.arange(20,60.1,10), 12)#; ccrs_plot(ax, 225, 245, 35, 60)
#ax.text(-0.13, 0.5, 'ERA5', va='center', ha='right', transform=ax.transAxes, weight='bold', fontsize=13, rotation=90)

syear = '2000'; eyear = '2022'
ax = plt.subplot(gs[0,1], projection=projection_map); ax.set_title('b', weight='bold', fontsize=20, x=-0.07, y=1.07)
ax.set_title('[2000-2022] PrA/AR$_{freq}$A', loc='left', fontsize=13); #ax.set_title('ERA5: 2000-2022', loc='right', fontsize=13)
FRQ_reg = data_domain(P2_FRQ_slope, domain); TPR_reg = xr.Dataset(); TPR_reg['slope'] = P2_TPR_slope; TPR_reg['p_value'] = P2_TPR_p_value; clim = data_domain(clim_STA, domain)
cb = reg_map(ax, TPR_reg, domain, cmap, levs, transform_map); cb.ax.set_title(r'(mm day$^{-1}$)', pad=13, fontsize=11)
levels = [-0.9, -0.6, -0.3, 0.3, 0.6, 0.9]; colors = ['blue', 'blue', 'blue', 'red', 'red', 'red']; styles = ['dashed', 'dashed', 'dashed', 'solid', 'solid', 'solid']
cs = ax.contour(FRQ_reg.lon, FRQ_reg.lat, FRQ_reg, levels=levels, colors=colors, linewidths=1, linestyles=styles, transform=transform_map)
ax.clabel(cs, inline=True, fontsize=11, fmt='%.1f', colors='k')
STL  = clim.lat[clim.argmax(axis=0)].sel(lon=slice(140,238)).rolling(lon=4, center=True).mean()
ax.plot(STL.lon, STL, c='k', linewidth=2.4, transform=transform_map)
ccrs_grid(ax, np.arange(120,240.1,30), np.arange(20,60.1,10), 12)#; ccrs_plot(ax, 215, 245, 33, 43)
#ax.text(-0.13, 0.5, 'ERA5', va='center', ha='right', transform=ax.transAxes, weight='bold', fontsize=13, rotation=90)

plt.tight_layout(w_pad=3, h_pad=3)
plt.show()
