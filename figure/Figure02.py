import numpy as np
import pandas as pd
import xarray as xr
import sacpy as scp

import re, os, sys, glob
sys.path.append('../library')
from cbasic import *
from linrm  import *
from cdata  import *
from cplot  import *

import cmocean
import cartopy.crs as ccrs
import matplotlib.gridspec
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

####################################################################################################################################
projection_map = ccrs.PlateCarree(central_longitude=180); transform_map = ccrs.PlateCarree()
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
          
####################################################################################################################################
fig = plt.figure(figsize=(13,7)); gs = gridspec.GridSpec(nrows=2, ncols=2, height_ratios=[1,1], width_ratios=[1,1])

syear = '1977'; eyear = '1999'
domain = [18, 62, 110, 250]; cmap = cmap_white_center(plt.cm.PiYG_r); levs = np.arange(-1.2,1.21,0.12)
ax = plt.subplot(gs[0,0], projection=projection_map); ax.set_title('a', weight='bold', fontsize=20, x=-0.07, y=1.07)
ax.set_title('[1977-1999] STKA', loc='left', fontsize=13)#; ax.set_title('D(0)J(+1)', loc='right', fontsize=13)
STA_reg = xr.Dataset(); STA_reg['slope'] = P1_STA_slope; STA_reg['p_value'] = P1_STA_p_value; clim = data_domain(clim_STA, domain)
cb = reg_map(ax, STA_reg, domain, cmap, levs, transform_map); cb.ax.set_title(r'(K m s$^{-1}$)', pad=13, fontsize=11)
cs = ax.contour(clim.lon, clim.lat, clim, levels=[4,8,12], colors='grey', linewidths=1.3, transform=transform_map)
ax.clabel(cs, inline=True, fontsize=13, fmt='%d', colors='k')
STL  = clim.lat[clim.argmax(axis=0)].sel(lon=slice(140,238)).rolling(lon=4, center=True).mean()
ax.plot(STL.lon, STL, c='k', linewidth=2.4, transform=transform_map)
ccrs_grid(ax, np.arange(120,240.1,30), np.arange(20,60.1,10), 12)
#ax.text(-0.13, 0.5, 'ERA5: 1977-1999', va='center', ha='right', transform=ax.transAxes, weight='bold', fontsize=13, rotation=90)

domain = [18, 62, 110, 250]; cmap = cmap_white_center(plt.cm.coolwarm); levs = np.arange(-2,2.1,0.2)
ax = plt.subplot(gs[0,1], projection=projection_map); ax.set_title('b', weight='bold', fontsize=20, x=-0.07, y=1.07)
ax.set_title('[1977-1999] SLPA/UV$_{925}$A', loc='left', fontsize=13)#; ax.set_title('D(0)J(+1)', loc='right', fontsize=13)
SLP_reg = xr.Dataset(); SLP_reg['slope'] = P1_SLP_slope; SLP_reg['p_value'] = P1_SLP_p_value
UWD_reg = xr.Dataset(); UWD_reg['slope'] = P1_UWD_slope; UWD_reg['p_value'] = P1_UWD_p_value
VWD_reg = xr.Dataset(); VWD_reg['slope'] = P1_VWD_slope; VWD_reg['p_value'] = P1_VWD_p_value
cb = reg_map(ax, SLP_reg, domain, cmap, levs, transform_map); cb.ax.set_title('(hPa)', pad=13, fontsize=11)
qv = quiver_reg(ax, UWD_reg, VWD_reg, domain, transform_map, step=4, scale=13, siglev=0.1)
ax.quiverkey(qv, X=0.85, Y=1.05, U=1, label=r'1 m s$^{-1}$', labelpos='E')
ccrs_grid(ax, np.arange(120,240.1,30), np.arange(20,60.1,10), 12); ccrs_plot(ax, 140, 175, 31, 45)

syear = '2000'; eyear = '2022'
domain = [18, 62, 110, 250]; cmap = cmap_white_center(plt.cm.PiYG_r); levs = np.arange(-1.2,1.21,0.12)
ax = plt.subplot(gs[1,0], projection=projection_map); ax.set_title('c', weight='bold', fontsize=20, x=-0.07, y=1.07)
ax.set_title('[2000-2022] STKA', loc='left', fontsize=13)#; ax.set_title('D(0)J(+1)', loc='right', fontsize=13)
STA_reg = xr.Dataset(); STA_reg['slope'] = P2_STA_slope; STA_reg['p_value'] = P2_STA_p_value; clim = data_domain(clim_STA, domain)
cb = reg_map(ax, STA_reg, domain, cmap, levs, transform_map); cb.ax.set_title(r'(K m s$^{-1}$)', pad=13, fontsize=11)
cs = ax.contour(clim.lon, clim.lat, clim, levels=[4,8,12], colors='grey', linewidths=1.3, transform=transform_map)
ax.clabel(cs, inline=True, fontsize=13, fmt='%d', colors='k')
STL  = clim.lat[clim.argmax(axis=0)].sel(lon=slice(140,238)).rolling(lon=4, center=True).mean()
ax.plot(STL.lon, STL, c='k', linewidth=2.4, transform=transform_map)
ccrs_grid(ax, np.arange(120,240.1,30), np.arange(20,60.1,10), 12)
#ax.text(-0.13, 0.5, 'ERA5: 2000-2022', va='center', ha='right', transform=ax.transAxes, weight='bold', fontsize=13, rotation=90)

domain = [18, 62, 110, 250]; cmap = cmap_white_center(plt.cm.coolwarm); levs = np.arange(-2,2.1,0.2)
ax = plt.subplot(gs[1,1], projection=projection_map); ax.set_title('d', weight='bold', fontsize=20, x=-0.07, y=1.07)
ax.set_title('[2000-2022] SLPA/UV$_{925}$A', loc='left', fontsize=13)#; ax.set_title('D(0)J(+1)', loc='right', fontsize=13)
SLP_reg = xr.Dataset(); SLP_reg['slope'] = P2_SLP_slope; SLP_reg['p_value'] = P2_SLP_p_value
UWD_reg = xr.Dataset(); UWD_reg['slope'] = P2_UWD_slope; UWD_reg['p_value'] = P2_UWD_p_value
VWD_reg = xr.Dataset(); VWD_reg['slope'] = P2_VWD_slope; VWD_reg['p_value'] = P2_VWD_p_value
cb = reg_map(ax, SLP_reg, domain, cmap, levs, transform_map); cb.ax.set_title('(hPa)', pad=13, fontsize=11)
qv = quiver_reg(ax, UWD_reg, VWD_reg, domain, transform_map, step=4, scale=13, siglev=0.1)
ax.quiverkey(qv, X=0.85, Y=1.05, U=1, label=r'1 m s$^{-1}$', labelpos='E')
ccrs_grid(ax, np.arange(120,240.1,30), np.arange(20,60.1,10), 12); ccrs_plot(ax, 140, 175, 31, 45)

plt.tight_layout(w_pad=3, h_pad=3)
plt.show()
