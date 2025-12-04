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
plt.rcParams['figure.dpi']  = 200
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['xtick.labelsize'] = 11
plt.rcParams['ytick.labelsize'] = 11

def grad_lat(data):
    dlat = data.lat.diff('lat').mean().values
    grad = data.diff('lat') / abs(dlat)
    lat_idx = grad.argmin(dim='lat')
    max_lat = data.lat.isel(lat=lat_idx)
    return max_lat
  
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
ds = xr.open_dataset("../process/Fig01_mvcorr_KE_KOE_SST.nc")
mvcorr_SST = ds['mvcorr']

ds = xr.open_dataset("../process/Fig01_mvcorr_KE_KOE_LHF.nc")
mvcorr_LHF = ds['mvcorr']

ds = xr.open_dataset("../process/Fig01_SST_clim.nc")
clim_SST = ds['clim']

ds = xr.open_dataset("../process/Fig01_KE_reg_SST_P1.nc")
P1_SST_slope = ds['SST_slope']; P1_SST_p_value = ds['SST_p_value']
ds = xr.open_dataset("../process/Fig01_KE_reg_LHF_P1.nc")
P1_LHF_slope = ds['LHF_slope']; P1_LHF_p_value = ds['LHF_p_value']

ds = xr.open_dataset("../process/Fig01_KE_reg_SST_P2.nc")
P2_SST_slope = ds['SST_slope']; P2_SST_p_value = ds['SST_p_value']
ds = xr.open_dataset("../process/Fig01_KE_reg_LHF_P2.nc")
P2_LHF_slope = ds['LHF_slope']; P2_LHF_p_value = ds['LHF_p_value']

####################################################################################################################################
fig = plt.figure(figsize=(13,12)); gs = gridspec.GridSpec(nrows=3, ncols=4, height_ratios=[0.7,1,1], width_ratios=[2.0,0.1,2.0,0.1])

# Fig. 1a
thr95 = 0.514; thr90 = 0.441; time = mvcorr_SST.time; mvcorr = mvcorr_SST
ax = plt.subplot(gs[0,0]); ax.set_title('a', weight='bold', fontsize=20, x=-0.07, y=1.07)
ax.set_title('SON KE SSH and subsequent DJ KOE SST', loc='left', fontsize=12.5)
ax.bar(time, mvcorr, color='white', edgecolor='C0')
mask90 = mvcorr > thr90; ax.bar(time[mask90], mvcorr[mask90], color='C0', alpha=0.5)
mask95 = mvcorr > thr95; ax.bar(time[mask95], mvcorr[mask95], color='C0', label='KOE SST')
ax.axhline(y=0, c='k', lw=0.7)#; ax.legend(loc='upper right', fontsize=12)
ax.set(xlim=[time.min()-1.5,time.max()+1.5], xticks=np.arange(1985,2015.1,5), ylim=[-0.7,1.0], yticks=np.arange(-0.6,0.91,0.3))
ax.set_xlabel('Central year of sliding window', fontsize=12)
ax.set_ylabel('Correlation coefficient', fontsize=12)

# Fig. 1b
thr95 = 0.514; thr90 = 0.441; time = mvcorr_LHF.time; mvcorr = mvcorr_LHF
ax = plt.subplot(gs[0,2]); ax.set_title('b', weight='bold', fontsize=20, x=-0.07, y=1.07)
ax.set_title('SON KE SSH and subsequent DJ KOE LHF', loc='left', fontsize=12.5)
ax.bar(time, mvcorr, color='white', edgecolor='C1')
mask90 = mvcorr > thr90; ax.bar(time[mask90], mvcorr[mask90], color='C1', alpha=0.5)
mask95 = mvcorr > thr95; ax.bar(time[mask95], mvcorr[mask95], color='C1', label='KOE LHF')
ax.axhline(y=0, c='k', lw=0.7)#; ax.legend(loc='upper left', fontsize=12)
ax.set(xlim=[time.min()-1.5,time.max()+1.5], xticks=np.arange(1985,2015.1,5), ylim=[-0.7,1.0], yticks=np.arange(-0.6,0.91,0.3))
ax.set_xlabel('Central year of sliding window', fontsize=12)
ax.set_ylabel('Correlation coefficient', fontsize=12)

syear = '1977'; eyear = '1999'
# Fig. 1c
domain = [-22, 62, 100, 290]; cmap = cmap_white_center(plt.cm.RdBu_r); levs = np.arange(-0.6,0.61,0.06)
ax = plt.subplot(gs[1,0:2], projection=projection_map); ax.set_title('c', weight='bold', fontsize=20, x=-0.07, y=1.07)
ax.set_title('[1977-1999] SSTA', loc='left', fontsize=13)#; ax.set_title('D(0)J(+1)', loc='right', fontsize=13)
SST_reg = xr.Dataset(); SST_reg['slope'] = P1_SST_slope; SST_reg['p_value'] = P1_SST_p_value; clim = data_domain(clim_SST, domain)
cb = reg_map(ax, SST_reg, domain, cmap, levs, transform_map); cb.ax.set_title('(°C)', pad=13, fontsize=11)
KEF  = grad_lat(clim.sel(lat=slice(25,50), lon=slice(142,175))).rolling(lon=8, center=True).mean()
ax.plot(KEF.lon, KEF, c='lime', linewidth=2.5, transform=transform_map)
ccrs_grid(ax, np.arange(120,270.1,30), np.arange(-20,60.1,20), 12); ccrs_plot(ax, 140, 175, 31, 45)
#ax.text(-0.13, 0.5, '1977-1999', va='center', ha='right', transform=ax.transAxes, weight='bold', fontsize=15, rotation=90)

# Fig. 1e
domain = [-22, 62, 100, 290]; cmap = cmap_white_center(plt.cm.RdBu_r); levs = np.arange(-24,24.1,2.4)
ax = plt.subplot(gs[1,2:4], projection=projection_map); ax.set_title('d', weight='bold', fontsize=20, x=-0.07, y=1.07)
ax.set_title('[1977-1999] LHFA', loc='left', fontsize=13)#; ax.set_title('D(0)J(+1)', loc='right', fontsize=13)
LHF_reg = xr.Dataset(); LHF_reg['slope'] = P1_LHF_slope; LHF_reg['p_value'] = P1_LHF_p_value
cb = reg_map(ax, LHF_reg, domain, cmap, levs, transform_map); cb.ax.set_title(r'(W m$^{-2}$)', pad=13, fontsize=11)
ax.plot(KEF.lon, KEF, c='lime', linewidth=2.5, transform=transform_map)
ccrs_grid(ax, np.arange(120,270.1,30), np.arange(-20,60.1,20), 12); ccrs_plot(ax, 140, 175, 31, 45)

syear = '2000'; eyear = '2022'
# Fig. 1d
domain = [-22, 62, 100, 290]; cmap = cmap_white_center(plt.cm.RdBu_r); levs = np.arange(-0.6,0.61,0.06)
ax = plt.subplot(gs[2,0:2], projection=projection_map); ax.set_title('e', weight='bold', fontsize=20, x=-0.07, y=1.07)
ax.set_title('[2000-2022] SSTA', loc='left', fontsize=13)#; ax.set_title('D(0)J(+1)', loc='right', fontsize=13)
SST_reg = xr.Dataset(); SST_reg['slope'] = P2_SST_slope; SST_reg['p_value'] = P2_SST_p_value; clim = data_domain(clim_SST, domain)
cb = reg_map(ax, SST_reg, domain, cmap, levs, transform_map); cb.ax.set_title('(°C)', pad=13, fontsize=11)
KEF  = grad_lat(clim.sel(lat=slice(25,50), lon=slice(142,175))).rolling(lon=8, center=True).mean()
ax.plot(KEF.lon, KEF, c='lime', linewidth=2.5, transform=transform_map)
ccrs_grid(ax, np.arange(120,270.1,30), np.arange(-20,60.1,20), 12); ccrs_plot(ax, 140, 175, 31, 45)
#ax.text(-0.13, 0.5, '2000-2022', va='center', ha='right', transform=ax.transAxes, weight='bold', fontsize=15, rotation=90)

# Fig. 1f
domain = [-22, 62, 100, 290]; cmap = cmap_white_center(plt.cm.RdBu_r); levs = np.arange(-24,24.1,2.4)
ax = plt.subplot(gs[2,2:4], projection=projection_map); ax.set_title('f', weight='bold', fontsize=20, x=-0.07, y=1.07)
ax.set_title('[2000-2022] LHFA', loc='left', fontsize=13)#; ax.set_title('D(0)J(+1)', loc='right', fontsize=13)
LHF_reg = xr.Dataset(); LHF_reg['slope'] = P2_LHF_slope; LHF_reg['p_value'] = P2_LHF_p_value
cb = reg_map(ax, LHF_reg, domain, cmap, levs, transform_map); cb.ax.set_title(r'(W m$^{-2}$)', pad=13, fontsize=11)
ax.plot(KEF.lon, KEF, c='lime', linewidth=2.5, transform=transform_map)
ccrs_grid(ax, np.arange(120,270.1,30), np.arange(-20,60.1,20), 12); ccrs_plot(ax, 140, 175, 31, 45)

plt.tight_layout(w_pad=3, h_pad=2)
plt.show()
