import numpy as np
import pandas as pd
import xarray as xr
import sacpy as scp

import re, os, sys, glob
sys.path.append('../library')
from cbasic import *
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
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12

def data_map(ax, data, domain, cmap, levs, transform_map):
    clim = data_domain(data, domain); cf = ccrs_contourf(ax, clim, levs, cmap, transform_map)
    #cb = plt.colorbar(cf, shrink=0.7, pad=0.07, orientation='vertical', ticks=levs[::5])
    ax.set_aspect('auto')
    return cf

def comp_map(ax, comp, pval, domain, cmap, levs, transform_map, siglev=0.1):
    comp = data_domain(comp, domain); cf = ccrs_contourf(ax, comp, levs, cmap, transform_map)
    #cb = plt.colorbar(cf, shrink=0.6, pad=0.1, orientation='horizontal', ticks=levs[::5])
    p_val = data_domain(pval, domain); ccrs_sig(p_val, transform_map, siglev)
    ax.set_aspect('auto')
    return cf

def ccrs_sig(sig, transform, siglev, hatch='////', color='k'):
    sig = plt.contourf(sig.lon, sig.lat, sig, levels=[0., siglev, 1.0], colors='none', transform=transform)
    for i, contour in enumerate(sig.collections):
        if i == 0: contour.set_hatch(hatch); contour.set_edgecolor(color); contour.set_linewidth(0.0)

####################################################################################################################################
ds = xr.open_dataset("../process/Fig05_E3SM_CLIM_STA_clim_P2.nc")
clim_STA_P2 = ds['clim']

ds = xr.open_dataset("../process/Fig05_E3SM_KE_TP_forcing_P2.nc")
SST_comp_TP_P2 = ds['SST_comp']; SST_pval_TP_P2 = ds['SST_p_value']
ds = xr.open_dataset("../process/Fig05_E3SM_KE_TP_comp_P2.nc")
STA_comp_TP_P2 = ds['STA_comp']; STA_pval_TP_P2 = ds['STA_p_value']
SLP_comp_TP_P2 = ds['SLP_comp']; SLP_pval_TP_P2 = ds['SLP_p_value']
UWD_comp_TP_P2 = ds['UWD_comp']; UWD_pval_TP_P2 = ds['UWD_p_value']
VWD_comp_TP_P2 = ds['VWD_comp']; VWD_pval_TP_P2 = ds['VWD_p_value']

ds = xr.open_dataset("../process/Fig05_E3SM_KE_NP_forcing_P2.nc")
SST_comp_NP_P2 = ds['SST_comp']; SST_pval_NP_P2 = ds['SST_p_value']
ds = xr.open_dataset("../process/Fig05_E3SM_KE_NP_comp_P2.nc")
STA_comp_NP_P2 = ds['STA_comp']; STA_pval_NP_P2 = ds['STA_p_value']
SLP_comp_NP_P2 = ds['SLP_comp']; SLP_pval_NP_P2 = ds['SLP_p_value']
UWD_comp_NP_P2 = ds['UWD_comp']; UWD_pval_NP_P2 = ds['UWD_p_value']
VWD_comp_NP_P2 = ds['VWD_comp']; VWD_pval_NP_P2 = ds['VWD_p_value']

ds = xr.open_dataset("../process/Fig05_E3SM_KE_TPNP_forcing_P2.nc")
SST_comp_TPNP_P2 = ds['SST_comp']; SST_pval_TPNP_P2 = ds['SST_p_value']
ds = xr.open_dataset("../process/Fig05_E3SM_KE_TPNP_comp_P2.nc")
STA_comp_TPNP_P2 = ds['STA_comp']; STA_pval_TPNP_P2 = ds['STA_p_value']
SLP_comp_TPNP_P2 = ds['SLP_comp']; SLP_pval_TPNP_P2 = ds['SLP_p_value']
UWD_comp_TPNP_P2 = ds['UWD_comp']; UWD_pval_TPNP_P2 = ds['UWD_p_value']
VWD_comp_TPNP_P2 = ds['VWD_comp']; VWD_pval_TPNP_P2 = ds['VWD_p_value']

####################################################################################################################################
fig = plt.figure(figsize=(16,10)); gs = gridspec.GridSpec(nrows=3, ncols=4, height_ratios=[1,1.1,1.1], width_ratios=[1,1,1,0.01])

projection_map = ccrs.Robinson(central_longitude=180); domain = [-35, 65, 0, 360]
cmap = cmap_white_center(plt.cm.RdBu_r); levs = np.arange(-2,2.1,0.2)

ax = plt.subplot(gs[0,0], projection=projection_map); ax.text(-0.09, 1.1, "a", fontsize=20, weight='bold', transform=ax.transAxes)
ax.set_title('Tropical Pacific minus CTRL\n\n', fontsize=14, x=0.5, weight='bold'); ax.set_title('          SST forcing from P2 (2000-2022)', loc='left', fontsize=13)
comp = SST_comp_TP_P2; pval = SST_pval_TP_P2
cf = comp_map(ax, comp.where(comp!=0), pval, domain, cmap, levs, transform_map)
ccrs_grid(ax, np.arange(20,360.1,80), np.arange(-30,60.1,30), 13)

ax = plt.subplot(gs[0,1], projection=projection_map); ax.text(-0.09, 1.1, "b", fontsize=20, weight='bold', transform=ax.transAxes)
ax.set_title('North Pacific minus CTRL\n\n', fontsize=14, x=0.5, weight='bold'); ax.set_title('          SST forcing from P2 (2000-2022)', loc='left', fontsize=13)
comp = SST_comp_NP_P2; pval = SST_pval_NP_P2
cf = comp_map(ax, comp.where(comp!=0), pval, domain, cmap, levs, transform_map)
ccrs_grid(ax, np.arange(20,360.1,80), np.arange(-30,60.1,30), 13)

ax = plt.subplot(gs[0,2], projection=projection_map); ax.text(-0.09, 1.1, "c", fontsize=20, weight='bold', transform=ax.transAxes)
ax.set_title('Tropical-North Pacific minus CTRL\n\n', fontsize=15, x=0.5, weight='bold'); ax.set_title('          SST forcing from P2 (2000-2022)', loc='left', fontsize=13)
comp = SST_comp_TPNP_P2; pval = SST_pval_TPNP_P2
cf = comp_map(ax, comp.where(comp!=0), pval, domain, cmap, levs, transform_map)
ccrs_grid(ax, np.arange(20,360.1,80), np.arange(-30,60.1,30), 13)

cax = fig.add_axes([0.95, 0.70, 0.006, 0.21])
cb = fig.colorbar(cf, cax=cax, orientation='vertical', ticks=levs[::5])
cb.set_label('(°C)', fontsize=13)

projection_map = ccrs.PlateCarree(central_longitude=180); domain = [18, 62, 110, 250]; bar = 60
cmap = cmap_white_center(plt.cm.PiYG_r); levs = np.arange(-2.4,2.41,0.24)*0.5

ax = plt.subplot(gs[1,0], projection=projection_map); ax.set_title('d', weight='bold', fontsize=20, x=-0.07, y=1.07)
ax.set_title('Tropical Pacific SST-forced ΔSTK', loc='left', fontsize=13)
comp = STA_comp_TP_P2; pval = STA_pval_TP_P2; clim = data_domain(clim_STA_P2, domain)
STL  = clim.lat[clim.argmax(axis=0)].sel(lon=slice(140,238)).rolling(lon=4, center=True).mean()
cf = comp_map(ax, comp.where(comp!=0)/2., pval, domain, cmap, levs, transform_map)
cs = ax.contour(clim.lon, clim.lat, clim, colors='grey', levels=[4,8,12], linewidths=1.3, transform=transform_map)
ax.clabel(cs, inline=True, fontsize=13, fmt='%d', colors='k')
ax.plot(STL.lon, STL, c='k', linewidth=2.4, transform=transform_map)
ccrs_grid(ax, np.arange(120,240.1,30), np.arange(20,60.1,10), 12)

ax = plt.subplot(gs[1,1], projection=projection_map); ax.set_title('e', weight='bold', fontsize=20, x=-0.07, y=1.07)
ax.set_title('North Pacific SST-forced ΔSTK', loc='left', fontsize=13)
comp = STA_comp_NP_P2; pval = STA_pval_NP_P2; clim = data_domain(clim_STA_P2, domain)
STL  = clim.lat[clim.argmax(axis=0)].sel(lon=slice(140,238)).rolling(lon=4, center=True).mean()
cf = comp_map(ax, comp.where(comp!=0)/2., pval, domain, cmap, levs, transform_map)
cs = ax.contour(clim.lon, clim.lat, clim, colors='grey', levels=[4,8,12], linewidths=1.3, transform=transform_map)
ax.clabel(cs, inline=True, fontsize=13, fmt='%d', colors='k')
ax.plot(STL.lon, STL, c='k', linewidth=2.4, transform=transform_map)
ccrs_grid(ax, np.arange(120,240.1,30), np.arange(20,60.1,10), 12)

ax = plt.subplot(gs[1,2], projection=projection_map); ax.set_title('f', weight='bold', fontsize=20, x=-0.07, y=1.07)
ax.set_title('Tropical-North Pacific SST-forced ΔSTK', loc='left', fontsize=13)
comp = STA_comp_TPNP_P2; pval = STA_pval_TPNP_P2; clim = data_domain(clim_STA_P2, domain)
STL  = clim.lat[clim.argmax(axis=0)].sel(lon=slice(140,238)).rolling(lon=4, center=True).mean()
cf = comp_map(ax, comp.where(comp!=0)/2., pval, domain, cmap, levs, transform_map)
cs = ax.contour(clim.lon, clim.lat, clim, colors='grey', levels=[4,8,12], linewidths=1.3, transform=transform_map)
ax.clabel(cs, inline=True, fontsize=13, fmt='%d', colors='k')
ax.plot(STL.lon, STL, c='k', linewidth=2.4, transform=transform_map)
ccrs_grid(ax, np.arange(120,240.1,30), np.arange(20,60.1,10), 12)

cax = fig.add_axes([0.95, 0.365, 0.006, 0.23])
cb = fig.colorbar(cf, cax=cax, orientation='vertical', ticks=levs[::5])
cb.set_label(r'(K m s$^{-1}$)', fontsize=13)

cmap = cmap_white_center(plt.cm.coolwarm); levs = np.arange(-4,4.1,0.4)*0.5

ax = plt.subplot(gs[2,0], projection=projection_map); ax.set_title('g', weight='bold', fontsize=20, x=-0.07, y=1.07)
ax.set_title('Tropical Pacific SST-forced ΔSLP/ΔUV$_{925}$', loc='left', fontsize=13)
comp  = SLP_comp_TP_P2;  pval = SLP_pval_TP_P2
compu = UWD_comp_TP_P2; pvalu = UWD_pval_TP_P2
compv = VWD_comp_TP_P2; pvalv = VWD_pval_TP_P2
cf = comp_map(ax, comp.where(comp!=0)/2., pval, domain, cmap, levs, transform_map)
qv = quiver_comp(ax, compu/2., pvalu, compv/2., pvalv, domain, transform_map, step=5, scale=13, siglev=0.1)
ccrs_grid(ax, np.arange(120,240.1,30), np.arange(20,60.1,10), 12); ccrs_plot(ax, 140, 175, 31, 45)

ax = plt.subplot(gs[2,1], projection=projection_map); ax.set_title('h', weight='bold', fontsize=20, x=-0.07, y=1.07)
ax.set_title('North Pacific SST-forced ΔSLP/ΔUV$_{925}$', loc='left', fontsize=13)
comp  = SLP_comp_NP_P2;  pval = SLP_pval_NP_P2
compu = UWD_comp_NP_P2; pvalu = UWD_pval_NP_P2
compv = VWD_comp_NP_P2; pvalv = VWD_pval_NP_P2
cf = comp_map(ax, comp.where(comp!=0)/2., pval, domain, cmap, levs, transform_map)
qv = quiver_comp(ax, compu/2., pvalu, compv/2., pvalv, domain, transform_map, step=5, scale=13, siglev=0.1)
ccrs_grid(ax, np.arange(120,240.1,30), np.arange(20,60.1,10), 12); ccrs_plot(ax, 140, 175, 31, 45)

ax = plt.subplot(gs[2,2], projection=projection_map); ax.set_title('i', weight='bold', fontsize=20, x=-0.07, y=1.07)
ax.set_title('Tropical-North Pacific SST-forced ΔSLP/ΔUV$_{925}$', loc='left', fontsize=13)
comp  = SLP_comp_TPNP_P2;  pval = SLP_pval_TPNP_P2
compu = UWD_comp_TPNP_P2; pvalu = UWD_pval_TPNP_P2
compv = VWD_comp_TPNP_P2; pvalv = VWD_pval_TPNP_P2
cf = comp_map(ax, comp.where(comp!=0)/2., pval, domain, cmap, levs, transform_map)
qv = quiver_comp(ax, compu/2., pvalu, compv/2., pvalv, domain, transform_map, step=5, scale=13, siglev=0.1)
ccrs_grid(ax, np.arange(120,240.1,30), np.arange(20,60.1,10), 12); ccrs_plot(ax, 140, 175, 31, 45)

cax = fig.add_axes([0.95, 0.033, 0.006, 0.23])
cb = fig.colorbar(cf, cax=cax, orientation='vertical', ticks=levs[::5])
cb.set_label('(hPa)', fontsize=13)
ax.quiverkey(qv, X=1.135, Y=1.115, U=1, label=r'1 m s$^{-1}$', labelpos='S')

plt.tight_layout(w_pad=5, h_pad=3)
plt.show()
