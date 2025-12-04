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

def data_map(ax, data, domain, cmap, levs, transform_map):
    clim = data_domain(data, domain); cf = ccrs_contourf(ax, clim, levs, cmap, transform_map)
    cb = plt.colorbar(cf, shrink=0.7, pad=0.07, orientation='vertical', ticks=levs[::5])
    ax.set_aspect('auto')
    return cb

def comp_map(ax, comp, pval, domain, cmap, levs, transform_map, siglev=0.1):
    comp = data_domain(comp, domain); cf = ccrs_contourf(ax, comp, levs, cmap, transform_map)
    cb = plt.colorbar(cf, shrink=0.7, pad=0.07, orientation='vertical', ticks=levs[::5])
    p_val = data_domain(pval, domain); ccrs_sig(p_val, transform_map, siglev)
    ax.set_aspect('auto')
    return cb
  
def sign_agreement(comp, threshold=8):
    nens, ny, nx = comp.shape

    sign_map = np.sign(comp)  # (ensemble, lat, lon)
    pos_count = (sign_map > 0).sum(axis=0)  # (lat, lon)
    neg_count = (sign_map < 0).sum(axis=0)  # (lat, lon)
  
    agreement = np.full((ny, nx), 1.0)
    agreement[np.where((pos_count >= threshold) | (neg_count >= threshold))] = 0.09
    agreement = np2xr(comp[0], agreement)

    return agreement

def ccrs_sig(sig, transform, siglev, hatch='////', color='k'):
    sig = plt.contourf(sig.lon, sig.lat, sig, levels=[0., siglev, 1.0], colors='none', transform=transform)
    for i, contour in enumerate(sig.collections):
        if i == 0: contour.set_hatch(hatch); contour.set_edgecolor(color); contour.set_linewidth(0.0)

def quiver_comp(ax, ucomp, upval, vcomp, vpval, domain, transform_map, siglev=0.1, step=4, scale=40):
    u_data = data_domain(ucomp, domain); u_pval = data_domain(upval, domain)
    v_data = data_domain(vcomp, domain); v_pval = data_domain(vpval, domain)
    u_sig = u_data.where(((u_pval < siglev) |  (v_pval < siglev)))
    v_sig = v_data.where(((u_pval < siglev) |  (v_pval < siglev)))
    ccrs_quiver(ax, u_data, v_data, 'silver', transform_map, step=step, scale=scale)
    quiver = ccrs_quiver(ax, u_sig,  v_sig,  'k', transform_map, step=step, scale=scale)
    ax.set_aspect('auto')
    return quiver

####################################################################################################################################
ds = xr.open_dataset("../process/Fig03_GOGA_STA_clim.nc")
clim_STA = ds['clim']

ds = xr.open_dataset("../process/Fig03_GOGA_KE_reg_STA_P1.nc")
P1_STA_slope = ds['STA_slope']; P1_STA_p_value = ds['STA_p_value']
ds = xr.open_dataset("../process/Fig03_GOGA_KE_reg_SLPUV_P1.nc")
P1_SLP_slope = ds['SLP_slope']; P1_SLP_p_value = ds['SLP_p_value']
P1_UWD_slope = ds['UWD_slope']; P1_UWD_p_value = ds['UWD_p_value']
P1_VWD_slope = ds['VWD_slope']; P1_VWD_p_value = ds['VWD_p_value']

ds = xr.open_dataset("../process/Fig03_GOGA_KE_reg_STA_P2.nc")
P2_STA_slope = ds['STA_slope']; P2_STA_p_value = ds['STA_p_value']
ds = xr.open_dataset("../process/Fig03_GOGA_KE_reg_SLPUV_P2.nc")
P2_SLP_slope = ds['SLP_slope']; P2_SLP_p_value = ds['SLP_p_value']
P2_UWD_slope = ds['UWD_slope']; P2_UWD_p_value = ds['UWD_p_value']
P2_VWD_slope = ds['VWD_slope']; P2_VWD_p_value = ds['VWD_p_value']

####################################################################################################################################
fig = plt.figure(figsize=(13,7)); gs = gridspec.GridSpec(nrows=2, ncols=2, height_ratios=[1,1], width_ratios=[1,1])

domain = [18, 62, 110, 250]
cmap = cmap_white_center(plt.cm.PiYG_r); levs = np.arange(-1.2,1.21,0.12)
ax = plt.subplot(gs[0,0], projection=projection_map); ax.set_title('a', weight='bold', fontsize=20, x=-0.07, y=1.07)
ax.set_title('[1977-1999] STKA', loc='left', fontsize=13)#; ax.set_title('D(0)J(+1)', loc='right', fontsize=13)
slope = P1_STA_slope; pval = sign_agreement(slope); clim = data_domain(clim_STA, domain)
cb = comp_map(ax, slope.mean(axis=0), pval, domain, cmap, levs, transform_map); cb.ax.set_title(r'(K m s$^{-1}$)', pad=13, fontsize=11)
cs = ax.contour(clim.lon, clim.lat, clim, levels=[4,8,12], colors='grey', linewidths=1.3, transform=transform_map)
ax.clabel(cs, inline=True, fontsize=13, fmt='%d', colors='k')
STL  = clim.lat[clim.argmax(axis=0)].sel(lon=slice(140,238)).rolling(lon=4, center=True).mean()
ax.plot(STL.lon, STL, c='k', linewidth=2.4, transform=transform_map)
ccrs_grid(ax, np.arange(120,240.1,30), np.arange(20,60.1,10), 12)#; ccrs_plot(ax, 140, 175, 31 ,45)
ax.text(-0.13, 0.5, 'AGCM GOGA', va='center', ha='right', transform=ax.transAxes, weight='bold', fontsize=13, rotation=90)

cmap = cmap_white_center(plt.cm.coolwarm); levs = np.arange(-2,2.1,0.2)
slope  = P1_SLP_slope; pval  = sign_agreement(slope)
slopeu = P1_UWD_slope; pvalu = sign_agreement(slopeu)
slopev = P1_VWD_slope; pvalv = sign_agreement(slopev)
ax = plt.subplot(gs[0,1], projection=projection_map); ax.set_title('b', weight='bold', fontsize=20, x=-0.07, y=1.07)
ax.set_title('[1977-1999] SLPA/UV$_{925}$A', loc='left', fontsize=13)#; ax.set_title('D(0)J(+1)', loc='right', fontsize=13)
cb = comp_map(ax, slope.mean(axis=0), pval, domain, cmap, levs, transform_map); cb.ax.set_title('(hPa)', pad=13, fontsize=11)
qv = quiver_comp(ax, slopeu.mean(axis=0), pvalu, slopev.mean(axis=0), pvalv, domain, transform_map, scale=6, step=5)
ax.quiverkey(qv, X=0.85, Y=1.05, U=0.5, label=r'0.5 m s$^{-1}$', labelpos='E')
ccrs_grid(ax, np.arange(120,240.1,30), np.arange(20,60.1,10), 12); ccrs_plot(ax, 140, 175, 31 ,45)

syear = '2000'; eyear = '2021'
cmap = cmap_white_center(plt.cm.PiYG_r); levs = np.arange(-1.2,1.21,0.12)
ax = plt.subplot(gs[1,0], projection=projection_map); ax.set_title('c', weight='bold', fontsize=20, x=-0.07, y=1.07)
ax.set_title('[2000-2021] STKA', loc='left', fontsize=13)#; ax.set_title('D(0)J(+1)', loc='right', fontsize=13)
slope = P2_STA_slope; pval = sign_agreement(slope); clim = data_domain(clim_STA, domain)
cb = comp_map(ax, slope.mean(axis=0), pval, domain, cmap, levs, transform_map); cb.ax.set_title(r'(K m s$^{-1}$)', pad=13, fontsize=11)
cs = ax.contour(clim.lon, clim.lat, clim, levels=[4,8,12], colors='grey', linewidths=1.3, transform=transform_map)
ax.clabel(cs, inline=True, fontsize=13, fmt='%d', colors='k')
STL  = clim.lat[clim.argmax(axis=0)].sel(lon=slice(140,238)).rolling(lon=4, center=True).mean()
ax.plot(STL.lon, STL, c='k', linewidth=2.4, transform=transform_map)
ccrs_grid(ax, np.arange(120,240.1,30), np.arange(20,60.1,10), 12)#; ccrs_plot(ax, 140, 175, 31 ,45)
ax.text(-0.13, 0.5, 'AGCM GOGA', va='center', ha='right', transform=ax.transAxes, weight='bold', fontsize=13, rotation=90)

cmap = cmap_white_center(plt.cm.coolwarm); levs = np.arange(-2,2.1,0.2)
ax = plt.subplot(gs[1,1], projection=projection_map); ax.set_title('d', weight='bold', fontsize=20, x=-0.07, y=1.07)
ax.set_title('[2000-2021] SLPA/UV$_{925}$A', loc='left', fontsize=13)#; ax.set_title('D(0)J(+1)', loc='right', fontsize=13)
slope  = P2_SLP_slope; pval  = sign_agreement(slope)
slopeu = P2_UWD_slope; pvalu = sign_agreement(slopeu)
slopev = P2_VWD_slope; pvalv = sign_agreement(slopev)
cb = comp_map(ax, slope.mean(axis=0), pval, domain, cmap, levs, transform_map); cb.ax.set_title('(hPa)', pad=13, fontsize=11)
qv = quiver_comp(ax, slopeu.mean(axis=0), pvalu, slopev.mean(axis=0), pvalv, domain, transform_map, scale=8, step=5)
ax.quiverkey(qv, X=0.85, Y=1.05, U=0.5, label=r'0.5 m s$^{-1}$', labelpos='E')
ccrs_grid(ax, np.arange(120,240.1,30), np.arange(20,60.1,10), 12); ccrs_plot(ax, 140, 175, 31 ,45)

plt.tight_layout(w_pad=3, h_pad=3)
plt.show()
