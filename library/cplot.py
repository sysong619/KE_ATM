#!/usr/bin/env python
# coding: utf-8

## Modules for Calculate netCDF
import cmocean
import numpy as np
import pandas as pd
import xarray as xr
import sacpy as scp

## Modules for plottings
from matplotlib import gridspec
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point
from matplotlib.ticker import FixedLocator

def data_domain(data, domain):
    latS = domain[0]; latN = domain[1]; lonL = domain[2]; lonR = domain[3]
    plot_data = data.sel(lat=slice(latS,latN), lon=slice(lonL,lonR))
    return plot_data

def cmap_white_center(colormap, c_range=12):
    cmap = colormap
    cmap_white = [cmap(i) for i in np.arange(cmap.N)]
    cen1 = int(cmap.N / 2); cen2 = int(cen1 - 1)

    for i in np.arange(cen2 - (c_range-1), cen1 + (c_range)):
        cmap_white[i] = (1.0, 1.0, 1.0, 1.0)
    cmap_white = mcolors.LinearSegmentedColormap.from_list('my_cmap', cmap_white, cmap.N)
    return cmap_white

def ccrs_plot(ax, xmin, xmax, ymin, ymax, color='k'):
    x = np.array([xmin, xmax, xmax, xmin, xmin]); y = np.array([ymin, ymin, ymax, ymax, ymin])
    ax.plot(x, y, c='w', linewidth=2, transform=ccrs.PlateCarree())
    ax.plot(x, y, c=color, linewidth=1, transform=ccrs.PlateCarree())

def ccrs_grid(ax, xloc, yloc, size, padx=3, pady=3):
    if np.max(xloc) > 180:
        xloc = (xloc + 180) % 360 -180
    gl = ax.gridlines(draw_labels=True, linewidth=0)
    gl.xlocator = FixedLocator(xloc); gl.ylocator = FixedLocator(yloc)
    gl.top_labels = False;   gl.xlabel_style = {'size': size}; gl.xpadding = padx
    gl.right_labels = False; gl.ylabel_style = {'size': size}; gl.ypadding = pady

def ccrs_contourf(ax, data, levs, cmap, transform, shrink=0.9, pad=0.1):
    cf = ax.contourf(data.lon, data.lat, data, levels=levs, cmap=cmap, extend='both', transform=transform)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    return cf

def ccrs_sig(sig, transform, siglev, hatch='///', color='k'):
    sig = plt.contourf(sig.lon, sig.lat, sig, levels=[0., siglev, 1.0], colors='none', transform=transform)
    for i, contour in enumerate(sig.collections):
        if i == 0: contour.set_hatch(hatch); contour.set_edgecolor(color); contour.set_linewidth(0.0)

def ccrs_quiver(ax, udata, vdata, color, transform_map, step=2, scale=360, width=0.006, fontsize=10, xset=1, yset=1):
    lon, lat = udata.lon, udata.lat
    qv = ax.quiver(
        lon[::step].values, lat[::step].values,
        udata[::step, ::step].values, vdata[::step, ::step].values,
        scale=scale, width=width, color=color, transform=transform_map
    )
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    return qv
          
