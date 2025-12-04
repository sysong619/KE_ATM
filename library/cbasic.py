#!/usr/bin/env python
# coding: utf-8

## Modules for Calculate netCDF
import numpy    as np
import xarray   as xr
import pandas   as pd
import sacpy    as scp
import netCDF4

## Modules for caculating statistics
from scipy   import stats, signal
from sklearn import linear_model
from eofs.standard import Eof

## Modules for plottings
import matplotlib.pyplot as plt
import matplotlib.colorbar as cb

EARTH_RADIUS = 6371000.0  # m

def np2xr(xr_array,np_array):
    n_dims = xr_array.shape
    tmp = np_array.shape

    if (n_dims == tmp):
        output = xr.DataArray( np_array,
                               dims = xr_array.dims,
                               coords = xr_array.coords, attrs = xr_array.attrs )
    else:
        print( 'Numpy array shape does not math with the Xarray' )

    return output

def norm(data):
    norm = (data - data.mean(dim='time'))/data.std(dim='time')
    return norm

def std(data):
    std = (data)/data.std(dim='time')
    return std

def season(data, months):
    return scp.XrTools.spec_moth_yrmean(data, months)

def calc_anomaly(var):
    clim_var = var.groupby("time.month").mean()
    anom_var = var.groupby("time.month") - clim_var
    dtr_anom_var = dtr(anom_var)
    return clim_var, anom_var, dtr_anom_var

def calc_anomaly_dclim(var1, var2):
    clim_var = var1.groupby("time.month").mean()
    anom_var = var2.groupby("time.month") - clim_var
    dtr_anom_var = dtr(anom_var)
    return clim_var, anom_var, dtr_anom_var

def extract_index(var,latS,latN,lonL,lonR):
    var_sel = var.sel(lat=slice(latS,latN),lon=slice(lonL,lonR))
    lat = var_sel['lat']; lon = var_sel['lon']
    weights_lat      = np.cos(np.deg2rad(lat))
    weights_lat.name = "weights"
    index_var = var_sel.weighted(weights_lat).mean(("lat","lon"))
    return index_var
