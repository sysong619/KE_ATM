import numpy as np
import xarray as xr

def lonflip(var, option):
    if option == 1:  # 0 to 360 -> -180 to 180
        var['lon'] = xr.where(var['lon'] > 180, var['lon'] - 360, var['lon'])
        var = var.sortby(var['lon'])  # 경도에 따라 정렬

    elif option == 2:  # -180 to 180 -> 0 to 360
        var['lon'] = xr.where(var['lon'] < 0, var['lon'] + 360, var['lon'])
        var = var.sortby(var['lon'])  # 경도에 따라 정렬

    elif option == 3:  # -180 to 180 -> -360 to 0
        var['lon'] = xr.where(var['lon'] > 0, var['lon'] - 360, var['lon'])
        var = var.sortby('lon')  # 경도에 따라 정렬

    elif option == 4: # 0 to 360 -> -360 to 0
        var['lon'] = xr.where(var['lon'] > 0, var['lon'] - 360, var['lon'])
        var = var.sortby(var['lon'])  # 경도에 따라 정렬

    return var
