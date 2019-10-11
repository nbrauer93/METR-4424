#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 13:58:19 2019

@author: noahbrauer
"""

############# Isentropic coordinate function is courtesy of MetPy #############


import cartopy.crs as ccrs
import cartopy.feature as cfeature

import matplotlib.pyplot as plt
from netCDF4 import Dataset, num2date, MFDataset
import numpy as np
from datetime import datetime
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
from matplotlib.ticker import MultipleLocator
import xarray as xr
import metpy.calc as mpcalc
from metpy.cbook import get_test_data
from metpy.plots import add_metpy_logo, add_timestamp
from metpy.units import units


#List the files; You will have to modify this...


hgt_file = 'hgt.201902.nc'
temp_file = 'air.201902.nc'
q_file = 'shum.201902.nc'
u_file = 'uwnd.201902.nc'
v_file = 'vwnd.201902.nc'

#Read in each file 

nc_hgt = Dataset(hgt_file, 'r')
nc_temp = Dataset(temp_file, 'r')
nc_q = Dataset(q_file, 'r')
nc_u = Dataset(u_file, 'r')
nc_v = Dataset(v_file, 'r')


#Read in file attributes and extract times (Feb 19,2019 is what I'm pulling here)

lat = nc_hgt.variables['lat'][:]
lon = nc_hgt.variables['lon'][:] 

narr = {}
time = nc_hgt.variables['time'][:]
timeUnits = nc_hgt.variables['time'].units
tmpDates = num2date(time,timeUnits,calendar='gregorian')
narr['date'] = np.asarray([datetime(d.year,d.month,d.day) for d in tmpDates])
narr['day'] = np.asarray([d.day for d in narr['date']])
narr['month'] = np.asarray([d.month for d in narr['date']])
narr['year'] = np.asarray([d.year for d in narr['date']])

#Assign time index for our time of interest; You will have to change this line. 

time_index = np.where(narr['day']==19)[0]

#Now read in meteorological data for this day
level = nc_hgt.variables['level'][:] #In hPa
z = nc_hgt.variables['hgt'][time_index,:,:,:] #In meters
temp = nc_temp.variables['air'][time_index,:,:,:] #In Kelvin
q = nc_q.variables['shum'][time_index,:,:,:]
u = nc_u.variables['uwnd'][time_index,:,:,:]
v = nc_v.variables['vwnd'][time_index,:,:,:]


#%%


#Select your time of choice (03 UTC on 2/19 is the default now; change the 1 to 2 if you want 06 UTC, etc.; times are in 3 hour increments)

temp = temp[1,:,:,:]
z = z[1,:,:,:]
q = q[1,:,:,:]
u = u[1,:,:,:]
v = v[1,:,:,:]

#Assign proper units to each variables
temp = temp* units.kelvin
level = level*units.hectopascal


#Define isentropic levels

isentlevs = [296.] * units.kelvin

#Now convert to isentropic coordinates

isentropic = mpcalc.isentropic_interpolation(isentlevs, level, temp, q, u, v, z, tmpk_out = True)



#Separate variables (don't worry, this has nothing to do with PDEs)

isentprs, isenttmp, isentspech, isentu, isentv, isenthgt = isentropic


#Convert m/s into knots for plotting


def ms_to_knot(wind):
    knots = wind/0.514
    return knots

u_knots = ms_to_knot(isentu)
v_knots = ms_to_knot(isentv)


#Calculate relative humidity from the specific humidity field; RH = spec_hum/saturation spec_hum * 100

isentrh = 100 * mpcalc.relative_humidity_from_specific_humidity(isentspech, isenttmp, isentprs)



#NaN all values where longitude>0; This is an issue with NARR data...


rh_nan = np.ones((1,277,349))*np.nan
u_knots_nan = np.ones((1,277,349))*np.nan
v_knots_nan = np.ones((1,277,349))*np.nan
isentprs_nan = np.ones((1,277,349))*np.nan

for i in range(rh_nan.shape[1]):
    for j in range(rh_nan.shape[2]):
        
        if lon[i,j]>=0:
            rh_nan[:,i,j] = np.nan
            u_knots_nan[:,i,j] = np.nan
            v_knots_nan[:,i,j] = np.nan
            isentprs_nan[:,i,j] = np.nan
            
        else:
            rh_nan[:,i,j] = isentrh[:,i,j]
            u_knots_nan[:,i,j] = u_knots[:,i,j]
            v_knots_nan[:,i,j] = v_knots[:,i,j]
            isentprs_nan[:,i,j] = isentprs[:,i,j]



    



#%% 
    
#Okay, now we can plot



# Set up our projection
crs = ccrs.LambertConformal(central_longitude=-100.0, central_latitude=45.0)

# Coordinates to limit map area
bounds = [(-120., -75., 25., 50.)]
# Choose a level to plot, in this case 296 K
level = 296

#Can tweak figure size if you would like
fig = plt.figure(figsize=(15., 8.))

ax = fig.add_subplot(1, 1, 1, projection=crs)
ax.set_extent(*bounds, crs=ccrs.PlateCarree())
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.75)
ax.add_feature(cfeature.STATES, linewidth=0.5)

# Plot the surface
clevisent = np.arange(0, 1000, 25)
cs = ax.contour(lon, lat, isentprs_nan[level, :, :], clevisent,
                colors='k', linewidths=1.0, linestyles='solid', transform=ccrs.PlateCarree())
ax.clabel(cs, fontsize=10, inline=1, inline_spacing=7,
          fmt='%i', rightside_up=True, use_clabeltext=True)

# Plot RH

cf = ax.contourf(lon, lat, rh_nan[level, :, :], range(10, 106, 5), cmap=plt.cm.gist_earth_r, transform=ccrs.PlateCarree())


#Add a colorbar
cb = fig.colorbar(cf, orientation='horizontal', extend='max', aspect=65, shrink=0.5, pad=0.05,
                  extendrect='True')
cb.set_label('Relative Humidity', size='x-large')


#Plot wind barbs on the isentropic surface
ax.barbs(lon, lat, u_knots_nan[level, :, :], v_knots_nan[level, :, :], length=6,regrid_shape=20, transform=ccrs.PlateCarree())
# Make a title
plt.title('296 K Pressure, Relative Humidity, Wind 2/19 03 UTC', size = 20)
fig.tight_layout()




