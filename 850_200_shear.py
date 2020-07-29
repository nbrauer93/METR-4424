#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 15:46:46 2020

@author: noahbrauer
"""

import matplotlib.pyplot as plt
from netCDF4 import Dataset, num2date, MFDataset
import netCDF4 as nc
import numpy as np
from datetime import datetime
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
from matplotlib.ticker import MultipleLocator




import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap


#Read in the file

file = 'shear.nc'
nc = Dataset(file ,'r')

#Extract the attributes from the file; This will vary depending on what you are plotting

latitude = nc.variables['latitude'][:]
longitude = nc.variables['longitude'][:]
level = nc.variables['level'][:]
u = nc.variables['u'][:]*1.94384 #Convert from m/s to knots
v = nc.variables['v'][:]*1.94384 #Convert fromn m/s to knots


u_200 = nc.variables['u'][0,0,:,:]*1.94384
v_200 = nc.variables['v'][0,0,:,:]*1.94384

u_850 = nc.variables['u'][0,1,:,:]*1.94384
v_850 = nc.variables['v'][0,1,:,:]*1.94384
#Calculate the magnitude of the wind
#%%
def wind_magnitude(u,v):
    
    magnitude = np.sqrt((u**2)+(v**2))
    
    return magnitude

wind_850 = wind_magnitude(u_850,v_850)
wind_200 = wind_magnitude(u_200,v_200)

wind_magnitude = wind_magnitude(u,v)

def wind_shear(wind_upper,wind_lower):
    shear = wind_upper-wind_lower
    
    return shear


shear_mag = wind_shear(wind_200,wind_850)



#Create a grid for plotting
lat2,lon2 = np.meshgrid(latitude,longitude)


#Mask out values of wind less than 50 knots

wind_magnitude[wind_magnitude<50] = np.nan
shear_mag[shear_mag<5] = np.nan

#%%

def plot_streamlines(lon_min,lon_max,lat_min,lat_max,min_value, max_value, value_interval, title_font_size,density,declutter = None):



    if declutter is None:
        declutter = 12
        
        
    title_name = '850-200 mb Shear, 200 mb Streamlines '
    time = input('Enter the analysis time:')    

    cmin = min_value; cmax = max_value; cint = value_interval; clevs = np.round(np.arange(cmin,cmax,cint),2)
   
    plt.figure(figsize=(10,10))
    
    xlim = np.array([lon_min,lon_max]); ylim = np.array([lat_min,lat_max])

    m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
    m.drawcoastlines(); m.drawstates(), m.drawcountries() 
    cs = m.contourf(lon2,lat2,shear_mag.T, clevs, cmap = 'YlOrRd', extend = 'max')
    cs2 = plt.streamplot(longitude,latitude,u_200,v_200, density = density, linewidth = 2, color = 'k')
    
    m.drawcounties()

    cbar = m.colorbar(cs,size='2%')
    cbar.ax.set_ylabel('[knots]',size=title_font_size)
    plt.title(str(title_name) + str(time),name='Calibri',size=title_font_size)
    
    plot = plt.show(block=False) 
    
    return plot



plot_wind = plot_streamlines(-100,-5,-20,65,5,100, 5, 20,2)
