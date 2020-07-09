#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 12:24:51 2020

@author: noahbrauer
"""
r"""
This script uses ERA-5 renalaysis data to plot relative vorticity, geopotential height, and wind at 500 mb. When selecting data, you will need to choose a level,
geopotential height, u-wind, v-wind, and relative vorticity. 
"""



#Import all libraries 

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

file = '500.nc'
nc = Dataset(file ,'r')

#Extract the attributes from the file; This will vary depending on what you are plotting

latitude = nc.variables['latitude'][:]
longitude = nc.variables['longitude'][:]
relative_vort = nc.variables['vo'][:]*10**5
z = nc.variables['z'][:]/100 #Convert to decameters
u = nc.variables['u'][:]*1.94384 #Convert from m/s to knots
v = nc.variables['v'][:]*1.94384 #Convert fromn m/s to knots


#Create a grid for plotting
lat2,lon2 = np.meshgrid(latitude,longitude)

#Mask out values of negative (anticyclonic) vorticity

relative_vort[relative_vort<=0] = np.nan


def plot_500_vort(time,lon_min,lon_max,lat_min,lat_max,title_font_size,declutter = None):
    
    r"""
    This function plots relative vorticity, gepotential height, and wind (knots) on a grid.
    
    Parameters:
    -----------
    time(int): Time index for the time being plotted
    lon_min,lon_max(float): Minimum and maximum longitude of the grid domain
    lat_min,lat_max(float): Minimum and maximum latitude of the grid domain
    title_font_size(float): Font size of the title and colorbar label
    declutter(int): Sets the declutter rate of the wind barbs; Greater number means lower barb density; Default is 12
    
    """


    if declutter is None:
        declutter = 12
        
        
    title_name = '500 mb Relative Vorticity, Geopotential height (dam), Winds (kt) '
    time = input('Enter the analysis time:')    

    cmin = 0; cmax = 26; cint = 1; clevs = np.round(np.arange(cmin,cmax,cint),2)
   
    plt.figure(figsize=(10,10))
    
    xlim = np.array([lon_min,lon_max]); ylim = np.array([lat_min,lat_max])

    m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
    m.drawcoastlines(); m.drawstates(), m.drawcountries()  
    cs = m.contourf(lon2,lat2,relative_vort[time,:,:].T,clevs,cmap='YlOrRd',extend='both') 
    cs2 = m.contour(lon2,lat2,z[time,:,:].T, colors = 'k')
    plt.clabel(cs2, fontsize=10, inline=1,fmt = '%1.0f')
    plt.barbs(lon2[::declutter,::declutter],lat2[::declutter,::declutter],u[time,::declutter,::declutter].T,v[time,::declutter,::declutter].T)
    m.drawcounties()

    cbar = m.colorbar(cs,size='2%')
    cbar.ax.set_ylabel(r'[$10^{-5}$ $s^{-1}$]',size=title_font_size)
    plt.title(str(title_name) + str(time),name='Calibri',size=title_font_size)
    
    plot = plt.show(block=False) 
    
    return plot



plot_vort = plot_500_vort(0,-120,-75,20,65,20)
