#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 15:21:25 2020

@author: noahbrauer
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from netCDF4 import Dataset, num2date, MFDataset
from datetime import datetime
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
from matplotlib.ticker import MultipleLocator
from scipy.ndimage import gaussian_filter



import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap
import pyart


file = 'pv.nc'
era = {}

nc = Dataset(file, 'r', unpack = True)

lat = nc.variables['latitude'][:]
lon = nc.variables['longitude'][:] - 360

def define_lat_lon_grid(lon_min,lon_max,lat_min,lat_max):
    
    r"""
    Function creates a lat-lon grid for plotting
    
    Parameters:
    -----------

    lon_min, lon_max(float): Minimum and maximum longitude coordinate on -180 to 180 degree grid
    lat_min, lat_max(float): Minimum and maximum latitude coordinate on -90 to 90 degree grid   
    
    
    """
    

    xlim = [lon_min,lon_max]; ylim = [lat_min,lat_max]
    ilat = np.where((lat>ylim[0])&(lat<ylim[1]))[0]
    ilon = np.where((lon>xlim[0])&(lon<xlim[1]))[0]

    latitude = lat[ilat]
    longitude = lon[ilon]


    lat2,lon2 = np.meshgrid(latitude,longitude)
    
    return ilat,ilon, latitude,longitude, lat2,lon2


grid_output = define_lat_lon_grid(-93,-85,36,41)

ilat = grid_output[0]
ilon = grid_output[1]
latitude = grid_output[2]
longitude = grid_output[3]
lat2 = grid_output[4]
lon2 = grid_output[5]


#%%
level = nc.variables['level'][:]*100
level_mb = level/100

time = nc.variables['time'][:]
timeUnits = nc.variables['time'].units
tmpDates = num2date(time,timeUnits,calendar='gregorian')
era['date'] = np.asarray([datetime(d.year,d.month,d.day) for d in tmpDates])
era['day'] = np.asarray([d.day for d in era['date']])
era['month'] = np.asarray([d.month for d in era['date']])
era['year'] = np.asarray([d.year for d in era['date']])


pv = nc.variables['pv'][:,:,ilat,ilon]
temp = nc.variables['t'][:,:,ilat,ilon]




#%%
#Calculate potential temperature from temperature

def potential_temp(temp,pres):
    r"""
    Function returns potential temperautre in Kelvin

    Parameters:
    -----------

    temp(float): Temperature value in Kelvin
    pres(float): Pressure value in pascals    
    """
    
    theta = temp*(100000/pres)**(287/1004)
    return theta


T,Z,I,J = temp.shape
tempsquish = temp.reshape(T,Z,I*J, order = 'F')


theta_squished = np.ones((tempsquish.shape[0],tempsquish.shape[1],tempsquish.shape[2]))*np.nan

for i in range(theta_squished.shape[0]):
    for j in range(theta_squished.shape[1]):
        for k in range(theta_squished.shape[2]):
            theta_squished[i,j,k] = potential_temp(tempsquish[i,j,k], level[j])
 

#Now reshape back into oringal form
        
theta = theta_squished.reshape(T,Z,I,J)        
#%%

#Define a constant latitude and take cross section along this; varies by longitude (-99.5 -93.5 )
#constant latitude is 30degN

def plot_xsect(latitude_index, latitude_value, time, label_size, title_size, sigma_value = None):
    r"""
    Outputs a longitude-height cross-section of potential vorticity and potential temperature with pressure as a vertical coordinate
    
    Parameters:
    -----------
    latitude_index(int): Latitude index such that cross-section is plotted on a constant latitude
    latitude_value(float): The value corresponding to the aforementioned latitude index
    time(int): Time index for the cross-section corresponding to a fixed time
    label_size(float): Label size for x and y axes and colorbar label
    title_size(float): Size of the title
    sigma_value(float): Determines degree of smoothing. If no value is input, default is 1.25
        
    """
        
    
    


    theta_bill = theta[:,:,latitude_index,:]
    pv_bill = pv[:,:,latitude_index,:]*10**6


   

    if sigma_value is None:
        sigma = 1.25
        
    theta_smooth = gaussian_filter(theta_bill, sigma)


#Define a time index (in 3hr increments starting at 00 UTC 6/19/2015)
    
    date_time = input('Enter the date and time of the image:')

    fig = plt.figure(figsize = [10,10])
   

    plt.ylim(level_mb[36],level_mb[11])
    clevs = np.arange(0,10,0.5)
    plt.contour(longitude,level_mb,pv_bill[time,:,:],clevs,colors = 'black')
    cp = plt.contourf(longitude, level_mb,pv_bill[time,:,:],clevs, cmap = 'pyart_NWSRef')
    clevs2 = np.arange(210,370,10)
    plt.contour(longitude,level_mb,theta_smooth[time,:,:],colors = 'red',linewidths = 2)
    cs = plt.contour(longitude,level_mb,theta_smooth[time,:,:], clevs2,colors = 'red',linewidths = 2)


    plt.clabel(cs,inline = 1, fontsize = 10, fmt='%4.0f')
    cbar = plt.colorbar(cp,ticks = clevs[::2], orientation = 'horizontal')
    cbar.set_label('PVU', size = label_size)

    plt.xlabel('Longitude', size = label_size)
    plt.ylabel('Pressure (mb)', size = label_size)
    plt.title(r'Potential Vorticity and $\Theta$ at Latitude ='+str(latitude_value)+'$^{o}$N ' + str(date_time), size = title_size)
    plt.xticks(fontsize = label_size)
    plt.yticks(fontsize = label_size)


    figure = plt.show()
    
    return figure
    
    
plot = plot_xsect(11,38,3,24,26)
