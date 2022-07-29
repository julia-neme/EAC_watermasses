# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 22:48:09 2018

@author: hea144
"""
import matplotlib.pyplot as plt
#import seawater as sw

from netCDF4 import *

import numpy as np
import os
import seawater.eos80 as sw # pip install seawater
import glob
#ds_underw = xar.open_dataset(r'..\..\underway\in2021_v03uwy.nc')

#file = r"c:\python_test\CTD_23_comparison\VL_70931_190923134329.vp2"
#file = r"..\rapidcast\VL_70931_210508024325.vp2"
#list_FILES = sorted(glob.glob(r'..\rapidcast\VL_70931_210520*vp2'))
list_FILES = sorted(glob.glob(r'..\DATA\RAPIDCAST\VL_*vp2'))
N_FILES = len(list_FILES)

plot_folder = '../FIGURES/'
 
def fun_read(file):
    with open(file,'r') as rapidctdfile:
        line = rapidctdfile.readline()
        print(line)
        while "[DATA]" not in line:
            line = rapidctdfile.readline()
            if "Latitude" in line:                
                    print(line)
                    latitude = float(line.split("=",1)[1])
            if "Longitude" in line:                
                    print(line)
                    longitude = float(line.split("=",1)[1])
        #read the data headers
        header = rapidctdfile.readline().split()
        #read the data units
        units = rapidctdfile.readline().split()
        print(header)
        print(units)
        dataline = rapidctdfile.readline()
        pressure = []
        temp = []
        conductivity = []
        rapid_file_salinity = []
        
       
        #expected data format
        #2019/09/15 03:46:16.358	2.869	2.885	21.322	49.694	35.340	1525.495	1024.677	731486
        
        rapidctd_dataset = []
        while dataline:
            datastrings = dataline.split()
            
            #print(datastrings)
            numbers = []
            for i in range(2,len(datastrings)):
                #print(dataline[i])
                numbers.append(float(datastrings[i]))
            pressure.append(numbers[1])
            temp.append(numbers[2])
            conductivity.append(numbers[3])
            rapid_file_salinity.append(numbers[4])
           # dataset.append(numbers)
            dataline = rapidctdfile.readline()
            
       # derive salinity     
        ref_cond = 42.914
        rt = [x/ref_cond for x in conductivity]    
        salinity = sw.salt(rt,temp,pressure)    
        # plt.plot(salinity,pressure)

        return pressure, temp, salinity, datastrings, header, latitude, longitude
    
# Loop
list_CAST = []
for f in range(N_FILES):
    CAST = {}
    file = list_FILES[f]
    print(file)        
    [pressure, temp, salinity, datastrings, header, latitude, longitude] = fun_read(file)
#    plt.plot(salinity,pressure)
    CAST['TEMP'] = temp
    CAST['SAL'] = salinity
    CAST['PRES'] = pressure
    list_CAST.append(CAST)
    
    
plt.figure(figsize=(10, 12))
for f in range(N_FILES):
    plt.plot(list_CAST[f]['TEMP'],list_CAST[f]['PRES'],label = datastrings[0])
plt.gca().invert_yaxis()
plt.title('Rapid cast Temperature, date ' + file[22:28],fontsize=22)          
plt.xlabel('Temperature [$^oC$]',fontsize=16)                  
plt.ylabel('Pressure dbar]',fontsize=16)    
plt.legend()              
# SAVE
plt.tight_layout()
plt.savefig(plot_folder + 'plot_Rapid_Cast_temp' + file[22:28] + '.png', bbox_inches='tight', pad_inches=0.5, dpi=300)
plt.show()


plt.figure(figsize=(10, 12))
for f in range(N_FILES):
    plt.plot(list_CAST[f]['SAL'],list_CAST[f]['PRES'])
plt.gca().invert_yaxis()
plt.title('Rapid cast Salinity, date ' + file[22:28],fontsize=16)          
plt.xlabel('Salinity',fontsize=13)                  
plt.ylabel('Pressure dbar]',fontsize=13)                  
# SAVE
plt.tight_layout()
plt.savefig(plot_folder + 'plot_Rapid_Cast_sal' + file[22:28] + '.png', bbox_inches='tight', pad_inches=0.5, dpi=300)
plt.show()


