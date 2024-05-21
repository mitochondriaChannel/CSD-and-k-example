#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 20 16:32:50 2024

@author: azulbrigante
"""

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from methods import  methods_CSD


def frame_coordinates(data, frame): 
    # get indexes where the row change frame
    shift_idx = np.unique(data[:,0],return_index=True)[1][1:]
    # split array to thoose indexes
    out = np.split(data,shift_idx) #return list of len=nro frames 
    coordenadas = out[frame -1][:,2:]
    return(coordenadas)

def open_file(expData):
    arch = np.loadtxt(expData) 
    last_frame = int(np.unique(arch[:,0])[-1])
    coord = frame_coordinates(arch, 0)
    Xit, Yit = {}, {}
    Time = np.unique(arch[:,1])
    
    for  frame in range(1, last_frame+1):
        coord = frame_coordinates(arch, frame)
        Xit[frame-1], Yit[frame-1]  = coord[:,0], coord[:,1]
    
    return Xit, Yit, Time/1000 #return time in seconds 

def index_2time(time, t):
    return np.argmin(abs(time-t))

#importmethods to calculate CSD
mCSD =  methods_CSD()

#open file
path =  'formas/XTP/05-06-17XTP/'
file = path + 'snake_img_35_m3.txt'
Xit, Yit, Time = open_file(file)
XitI, YitI = mCSD.expDataInterpolator(Xit, Yit, Time) #interpolate experimental data so the mitochondria have the same amount of beads per frame

#calculate CSD and k
CSD =  mCSD.CSD(XitI, YitI)
k = np.diff(CSD)
#detect events
outlyers_index, detected_events_indexes =  mCSD.CSD_events(CSD)


q1 = np.quantile(k, 0.25)    
# finding the 3rd quartile
q3 = np.quantile(k, 0.75)
# finding the iqr region
iqr = q3-q1

# finding upper  whisker
upper_bound = q3+(4*iqr)

#color palette for events
pallette_outliers = sns.color_palette('Spectral')

#plot K and CSD alongside the indentefied events
fig, axs =plt.subplots(2, 1,  figsize = (6, 5))
ax = axs[0]
#ax.scatter(Time[:-1], k,  color = '#88A635', edgecolor =  '#88A635', alpha = 0.7)
ax.scatter(Time[:-1], k,  color = 'w', edgecolor =  'k')
ax.hlines(upper_bound , 0, Time[-1], color = 'r', linestyle = '--')

for c, j_in in enumerate( detected_events_indexes):
    if len(j_in)==1:
        ax.scatter(Time[:-1][j_in[0]], k[j_in[0] ], s=60, edgecolor =  'w',  color = 'w')
        ax.scatter(Time[:-1][j_in[0]], k[j_in[0] ], alpha = 0.7, s= 60, edgecolor =  pallette_outliers[c],  color = pallette_outliers[c])
    else:
        ax.scatter(Time[:-1][j_in[0]:j_in[1]+1], k[j_in[0]:j_in[1]+1],  s = 60, edgecolor =  'w',  color = 'w')
        ax.scatter(Time[:-1][j_in[0]:j_in[1]+1], k[j_in[0]:j_in[1]+1],  s = 60, alpha = 0.7, edgecolor =  pallette_outliers[c],   color = pallette_outliers[c])


ax.set_ylabel(r'K [$\mu$m$^2$]', fontsize = 20)
ax.tick_params(axis = 'both', labelsize = 15)
ax.set_xlim(0, Time[-1])


ax = axs[1]
ax.scatter(Time, CSD,  color = 'w', edgecolor =  'k')
for c, j_in in enumerate( detected_events_indexes):
    if len(j_in)==1:
        ax.scatter(Time[:-1][j_in[0]], CSD[j_in[0] ], s=60, edgecolor =  'w',   color = 'w')
        ax.scatter(Time[:-1][j_in[0]], CSD[j_in[0] ], s = 60, alpha = 0.7, edgecolor =  pallette_outliers[c],  color = pallette_outliers[c])
    else:
        ax.scatter(Time[:-1][j_in[0]:j_in[1]+1], CSD[j_in[0]:j_in[1]+1],   s = 60,  color = 'w', edgecolor = 'w')
        ax.scatter(Time[:-1][j_in[0]:j_in[1]+1], CSD[j_in[0]:j_in[1]+1], s = 60, alpha = 0.7, edgecolor =  pallette_outliers[c],   color = pallette_outliers[c])

ax.tick_params(axis = 'both', labelsize = 15)
ax.set_xlabel('time [s]', fontsize = 20)
ax.set_ylabel(r'CSD [$\mu$m$^2$]', fontsize = 20)
ax.set_xlim(0, Time[-1])

plt.tight_layout()












