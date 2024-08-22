#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 18:41:11 2024

@author: usuario
"""
import numpy as np
from  itertools import groupby 

import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams['mathtext.default'] = 'regular' 

from methods import AnalysisMethods, MitoAnalysis

mAnalysis = AnalysisMethods()



def plot_Kdist(expData, ax ,  outlier_threshold=4):
    ## plot parameters
    param_no_outliers = {'s':50, 'facecolors':'lightgray', 'edgecolors': 'dimgray','linewidth':0.8, 'alpha':0.8}
    param_outliers = {'s': 70,  'linewidth':0.8, 'alpha':0.7}
    pal_outliers = sns.color_palette(['#a52c60','#fb9b06', '#6d6dda', '#a0da39' ,'#1fa187'])
    
    ## input: mito information and data
    mito = MitoAnalysis(expData)   
    Time = mito.Time
    K = mito.Ks_mito()
    num_events = mito.n_events(outlier_threshold)
    
    ax.scatter(Time[:-1], K, **param_no_outliers)
    ax.set_ylabel(r' K [$\mu m^2$]', fontsize=20)
    
    ## detect outliers and group in events
    # finding the iqr region
    q1, q3 = np.quantile(K, 0.25) , np.quantile(K, 0.75)     
    iqr = q3-q1     
    # finding upper and lower whiskers
    upper_bound = q3+(outlier_threshold*iqr) #outlier threshold
    ax.axhline(y=upper_bound, ls='--', color='royalblue', zorder=1) 
    
    ind_outliers, outliers = mAnalysis.outliers_Kj( mito.Xint, mito.Yint, outlier_threshold)

    #group outliers in events for plotting
    index_True = []
    size_True = []
    index = 0
    for k, g in groupby(ind_outliers): #k is True/False, g is the group of True/False
        a = len(list(g)) #length of each group
        
        if k: 
            index_True.append(index) #first apparition
            size_True.append(a)
        index = index + a 
    
    col_outlier = sns.color_palette(pal_outliers, n_colors=len(index_True))
    #plot events with different colors
    for N in range(len(index_True)):
        ax.scatter(Time[:-1][index_True[N]:index_True[N]+size_True[N]], 
                   K[index_True[N]:index_True[N]+size_True[N]], 
                   color=col_outlier[N], **param_outliers)    
    
    return num_events


def plot_CSD(expData, ax , outlier_threshold=1.5):
    ## plot parameters
    param_no_outliers = {'s':50, 'facecolors':'lightgray', 'edgecolors': 'dimgray','linewidth':0.8, 'alpha':0.8}
    param_outliers = {'s': 70,  'linewidth':0.8, 'alpha':0.7}
    pal_outliers = sns.color_palette(['#a52c60','#fb9b06', '#6d6dda', '#a0da39' ,'#1fa187'])
    
    # experimental data and analysis
    mito = MitoAnalysis(expData)
    Time = mito.Time
    CSD = mito.CSD()
    num_events = mito.n_events(outlier_threshold)
    
    ax.scatter(Time, CSD, **param_no_outliers)
    ax.set_ylabel(r' CSD [$\mu m^2$]', fontsize=20)
    
    
    ind_outliers, outliers = mAnalysis.outliers_Kj( mito.Xint, mito.Yint, outlier_threshold)
    
   #group outliers in events for plotting
    index_True = []
    size_True = []
    index = 0
    for k, g in groupby(ind_outliers): 
        a = len(list(g)) 
        if k: 
            index_True.append(index) 
            size_True.append(a)
        index = index + a 
    
    col_outlier = sns.color_palette(pal_outliers, n_colors=len(index_True))
    
    for N in range(len(index_True)):
        ax.scatter(Time[1:][index_True[N]:index_True[N]+size_True[N]], 
                   CSD[1:][index_True[N]:index_True[N]+size_True[N]], 
                   color=col_outlier[N], **param_outliers)
    
    return num_events



file ='snake_fig3b.txt'


fig, ax = plt.subplots(2,1, figsize=(8,6),sharex=True)

Nevents = plot_Kdist(file,  outlier_threshold=4 , ax=ax[0])
ax[0].tick_params(axis='both', labelsize=16)

plot_CSD(file, outlier_threshold=4, ax=ax[1])
ax[1].set_xticks([0,50,100,150,200,250],[0,50,100,150,200,250], fontsize=16);
ax[1].set_xlabel('time [s]', fontsize=20);
ax[1].tick_params(axis='both', labelsize=16)
ax[1].set_ylim(-0.10,4)

print('number of events: ', Nevents)