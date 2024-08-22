#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 11:15:55 2023

@author: agus

Clase para analizar las simulaciones 

"""
import seaborn as sns

import numpy as np
import pandas as pd

from  itertools import groupby  






class AnalysisMethods:
    
    def __init__(self):
        return    

    def Ks(self, X, Y): 
        K = np.zeros(len(X[0,:])-1)
        Nseg = len(X[:,0])
        for i in range(1, len(X[0,:])-1): 
            Ki = (X[:, i] -X[:, i-1])**2 + (Y[:, i] -Y[:, i-1])**2 
            K[i] = np.sum(Ki)/Nseg
        return K


    def CSD(self, X, Y): 
        A = np.zeros(len(X[0,:]))
        Ki = self.Ks( X, Y)
        A[1:] = np.cumsum(Ki)
        return A
    
    def outliers_Kj(self, X, Y, crit_outlier=1.5):
        Kj = self.Ks( X, Y)       
        # finding the 1st and 3rd quartile
        q1, q3 = np.quantile(Kj, 0.25) , np.quantile(Kj, 0.75)     
        # finding the iqr region
        iqr = q3-q1     
        # finding upper and lower whiskers
        upper_bound = q3+(crit_outlier*iqr)
        lower_bound = q1-(crit_outlier*iqr)
        
        ind_outliers = (Kj <= lower_bound) | (Kj >= upper_bound) 
        
        # Find the indices where the array is True
        true_indices = np.argwhere(ind_outliers)
        # Iterate over the true indices and merge groups separated by one false value
        for i in range(1, len(true_indices)):
            if np.all(true_indices[i] - true_indices[i - 1] == 2):
                ind_outliers[true_indices[i]-1] = True
        
        outliers = Kj[ind_outliers]  
     
        return(ind_outliers, outliers)
    
    def Kj_filtrado(self, X, Y, crit_outlier=1.5):
        Kj = self.Ks( X, Y)
           
        ind_outliers, outliers = self.outliers_Kj(X, Y, crit_outlier)
     
        if len(outliers)>1:
            Kj_filt = Kj[ ~ ind_outliers]
        else:
            Kj_filt = Kj
        return(Kj_filt)

    
    def Nevents(self, X, Y, crit_outlier=1.5):
        ind_outliers, outliers = self.outliers_Kj(X, Y,crit_outlier)
        try:
            num_events = pd.Series([k for k, g in groupby(ind_outliers)]).value_counts()[True]
        except KeyError:
            num_events = 0
            
        return(num_events)


    def length(self, x, y):
        l = np.diff(x)**2 + np.diff(y)**2
        return sum(np.sqrt(l))
    
   




class MitoAnalysis(AnalysisMethods):

    def __init__(self, track):
        
        self.Xexp, self.Yexp, self.Time =  self.open_file_exp(track)

        self.DT_samp = self.Time[1] - self.Time[0]
        self.L = np.array([self.length(self.Xexp[ i], self.Yexp[i]) for i in range(len(self.Time))])
        
        self.Xint, self.Yint = self.expDataInterpolator() 
        self.Xcm = np.array([np.mean(self.Xexp[i]) for i  in range(len(self.Time))])
        self.Ycm = np.array([np.mean(self.Yexp[i]) for i  in range(len(self.Time))])


        return   
    
    def open_file_exp(self, expData):
        arch = np.loadtxt(expData) 
        last_frame = int(np.unique(arch[:,0])[-1])
        coord = self.coord_frame(arch, 0)
        Xit, Yit = {}, {}
        Time = np.unique(arch[:,1])
        
        if Time[1]>800: 
            for  frame in range(1, last_frame+1):
                coord = self.coord_frame(arch, frame)
                Xit[frame-1], Yit[frame-1]  = coord[:,0], coord[:,1]
            Time_F = Time
        else: #resample if sampling time < 800ms
            for  frame in range(1, last_frame+1,2):
                coord = self.coord_frame(arch, frame)
                Xit[int((frame-1)/2)], Yit[int((frame-1)/2)]  = coord[:,0], coord[:,1]
            Time_F = Time[::2]
                       
        return Xit, Yit, Time_F/1000

    def coord_frame(self, data, frame): 
        shift_idx = np.unique(data[:,0],return_index=True)[1][1:]
        
        out = np.split(data,shift_idx) 
        
        coord = out[frame -1][:,2:]
        return(coord)
    
    def arc_length(self, DS):
        partial_sum = 0
        S = np.zeros(len(DS)) #+1)
        #S[0] = 0
        for k in range(1,len(DS)):
            partial_sum = partial_sum + DS[k]
            S[k] = partial_sum
        return(S)

    def coord_change(self, coord_x, coord_y):
        n = len(coord_x) -1
        ds = np.zeros(n) 
        dt = np.zeros(n) 
        for k in range(n):
            ds[k] = np.sqrt( (coord_x[k+1]-coord_x[k])**2 + (coord_y[k+1]-coord_y[k])**2) 
            dt[k] = np.arctan2( (coord_y[k+1]-coord_y[k]),(coord_x[k+1]-coord_x[k]))
        #check if angle is continuous
        for j in range(len(dt)-1):   
            if  abs(dt[j+1]-dt[j])> np.pi: 
                dt[j+1] = dt[j+1] + np.sign(dt[j])*2*np.pi     
        S = self.arc_length(ds)
        return S, dt
    
    def coord_change2(self, s, t, x0, y0):
        X , Y = np.zeros(len(t)+1), np.zeros(len(t)+1)
        X[0] , Y[0] = x0, y0 
        ds_s =  np.mean(np.diff(s))
        s2 = np.concatenate([s, [s[-1] + ds_s]])
    
        ds = np.diff(s2)
        for i in range(1,len(X)):
            X[i] = X[i-1] + ds[i-1]*np.cos(t[i-1])
            Y[i] = Y[i-1] + ds[i-1]*np.sin(t[i-1])
    
        return (X, Y)
    
    def expDataInterpolator(self):
        len_interp = int(np.mean([len(self.Xexp[i]) for i in self.Xexp.keys()]))
        XitExpI, YitExpI = np.zeros((len_interp+1, len(self.Time))), np.zeros((len_interp+1, len(self.Time)))
        for i in range(len(self.Time)):
            s, t = self.coord_change(self.Xexp[i], self.Yexp[i])
            regular_s = np.linspace(np.min(s), np.max(s), num=len_interp)
            regular_t = np.interp(regular_s, s, t)
            
            XitExpI[:, i], YitExpI[:, i] = self.coord_change2(regular_s, regular_t, self.Xexp[i][0], self.Yexp[i][0])
    
        return XitExpI, YitExpI 
    
    def Ks_mito(self): 
        K = super().Ks(self.Xint,self.Yint)    
        return K

    def CSD(self): 
        A = super().CSD(self.Xint, self.Yint)
        return A
    
    def Kj_fil(self,crit_outlier=1.5):
        Kj_filt = super().Kj_filtrado(self.Xint, self.Yint, crit_outlier)
        return(Kj_filt)

    def Kmean(self, crit_outlier=1.5):
        Kj = self.Kj_fil( crit_outlier )
        return(np.mean(Kj))
    
    def n_events(self, crit_outlier=1.5):
        num_events = super().Nevents(self.Xint, self.Yint, crit_outlier)
        return(num_events)


