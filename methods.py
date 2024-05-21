#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 20 17:41:41 2024

@author: azulbrigante
"""


import numpy as np
from metodCurv import metodos_curvaturas
from scipy.optimize import curve_fit
import math


def index_2time(time, t):
    return np.argmin(abs(time-t))

def length(x, y):
    l = np.diff(x)**2 + np.diff(y)**2
    return sum(np.sqrt(l))

def prop(x, a):
    return a*x
    
def linear(x, a, b):
    return a*x+b

def jump_extr_detector(outliers_index, threshold = 2):
    
    if len(outliers_index)>1:    
        
        index_j = [[outliers_index[0]]]
        i = np.where(threshold<np.diff(outliers_index))[0]
        
        if len(i)>0:
            index_j.append(outliers_index[1:][i]) 
            index_j.append(outliers_index[:-1][i]) 
            
            index_j.append([outliers_index[-1]])
            
            j_index = np.sort(list(set(np.concatenate(index_j))))
            
        else:
            j_index  = [outliers_index[0], outliers_index[-1]]
        
    else:
        j_index =  outliers_index
        
    return j_index 

def jump_detector(outliers_index):
    j_index =  jump_extr_detector(outliers_index)
    jumps = []
    for i in range(len(j_index)):
        if  j_index[i]+1 in outliers_index:
            jumps.append([j_index[i],  j_index[i+1]]) #see if the event is larger than a point
            
        elif not  j_index[i]+1 in outliers_index and  not j_index[i]-1 in outliers_index:
            jumps.append([j_index[i]])
            
    return jumps, j_index
          
def long_arco( DS):
    suma_parcial = 0
    S = np.zeros(len(DS)) #+1)
    #S[0] = 0
    for k in range(1,len(DS)):
        suma_parcial = suma_parcial + DS[k]
        S[k] = suma_parcial
    return(S)
    
#change from xy coordinates  to tangent angle (theta) as function of S 
def change_coord(coord_x, coord_y):
    n = len(coord_x) -1
    ds = np.zeros(n) #lengths of segments vector
    dt = np.zeros(n) #angles vector
    for k in range(n):
        ds[k] = math.sqrt( (coord_x[k+1]-coord_x[k])**2 + (coord_y[k+1]-coord_y[k])**2) 
        dt[k] = math.atan2( (coord_y[k+1]-coord_y[k]),(coord_x[k+1]-coord_x[k]))
    #check if the angle is continuous
    for j in range(len(dt)-1):   
        if  abs(dt[j+1]-dt[j])> math.pi: #dt[j+1]*dt[j] < 0 and
            dt[j+1] = dt[j+1] + np.sign(dt[j])*2*math.pi 
    
    S = long_arco(ds)
    return S, dt

#change from stheta coordinates to xy 
def change_coord2( s, t, x0, y0):
    X , Y = np.zeros(len(t)+1), np.zeros(len(t)+1)
    X[0] , Y[0] = x0, y0 
    ds_s =  np.mean(np.diff(s))
    s2 = np.concatenate([s, [s[-1] + ds_s]])
    
    ds = np.diff(s2)
    for i in range(1,len(X)):
        X[i] = X[i-1] + ds[i-1]*np.cos(t[i-1])
        Y[i] = Y[i-1] + ds[i-1]*np.sin(t[i-1])

    return (X, Y)

    
class methods_CSD:
    
    def __init__(self):
        return    

  
    def CSD(self, X, Y): #calcula el area desplazada por unidad de tiempo
        A = np.zeros_like(X[0, :])
        for i in range(1, len(X[0, :])): 
            sumi = (X[:, i] -X[:, i-1])**2 + (Y[:, i] -Y[:, i-1])**2
            
            A[i] = A[i-1] + np.sum(sumi) #le saco el multiplicar por el Ds
        
        N = len(X[:, 0])
        return A/N


    def expDataInterpolator(self, XitExp, YitExp, TimeExp):
        len_interp = int(np.mean([len(XitExp[i]) for i in XitExp.keys()]))
        XitExpI, YitExpI = np.zeros((len_interp+1, len(TimeExp))), np.zeros((len_interp+1, len(TimeExp)))
        for i in range(len(TimeExp)):
            #data exp interpolated
            X = XitExp[i]
            Y=  YitExp[i]
            
            s, t = change_coord(X, Y)
            regular_s = np.linspace(np.min(s), np.max(s), num=len_interp)
            regular_t = np.interp(regular_s, s, t)
            XitExpI[:, i], YitExpI[:, i] = change_coord2(regular_s, regular_t, X[0], Y[0])
            
        return XitExpI, YitExpI 

    
 

    def CSD_events(self, CSD):      
        k = CSD[1:] - CSD[:-1]
        
        # finding the 1st quartile
        q1 = np.quantile(k, 0.25)    
        # finding the 3rd quartile
        q3 = np.quantile(k, 0.75)
         
        # finding the iqr region
        iqr = q3-q1
         
        # finding upper and lower whiskers
        upper_bound = q3+(4*iqr)
    
        outlyers_indexes = np.where((k >= upper_bound))[0]
        jumps_indexes, j_indexes = jump_detector(outlyers_indexes)
        
        
        return outlyers_indexes, jumps_indexes
    
    
    
    
    
    
    
    
  
    
    