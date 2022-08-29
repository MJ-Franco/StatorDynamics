#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 14:00:05 2021

@author: mariajose

Representation of bar plots
"""


import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
#from scipy.optimize import curve_fit


l=3000 #This numbers limits the number of steps we represent
nmax=13

#%%

files2=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/GTraces_Res_J=5.00  _mu=-5.01  .dat'] #Average of resurrection simulations
data2=[]
#files5=['/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/GTraces_Stall_Averages.dat']
#data5=[]

files1=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/GTraces_Res_Variances.dat'] #Average of resurrection simulations
data1=[]
#files3=['/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/GTraces_Stall_Variances.dat']
#data3=[]

for data_file in files2:
    data2.append(np.loadtxt(data_file))
        
#for data_file in files5:
#    data5.append(np.loadtxt(data_file))
    
for data_file in files1:
    data1.append(np.loadtxt(data_file))
        
#for data_file in files3:
#    data3.append(np.loadtxt(data_file))
    
    
Res_avg = data2[0]
#Stall_avg = data5[0]

Res_var = data1[0]
#Stall_var = data3[0]


#%%

Res_avg_r = np.around(Res_avg[:,1]*nmax,1)
#Stall_avg_r = np.around(Stall_avg[:,1]*nmax,1)

sum(Res_avg[:,1])/len(Res_avg)

#Res_std = sum(np.sqrt(Res_var[:,1]*nmax))/len(Res_var)
#Stall_std = sum(np.sqrt(Stall_var[:,1]*nmax))/len(Stall_var)

u_Res, co_Res = np.unique(Res_avg_r,return_counts=True)
#u_Stall, co_Stall = np.unique(Stall_avg_r,return_counts=True)

plt.rcParams["figure.figsize"] = [8.0,6.0]
plt.title('$<\phi>= 0.66$',size=16)
plt.xlabel('Avg steady states (N=13)',size=15)
plt.bar(u_Res,co_Res/len(Res_avg_r),width=1,label='J=5.00, $\mu$=-4.79')
plt.vlines(8.58, 0, 0.25, linestyles='dashed')

plt.legend(prop={'size':12})

#plt.title(r'Stall $<N>=8 \pm$ {:1.2f}'.format(Stall_std))
#plt.bar(u_Stall,co_Stall/len(Stall_avg_r),width=0.5)
#plt.vlines(8, 0, 0.18, linestyles='dashed')



