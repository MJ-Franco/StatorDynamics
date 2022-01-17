#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 12 12:07:02 2021

Analysis of dwell times for J=0 kT

@author: mariajose
"""

import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from scipy.optimize import curve_fit
from numba import jit


#%% Data load

files1=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/GTraces_Res_J=0.00  _mu=-1.61  .dat'] #All resurrection simulations
data1=[]
files2=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/GTraces_Res_J=0.00  _mu= 1.61  .dat'] #All resurrection simulations
data2=[]
files3=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/GTraces_Res_J=0.00  _mu= 0.00  .dat'] #All resurrection simulations
data3=[]

for data_file in files1:
    data1.append(np.loadtxt(data_file))
    
for data_file in files2:
    data2.append(np.loadtxt(data_file))
    
for data_file in files3:
    data3.append(np.loadtxt(data_file))
    
All_Res_mu_m69=data1[0]
All_Res_mu_69=data2[0]
All_Res_mu_0=data3[0]

#%%

sim = 10000
stp = 3000

nmax=13

k = np.arange(0,11,1)


#%% Dwell time for koff = 10kon variable definition

dwell_k0_mu_m69 = np.zeros(sim)
dwell_k1_mu_m69 = np.zeros(sim)
dwell_k2_mu_m69 = np.zeros(sim)
dwell_k3_mu_m69 = np.zeros(sim)
dwell_k4_mu_m69 = np.zeros(sim)
dwell_k5_mu_m69 = np.zeros(sim)
dwell_k6_mu_m69 = np.zeros(sim)
dwell_k7_mu_m69 = np.zeros(sim)
dwell_k8_mu_m69 = np.zeros(sim)
dwell_k9_mu_m69 = np.zeros(sim)
dwell_k10_mu_m69 = np.zeros(sim)


start2_m69 = np.zeros(sim)
start3_m69 = np.zeros(sim)
start4_m69 = np.zeros(sim)
start5_m69 = np.zeros(sim)
start6_m69 = np.zeros(sim)
start7_m69 = np.zeros(sim)
start8_m69 = np.zeros(sim)
start9_m69 = np.zeros(sim)
start10_m69 = np.zeros(sim)

#%% Dwell time for koff = 10kon calculation


for i in range(sim):
    for j in range(stp*i,stp*(i+1)-1):
        if All_Res_mu_m69[j,2]*nmax<1.:
            dwell_k0_mu_m69[i] = dwell_k0_mu_m69[i] + 1.
        else:
            break

for i in range(sim):     
    for j in range(int(stp*i+dwell_k0_mu_m69[i]),stp*(i+1)-1):
        if All_Res_mu_m69[j,2]*nmax==1.:
            dwell_k1_mu_m69[i] = dwell_k1_mu_m69[i] + 1.
        else:
            break
        
for i in range(sim):
    for j in range(int(stp*i+dwell_k1_mu_m69[i]),stp*(i+1)-1):
        if All_Res_mu_m69[j,2]*nmax==2.:
            start2_m69[i] = j           
            break        
for i in range(sim):
    for j in range(int(start2_m69[i]),stp*(i+1)-1):
            if All_Res_mu_m69[j,2]*nmax==2.:
                dwell_k2_mu_m69[i] = dwell_k2_mu_m69[i] + 1.
            else:
                break
            
for i in range(sim):
    for j in range(int(stp*i+dwell_k2_mu_m69[i]),stp*(i+1)-1):
        if All_Res_mu_m69[j,2]*nmax==3.:
            start3_m69[i] = j
            break
for i in range(sim):
    for j in range(int(start3_m69[i]),stp*(i+1)-1):
        if All_Res_mu_m69[j,2]*nmax==3.:
            dwell_k3_mu_m69[i] = dwell_k3_mu_m69[i] + 1.
        else:
            break

for i in range(sim):
    for j in range(int(stp*i+dwell_k3_mu_m69[i]),stp*(i+1)-1):
        if All_Res_mu_m69[j,2]*nmax==4.:
            start4_m69[i] = j
            break
for i in range(sim):
    for j in range(int(start4_m69[i]),stp*(i+1)-1):
        if All_Res_mu_m69[j,2]*nmax==4.:
            dwell_k4_mu_m69[i] = dwell_k4_mu_m69[i] + 1.
        else:
            break
        
for i in range(sim):
    for j in range(int(stp*i+dwell_k4_mu_m69[i]),stp*(i+1)-1):
        if All_Res_mu_m69[j,2]*nmax==5.:
            start5_m69[i] = j
            break
for i in range(sim):
    for j in range(int(start5_m69[i]),stp*(i+1)-1):
        if All_Res_mu_m69[j,2]*nmax==5.:
            dwell_k5_mu_m69[i] = dwell_k5_mu_m69[i] + 1.
        else:
            break
        
for i in range(sim):
    for j in range(int(stp*i+dwell_k5_mu_m69[i]),stp*(i+1)-1):
        if All_Res_mu_m69[j,2]*nmax==6.:
            start6_m69[i] = j
            break
for i in range(sim):
    for j in range(int(start6_m69[i]),stp*(i+1)-1):
        if All_Res_mu_m69[j,2]*nmax==6.:
            dwell_k6_mu_m69[i] = dwell_k6_mu_m69[i] + 1.
        else:
            break


for i in range(sim):
    for j in range(int(stp*i+dwell_k6_mu_m69[i]),stp*(i+1)-1):
        if All_Res_mu_m69[j,2]*nmax==7.:
            start7_m69[i] = j
            break
for i in range(sim):
    for j in range(int(start7_m69[i]),stp*(i+1)-1):
        if All_Res_mu_m69[j,2]*nmax==7.:
            dwell_k7_mu_m69[i] = dwell_k7_mu_m69[i] + 1.
        else:
            break    
              
for i in range(sim):
    for j in range(int(stp*i+dwell_k7_mu_m69[i]),stp*(i+1)-1):
        if All_Res_mu_m69[j,2]*nmax==8.:
            start8_m69[i] = j
            break
for i in range(sim):
    for j in range(int(start8_m69[i]),stp*(i+1)-1):
        if All_Res_mu_m69[j,2]*nmax==8.:
            dwell_k8_mu_m69[i] = dwell_k8_mu_m69[i] + 1.
        else:
            break 

for i in range(sim):
    for j in range(int(stp*i+dwell_k8_mu_m69[i]),stp*(i+1)-1):
        if All_Res_mu_m69[j,2]*nmax==9.:
            start9_m69[i] = j
            break
for i in range(sim):
    for j in range(int(start9_m69[i]),stp*(i+1)-1):
        if All_Res_mu_m69[j,2]*nmax==9.:
            dwell_k9_mu_m69[i] = dwell_k9_mu_m69[i] + 1.
        else:
            break 
   
    
for i in range(sim):
    for j in range(int(stp*i+dwell_k9_mu_m69[i]),stp*(i+1)-1):
        if All_Res_mu_m69[j,2]*nmax==10.:
            start10_m69[i] = j
            break
for i in range(sim):
    for j in range(int(start10_m69[i]),stp*(i+1)-1):
        if All_Res_mu_m69[j,2]*nmax==10.:
            dwell_k10_mu_m69[i] = dwell_k10_mu_m69[i] + 1.
        else:
            break    


  
#%% Dwell time for kon = 10koff variable definition

dwell_k0_mu_69 = np.zeros(sim)
dwell_k1_mu_69 = np.zeros(sim)
dwell_k2_mu_69 = np.zeros(sim)
dwell_k3_mu_69 = np.zeros(sim)
dwell_k4_mu_69 = np.zeros(sim)
dwell_k5_mu_69 = np.zeros(sim)
dwell_k6_mu_69 = np.zeros(sim)
dwell_k7_mu_69 = np.zeros(sim)
dwell_k8_mu_69 = np.zeros(sim)
dwell_k9_mu_69 = np.zeros(sim)
dwell_k10_mu_69 = np.zeros(sim)

start2_69 = np.zeros(sim)
start3_69 = np.zeros(sim)
start4_69 = np.zeros(sim)
start5_69 = np.zeros(sim)
start6_69 = np.zeros(sim)
start7_69 = np.zeros(sim)
start8_69 = np.zeros(sim)
start9_69 = np.zeros(sim)
start10_69 = np.zeros(sim)



#%% Dwell time for kon = 10koff calculation


for i in range(sim):
    for j in range(stp*i,stp*(i+1)-1):
        if All_Res_mu_69[j,2]*nmax<1.:
            dwell_k0_mu_69[i] = dwell_k0_mu_69[i] + 1.
        else:
            break

for i in range(sim):     
    for j in range(int(stp*i+dwell_k0_mu_69[i]),stp*(i+1)-1):
        if All_Res_mu_69[j,2]*nmax==1.:
            dwell_k1_mu_69[i] = dwell_k1_mu_69[i] + 1.
        else:
            break
        
for i in range(sim):
    for j in range(int(stp*i+dwell_k1_mu_69[i]),stp*(i+1)-1):
        if All_Res_mu_69[j,2]*nmax==2.:
            start2_69[i] = j           
            break        
for i in range(sim):
    for j in range(int(start2_69[i]),stp*(i+1)-1):
            if All_Res_mu_69[j,2]*nmax==2.:
                dwell_k2_mu_69[i] = dwell_k2_mu_69[i] + 1.
            else:
                break
            
for i in range(sim):
    for j in range(int(stp*i+dwell_k2_mu_69[i]),stp*(i+1)-1):
        if All_Res_mu_69[j,2]*nmax==3.:
            start3_69[i] = j
            break
for i in range(sim):
    for j in range(int(start3_69[i]),stp*(i+1)-1):
        if All_Res_mu_69[j,2]*nmax==3.:
            dwell_k3_mu_69[i] = dwell_k3_mu_69[i] + 1.
        else:
            break

for i in range(sim):
    for j in range(int(stp*i+dwell_k3_mu_69[i]),stp*(i+1)-1):
        if All_Res_mu_69[j,2]*nmax==4.:
            start4_69[i] = j
            break
for i in range(sim):
    for j in range(int(start4_69[i]),stp*(i+1)-1):
        if All_Res_mu_69[j,2]*nmax==4.:
            dwell_k4_mu_69[i] = dwell_k4_mu_69[i] + 1.
        else:
            break
        
for i in range(sim):
    for j in range(int(stp*i+dwell_k4_mu_69[i]),stp*(i+1)-1):
        if All_Res_mu_69[j,2]*nmax==5.:
            start5_69[i] = j
            break
for i in range(sim):
    for j in range(int(start5_69[i]),stp*(i+1)-1):
        if All_Res_mu_69[j,2]*nmax==5.:
            dwell_k5_mu_69[i] = dwell_k5_mu_69[i] + 1.
        else:
            break
        
for i in range(sim):
    for j in range(int(stp*i+dwell_k5_mu_69[i]),stp*(i+1)-1):
        if All_Res_mu_69[j,2]*nmax==6.:
            start6_69[i] = j
            break
for i in range(sim):
    for j in range(int(start6_69[i]),stp*(i+1)-1):
        if All_Res_mu_69[j,2]*nmax==6.:
            dwell_k6_mu_69[i] = dwell_k6_mu_69[i] + 1.
        else:
            break
        
for i in range(sim):
    for j in range(int(stp*i+dwell_k6_mu_69[i]),stp*(i+1)-1):
        if All_Res_mu_69[j,2]*nmax==7.:
            start7_69[i] = j
            break
for i in range(sim):
    for j in range(int(start7_69[i]),stp*(i+1)-1):
        if All_Res_mu_69[j,2]*nmax==7.:
            dwell_k7_mu_69[i] = dwell_k7_mu_69[i] + 1.
        else:
            break
        
for i in range(sim):
    for j in range(int(stp*i+dwell_k7_mu_69[i]),stp*(i+1)-1):
        if All_Res_mu_69[j,2]*nmax==8.:
            start8_69[i] = j
            break
for i in range(sim):
    for j in range(int(start8_69[i]),stp*(i+1)-1):
        if All_Res_mu_69[j,2]*nmax==8.:
            dwell_k8_mu_69[i] = dwell_k8_mu_69[i] + 1.
        else:
            break        
        
for i in range(sim):
    for j in range(int(stp*i+dwell_k8_mu_69[i]),stp*(i+1)-1):
        if All_Res_mu_69[j,2]*nmax==9.:
            start9_69[i] = j
            break
for i in range(sim):
    for j in range(int(start9_69[i]),stp*(i+1)-1):
        if All_Res_mu_69[j,2]*nmax==9.:
            dwell_k9_mu_69[i] = dwell_k9_mu_69[i] + 1.
        else:
            break        

for i in range(sim):
    for j in range(int(stp*i+dwell_k9_mu_69[i]),stp*(i+1)-1):
        if All_Res_mu_69[j,2]*nmax==10.:
            start10_69[i] = j
            break
for i in range(sim):
    for j in range(int(start10_69[i]),stp*(i+1)-1):
        if All_Res_mu_69[j,2]*nmax==10.:
            dwell_k10_mu_69[i] = dwell_k10_mu_69[i] + 1.
        else:
            break

        
#%% Dwell time for kon = koff variable declaration

dwell_k0_mu_0 = np.zeros(sim)
dwell_k1_mu_0 = np.zeros(sim)
dwell_k2_mu_0 = np.zeros(sim)
dwell_k3_mu_0 = np.zeros(sim)
dwell_k4_mu_0 = np.zeros(sim)
dwell_k5_mu_0 = np.zeros(sim)
dwell_k6_mu_0 = np.zeros(sim)
dwell_k7_mu_0 = np.zeros(sim)
dwell_k8_mu_0 = np.zeros(sim)
dwell_k9_mu_0 = np.zeros(sim)
dwell_k10_mu_0 = np.zeros(sim)

start2_0 = np.zeros(sim)
start3_0 = np.zeros(sim)
start4_0 = np.zeros(sim)
start5_0 = np.zeros(sim)
start6_0 = np.zeros(sim)
start7_0 = np.zeros(sim)
start8_0 = np.zeros(sim)
start9_0 = np.zeros(sim)
start10_0 = np.zeros(sim)

stop0_0 = np.zeros(sim)
stop1_0 = np.zeros(sim)
stop2_0 = np.zeros(sim)
stop3_0 = np.zeros(sim)
stop4_0 = np.zeros(sim)
stop5_0 = np.zeros(sim)
stop6_0 = np.zeros(sim)
stop7_0 = np.zeros(sim)
stop8_0 = np.zeros(sim)
stop9_0 = np.zeros(sim)


#%% Dwell time for kon = koff calculation


for i in range(sim):
    for j in range(stp*i,stp*(i+1)-1):
        if All_Res_mu_0[j,2]*nmax==0.:
            dwell_k0_mu_0[i] = dwell_k0_mu_0[i] + 1.
        else:
            stop0_0[i] = j
            break

for i in range(sim):     
    for j in range(int(stop0_0[i]),stp*(i+1)-1):
        if All_Res_mu_0[j,2]*nmax==1.:
            dwell_k1_mu_0[i] = dwell_k1_mu_0[i] + 1.
        else:
            stop1_0[i] = j
            break
        
for i in range(sim):
    for j in range(int(stop1_0[i]),stp*(i+1)-1):
        if All_Res_mu_0[j,2]*nmax==2.:
            start2_0[i] = j           
            break        
for i in range(sim):
    for j in range(int(start2_0[i]),stp*(i+1)-1):
            if All_Res_mu_0[j,2]*nmax==2.:
                dwell_k2_mu_0[i] = dwell_k2_mu_0[i] + 1.
            else:
                stop2_0[i] = j
                break
            
for i in range(sim):
    for j in range(int(stop2_0[i]),stp*(i+1)-1):
        if All_Res_mu_0[j,2]*nmax==3.:
            start3_0[i] = j
            break
for i in range(sim):
    for j in range(int(start3_0[i]),stp*(i+1)-1):
        if All_Res_mu_0[j,2]*nmax==3.:
            dwell_k3_mu_0[i] = dwell_k3_mu_0[i] + 1.
        else:
            stop3_0[i] = j
            break

for i in range(sim):
    for j in range(int(stop3_0[i]),stp*(i+1)-1):
        if All_Res_mu_0[j,2]*nmax==4.:
            start4_0[i] = j
            break
for i in range(sim):
    for j in range(int(start4_0[i]),stp*(i+1)-1):
        if All_Res_mu_0[j,2]*nmax==4.:
            dwell_k4_mu_0[i] = dwell_k4_mu_0[i] + 1.
        else:
            stop4_0[i] = j
            break
        
for i in range(sim):
    for j in range(int(stop4_0[i]),stp*(i+1)-1):
        if All_Res_mu_0[j,2]*nmax==5.:
            start5_0[i] = j
            break
for i in range(sim):
    for j in range(int(start5_0[i]),stp*(i+1)-1):
        if All_Res_mu_0[j,2]*nmax==5.:
            dwell_k5_mu_0[i] = dwell_k5_mu_0[i] + 1.
        else:
            stop5_0[i] = j
            break
        
for i in range(sim):
    for j in range(int(stop5_0[i]),stp*(i+1)-1):
        if All_Res_mu_0[j,2]*nmax==6.:
            start6_0[i] = j
            break
for i in range(sim):
    for j in range(int(start6_0[i]),stp*(i+1)-1):
        if All_Res_mu_0[j,2]*nmax==6.:
            dwell_k6_mu_0[i] = dwell_k6_mu_0[i] + 1.
        else:
            stop6_0[i] = j
            break

for i in range(sim):
    for j in range(int(stop6_0[i]),stp*(i+1)-1):
        if All_Res_mu_0[j,2]*nmax==7.:
            start7_0[i] = j
            break
for i in range(sim):
    for j in range(int(start7_0[i]),stp*(i+1)-1):
        if All_Res_mu_0[j,2]*nmax==7.:
            dwell_k7_mu_0[i] = dwell_k7_mu_0[i] + 1.
        else:
            stop7_0[i] = j
            break
        
for i in range(sim):
    for j in range(int(stop7_0[i]),stp*(i+1)-1):
        if All_Res_mu_0[j,2]*nmax==8.:
            start8_0[i] = j
            break
for i in range(sim):
    for j in range(int(start8_0[i]),stp*(i+1)-1):
        if All_Res_mu_0[j,2]*nmax==8.:
            dwell_k8_mu_0[i] = dwell_k8_mu_0[i] + 1.
        else:
            stop8_0[i] = j
            break
        
for i in range(sim):
    for j in range(int(stop8_0[i]),stp*(i+1)-1):
        if All_Res_mu_0[j,2]*nmax==9.:
            start9_0[i] = j
            break
for i in range(sim):
    for j in range(int(start9_0[i]),stp*(i+1)-1):
        if All_Res_mu_0[j,2]*nmax==9.:
            dwell_k9_mu_0[i] = dwell_k9_mu_0[i] + 1.
        else:
            stop9_0[i] = j
            break
        
for i in range(sim):
    for j in range(int(stop9_0[i]),stp*(i+1)-1):
        if All_Res_mu_0[j,2]*nmax==10.:
            start10_0[i] = j
            break
for i in range(sim):
    for j in range(int(start10_0[i]),stp*(i+1)-1):
        if All_Res_mu_0[j,2]*nmax==10.:
            dwell_k10_mu_0[i] = dwell_k10_mu_0[i] + 1.
        else:
            break
 

           
#%% Average and sd for koff = 10kon

dwell_avg_k0_mu_m69 = sum(dwell_k0_mu_m69)/sim
dwell_avg2_k0_mu_m69 = sum(dwell_k0_mu_m69**2.)/sim
dwell_sd_k0_mu_m69 = np.sqrt(dwell_avg2_k0_mu_m69-dwell_avg_k0_mu_m69**2)

dwell_avg_k1_mu_m69 = sum(dwell_k1_mu_m69)/sim
dwell_avg2_k1_mu_m69 = sum(dwell_k1_mu_m69**2.)/sim
dwell_sd_k1_mu_m69 = np.sqrt(dwell_avg2_k1_mu_m69-dwell_avg_k1_mu_m69**2)

dwell_avg_k2_mu_m69 = sum(dwell_k2_mu_m69)/sim
dwell_avg2_k2_mu_m69 = sum(dwell_k2_mu_m69**2.)/sim
dwell_sd_k2_mu_m69 = np.sqrt(dwell_avg2_k2_mu_m69-dwell_avg_k2_mu_m69**2)   

dwell_avg_k3_mu_m69 = sum(dwell_k3_mu_m69)/sim
dwell_avg2_k3_mu_m69 = sum(dwell_k3_mu_m69**2.)/sim
dwell_sd_k3_mu_m69 = np.sqrt(dwell_avg2_k3_mu_m69-dwell_avg_k3_mu_m69**2)    

dwell_avg_k4_mu_m69 = sum(dwell_k4_mu_m69)/sim
dwell_avg2_k4_mu_m69 = sum(dwell_k4_mu_m69**2.)/sim
dwell_sd_k4_mu_m69 = np.sqrt(dwell_avg2_k4_mu_m69-dwell_avg_k4_mu_m69**2) 

dwell_avg_k5_mu_m69 = sum(dwell_k5_mu_m69)/sim
dwell_avg2_k5_mu_m69 = sum(dwell_k5_mu_m69**2.)/sim
dwell_sd_k5_mu_m69 = np.sqrt(dwell_avg2_k5_mu_m69-dwell_avg_k5_mu_m69**2)

dwell_avg_k6_mu_m69 = sum(dwell_k6_mu_m69)/sim
dwell_avg2_k6_mu_m69 = sum(dwell_k6_mu_m69**2.)/sim
dwell_sd_k6_mu_m69 = np.sqrt(dwell_avg2_k6_mu_m69-dwell_avg_k6_mu_m69**2)

dwell_avg_k7_mu_m69 = sum(dwell_k7_mu_m69)/sim
dwell_avg2_k7_mu_m69 = sum(dwell_k7_mu_m69**2.)/sim
dwell_sd_k7_mu_m69 = np.sqrt(dwell_avg2_k7_mu_m69-dwell_avg_k7_mu_m69**2)

dwell_avg_k8_mu_m69 = sum(dwell_k8_mu_m69)/sim
dwell_avg2_k8_mu_m69 = sum(dwell_k8_mu_m69**2.)/sim
dwell_sd_k8_mu_m69 = np.sqrt(dwell_avg2_k8_mu_m69-dwell_avg_k8_mu_m69**2)

dwell_avg_k9_mu_m69 = sum(dwell_k9_mu_m69)/sim
dwell_avg2_k9_mu_m69 = sum(dwell_k9_mu_m69**2.)/sim
dwell_sd_k9_mu_m69 = np.sqrt(dwell_avg2_k9_mu_m69-dwell_avg_k9_mu_m69**2)

dwell_avg_k10_mu_m69 = sum(dwell_k10_mu_m69)/sim
dwell_avg2_k10_mu_m69 = sum(dwell_k10_mu_m69**2.)/sim
dwell_sd_k10_mu_m69 = np.sqrt(dwell_avg2_k10_mu_m69-dwell_avg_k10_mu_m69**2)



#%% Average and sd of kon = 10koff

dwell_avg_k0_mu_69 = sum(dwell_k0_mu_69)/sim
dwell_avg2_k0_mu_69 = sum(dwell_k0_mu_69**2.)/sim
dwell_sd_k0_mu_69 = np.sqrt(dwell_avg2_k0_mu_69-dwell_avg_k0_mu_69**2)

dwell_avg_k1_mu_69 = sum(dwell_k1_mu_69)/sim
dwell_avg2_k1_mu_69 = sum(dwell_k1_mu_69**2.)/sim
dwell_sd_k1_mu_69 = np.sqrt(dwell_avg2_k1_mu_69-dwell_avg_k1_mu_69**2)

dwell_avg_k2_mu_69 = sum(dwell_k2_mu_69)/sim
dwell_avg2_k2_mu_69 = sum(dwell_k2_mu_69**2.)/sim
dwell_sd_k2_mu_69 = np.sqrt(dwell_avg2_k2_mu_69-dwell_avg_k2_mu_69**2)   

dwell_avg_k3_mu_69 = sum(dwell_k3_mu_69)/sim
dwell_avg2_k3_mu_69 = sum(dwell_k3_mu_69**2.)/sim
dwell_sd_k3_mu_69 = np.sqrt(dwell_avg2_k3_mu_69-dwell_avg_k3_mu_69**2)    

dwell_avg_k4_mu_69 = sum(dwell_k4_mu_69)/sim
dwell_avg2_k4_mu_69 = sum(dwell_k4_mu_69**2.)/sim
dwell_sd_k4_mu_69 = np.sqrt(dwell_avg2_k4_mu_69-dwell_avg_k4_mu_69**2) 

dwell_avg_k5_mu_69 = sum(dwell_k5_mu_69)/sim
dwell_avg2_k5_mu_69 = sum(dwell_k5_mu_69**2.)/sim
dwell_sd_k5_mu_69 = np.sqrt(dwell_avg2_k5_mu_69-dwell_avg_k5_mu_69**2)

dwell_avg_k6_mu_69 = sum(dwell_k6_mu_69)/sim
dwell_avg2_k6_mu_69 = sum(dwell_k6_mu_69**2.)/sim
dwell_sd_k6_mu_69 = np.sqrt(dwell_avg2_k6_mu_69-dwell_avg_k6_mu_69**2)

dwell_avg_k7_mu_69 = sum(dwell_k7_mu_69)/sim
dwell_avg2_k7_mu_69 = sum(dwell_k7_mu_69**2.)/sim
dwell_sd_k7_mu_69 = np.sqrt(dwell_avg2_k7_mu_69-dwell_avg_k7_mu_69**2)

dwell_avg_k8_mu_69 = sum(dwell_k8_mu_69)/sim
dwell_avg2_k8_mu_69 = sum(dwell_k8_mu_69**2.)/sim
dwell_sd_k8_mu_69 = np.sqrt(dwell_avg2_k8_mu_69-dwell_avg_k8_mu_69**2)    

dwell_avg_k9_mu_69 = sum(dwell_k9_mu_69)/sim
dwell_avg2_k9_mu_69 = sum(dwell_k9_mu_69**2.)/sim
dwell_sd_k9_mu_69 = np.sqrt(dwell_avg2_k9_mu_69-dwell_avg_k9_mu_69**2)

dwell_avg_k10_mu_69 = sum(dwell_k10_mu_69)/sim
dwell_avg2_k10_mu_69 = sum(dwell_k10_mu_69**2.)/sim
dwell_sd_k10_mu_69 = np.sqrt(dwell_avg2_k10_mu_69-dwell_avg_k10_mu_69**2)  


#%% Average and sd of kon = koff

dwell_avg_k0_mu_0 = sum(dwell_k0_mu_0)/sim
dwell_avg2_k0_mu_0 = sum(dwell_k0_mu_0**2.)/sim
dwell_sd_k0_mu_0 = np.sqrt(dwell_avg2_k0_mu_0-dwell_avg_k0_mu_0**2)

dwell_avg_k1_mu_0 = sum(dwell_k1_mu_0)/sim
dwell_avg2_k1_mu_0 = sum(dwell_k1_mu_0**2.)/sim
dwell_sd_k1_mu_0 = np.sqrt(dwell_avg2_k1_mu_0-dwell_avg_k1_mu_0**2)

dwell_avg_k2_mu_0 = sum(dwell_k2_mu_0)/sim
dwell_avg2_k2_mu_0 = sum(dwell_k2_mu_0**2.)/sim
dwell_sd_k2_mu_0 = np.sqrt(dwell_avg2_k2_mu_0-dwell_avg_k2_mu_0**2)   

dwell_avg_k3_mu_0 = sum(dwell_k3_mu_0)/sim
dwell_avg2_k3_mu_0 = sum(dwell_k3_mu_0**2.)/sim
dwell_sd_k3_mu_0 = np.sqrt(dwell_avg2_k3_mu_0-dwell_avg_k3_mu_0**2)    

dwell_avg_k4_mu_0 = sum(dwell_k4_mu_0)/sim
dwell_avg2_k4_mu_0 = sum(dwell_k4_mu_0**2.)/sim
dwell_sd_k4_mu_0 = np.sqrt(dwell_avg2_k4_mu_0-dwell_avg_k4_mu_0**2) 

dwell_avg_k5_mu_0 = sum(dwell_k5_mu_0)/sim
dwell_avg2_k5_mu_0 = sum(dwell_k5_mu_0**2.)/sim
dwell_sd_k5_mu_0 = np.sqrt(dwell_avg2_k5_mu_0-dwell_avg_k5_mu_0**2)

dwell_avg_k6_mu_0 = sum(dwell_k6_mu_0)/sim
dwell_avg2_k6_mu_0 = sum(dwell_k6_mu_0**2.)/sim
dwell_sd_k6_mu_0 = np.sqrt(dwell_avg2_k6_mu_0-dwell_avg_k6_mu_0**2) 

dwell_avg_k7_mu_0 = sum(dwell_k7_mu_0)/sim
dwell_avg2_k7_mu_0 = sum(dwell_k7_mu_0**2.)/sim
dwell_sd_k7_mu_0 = np.sqrt(dwell_avg2_k7_mu_0-dwell_avg_k7_mu_0**2)

dwell_avg_k8_mu_0 = sum(dwell_k8_mu_0)/sim
dwell_avg2_k8_mu_0 = sum(dwell_k8_mu_0**2.)/sim
dwell_sd_k8_mu_0 = np.sqrt(dwell_avg2_k8_mu_0-dwell_avg_k8_mu_0**2)  

dwell_avg_k9_mu_0 = sum(dwell_k9_mu_0)/sim
dwell_avg2_k9_mu_0 = sum(dwell_k9_mu_0**2.)/sim
dwell_sd_k9_mu_0 = np.sqrt(dwell_avg2_k9_mu_0-dwell_avg_k9_mu_0**2)

dwell_avg_k10_mu_0 = sum(dwell_k10_mu_0)/sim
dwell_avg2_k10_mu_0 = sum(dwell_k10_mu_0**2.)/sim
dwell_sd_k10_mu_0 = np.sqrt(dwell_avg2_k10_mu_0-dwell_avg_k10_mu_0**2)  





#%%

k=np.arange(0,11,1)

dwell_avg_mu_m69 = np.array([dwell_avg_k0_mu_m69,dwell_avg_k1_mu_m69,dwell_avg_k2_mu_m69,dwell_avg_k3_mu_m69,
                              dwell_avg_k4_mu_m69,dwell_avg_k5_mu_m69,dwell_avg_k6_mu_m69,dwell_avg_k7_mu_m69,
                              dwell_avg_k8_mu_m69,dwell_avg_k9_mu_m69,dwell_avg_k10_mu_m69])/dwell_avg_k0_mu_m69

dwell_sd_mu_m69 = np.array([dwell_sd_k0_mu_m69,dwell_sd_k1_mu_m69,dwell_sd_k2_mu_m69,dwell_sd_k3_mu_m69,
                             dwell_sd_k4_mu_m69,dwell_sd_k5_mu_m69,dwell_sd_k6_mu_m69,dwell_sd_k7_mu_m69,
                             dwell_sd_k8_mu_m69,dwell_sd_k9_mu_m69,dwell_sd_k10_mu_m69])/dwell_sd_k0_mu_m69


dwell_avg_mu_69 = np.array([dwell_avg_k0_mu_69,dwell_avg_k1_mu_69,dwell_avg_k2_mu_69,dwell_avg_k3_mu_69,
                              dwell_avg_k4_mu_69,dwell_avg_k5_mu_69,dwell_avg_k6_mu_69,dwell_avg_k7_mu_69,
                              dwell_avg_k8_mu_69,dwell_avg_k9_mu_69,dwell_avg_k10_mu_69])/dwell_avg_k0_mu_69

dwell_sd_mu_69 = np.array([dwell_sd_k0_mu_69,dwell_sd_k1_mu_69,dwell_sd_k2_mu_69,dwell_sd_k3_mu_69,##
                              dwell_sd_k4_mu_69,dwell_sd_k5_mu_69,dwell_sd_k6_mu_69,dwell_sd_k7_mu_69,
                              dwell_sd_k8_mu_69,dwell_sd_k9_mu_69,dwell_sd_k10_mu_69])/dwell_sd_k0_mu_69


dwell_avg_mu_0 = np.array([dwell_avg_k0_mu_0,dwell_avg_k1_mu_0,dwell_avg_k2_mu_0,dwell_avg_k3_mu_0,
                              dwell_avg_k4_mu_0,dwell_avg_k5_mu_0,dwell_avg_k6_mu_0,dwell_avg_k7_mu_0,
                              dwell_avg_k8_mu_0,dwell_avg_k9_mu_0,dwell_avg_k10_mu_0])/dwell_avg_k0_mu_0

dwell_sd_mu_0 = np.array([dwell_sd_k0_mu_0,dwell_sd_k1_mu_0,dwell_sd_k2_mu_0,dwell_sd_k3_mu_0,
                             dwell_sd_k4_mu_0,dwell_sd_k5_mu_0,dwell_sd_k6_mu_0,dwell_sd_k7_mu_0,
                             dwell_sd_k8_mu_0,dwell_sd_k9_mu_0,dwell_sd_k10_mu_0])/dwell_sd_k0_mu_0


#%%

plt.rcParams["figure.figsize"] = [8.0,6.0]
plt.grid()
plt.title('$<\phi>=0.17$ ($k_{on}=1/5 k_{off}$)',size=17)
plt.xlabel('k',size=14)
plt.ylabel('<n(k)>/<n(0)>',size=14)
plt.errorbar(k,dwell_avg_mu_m69,dwell_sd_mu_m69,linestyle='',marker='o',capsize=3,label='J=0, $\mu=-1.61$')

plt.legend(prop={'size':15})


#%%   



plt.rcParams["figure.figsize"] = [8.0,6.0]
plt.grid()
plt.title('$<\phi>=0.83$ ($k_{on}=5 k_{off}$)',size=17)
plt.xlabel('k',size=14)
plt.ylabel('<n(k)>/<n(0)>',size=14)
#plt.yticks(np.arange(-0.5, 4.00, 0.50))
plt.errorbar(k,dwell_avg_mu_69,dwell_sd_mu_69,linestyle='',marker='o',capsize=3,label='J=0, $\mu=1.61$')

plt.legend(prop={'size':15})

#%%

plt.rcParams["figure.figsize"] = [8.0,6.0]
plt.grid()
plt.title('$<\phi>=0.5$ ($k_{on}=k_{off}$)',size=17)
plt.xlabel('k',size=14)
plt.ylabel('<n(k)>/<n(0)>',size=14)
plt.errorbar(k,dwell_avg_mu_0,dwell_sd_mu_0,linestyle='',marker='o',capsize=3,label='J=0, $\mu=0$')

plt.legend(prop={'size':15})


#%%


#%%


u_dwell_mu_0, co_dwell_mu_0 = np.unique(dwell_k0_mu_m69,return_counts=True)


def f(t,a,b):
    return a*np.exp(-b*t)

cn = np.linspace(0,np.max(u_dwell_mu_0),1000)

po_dwell_mu_0, pc_dwell_mu_0 = curve_fit(f,u_dwell_mu_0-1,co_dwell_mu_0/sim)
model_dwell_mu_0 = f(cn,*po_dwell_mu_0)

plt.title('P(k=1), $k_{on}=k_{off}$',size=15)
plt.ylabel('Frecuency')
plt.xlabel('n')
plt.bar(u_dwell_mu_0-1,co_dwell_mu_0/sim,color='blue',alpha=0.5)
plt.plot(cn,model_dwell_mu_0,color='red',label='$ae^{-bn}$')
plt.plot([],[],linestyle='',label='a=%.2f'%(po_dwell_mu_0[0]))
plt.plot([],[],linestyle='',label='b=%.2f'%(po_dwell_mu_0[1]))

#plt.bar(u_dwell_k0_mu_m69,co_dwell_k0_mu_m69/sim,color='green',label='koff=2kon',zorder=0)

plt.legend(prop={'size': 15})

a_0 = po_dwell_mu_0[0]
b_0 = po_dwell_mu_0[1]

#%%

a = np.array([a_0,a_1,a_2,a_3,a_4,a_5,a_6,a_7,a_8,a_9,a_10])

b = np.array([b_0,b_1,b_2,b_3,b_4,b_5,b_6,b_7,b_8,b_9,b_10])

plt.title('$J=0$, $\mu=$, ($k_{on}=1/5k_{off}$)',size=17)
plt.xlabel('k',size=14)
plt.plot([],[],linestyle='',label='$ae^{bn}$')
plt.plot(k,a,label='a',linestyle='--',marker='o')
plt.plot(k,b,label='b',linestyle='--',marker='o')

plt.legend(prop={'size': 15})


#%%

'''

u_dwell_mu_69, co_dwell_mu_69 = np.unique(dwell_k10_mu_69,return_counts=True)

cn = np.linspace(0,np.max(u_dwell_mu_69),1000)

po_dwell_mu_69, pc_dwell_mu_69 = curve_fit(f,u_dwell_mu_69-1,co_dwell_mu_69/sim)
model_dwell_mu_69 = f(cn,*po_dwell_mu_69)


plt.title('P(k=0), $J=0$, $\mu=1.61$, ($k_{on}=5k_{off}$)',size=15)
plt.ylabel('Frecuency')
plt.xlabel('n')
plt.bar(u_dwell_mu_69-1,co_dwell_mu_69/sim,color='blue',alpha=0.5)
plt.plot(cn,model_dwell_mu_69,color='red',label='$ae^{-bn}$')
plt.plot([],[],linestyle='',label='a=%.2f'%(po_dwell_mu_69[0]))
plt.plot([],[],linestyle='',label='b=%.2f'%(po_dwell_mu_69[1]))

plt.legend(prop={'size': 15})

a_10 = po_dwell_mu_69[0]
b_10 = po_dwell_mu_69[1]

#%%

a = np.array([a_0,a_1,a_2,a_3,a_4,a_5,a_6,a_7,a_8,a_9,a_10])

b = np.array([b_0,b_1,b_2,b_3,b_4,b_5,b_6,b_7,b_8,b_9,b_10])

plt.title('$J=0$, $\mu=1.69$, ($k_{on}=5k_{off}$)',size=17)
plt.xlabel('k',size=14)
plt.plot([],[],linestyle='',label='$ae^{bn}$')
plt.plot(k,a,label='a',linestyle='--',marker='o')
plt.plot(k,b,label='b',linestyle='--',marker='o')

plt.legend(prop={'size': 15})


#%%


u_dwell_mu_m69, co_dwell_mu_m69 = np.unique(dwell_k10_mu_m69,return_counts=True)

cn = np.linspace(0,np.max(u_dwell_mu_m69),1000)

po_dwell_mu_m69, pc_dwell_mu_m69 = curve_fit(f,u_dwell_mu_m69-1,co_dwell_mu_m69/sim)
model_dwell_mu_m69 = f(cn,*po_dwell_mu_m69)


plt.title('P(k=4), $J=0$, $\mu=-1.61$, ($k_{on}=1/5k_{off}$)',size=15)
plt.ylabel('Frecuency')
plt.xlabel('n')
plt.bar(u_dwell_mu_m69-1,co_dwell_mu_m69/sim,color='blue',alpha=0.5)
plt.plot(cn,model_dwell_mu_m69,color='red',label='$ae^{-bn}$')
plt.plot([],[],linestyle='',label='a=%.2f'%(po_dwell_mu_m69[0]))
plt.plot([],[],linestyle='',label='b=%.2f'%(po_dwell_mu_m69[1]))

plt.legend(prop={'size': 15})

a_10 = po_dwell_mu_m69[0]
b_10 = po_dwell_mu_m69[1]

#%%

a = np.array([a_0,a_1,a_2,a_3,a_4,a_5,a_6,a_7,a_8,a_9,a_10])

b = np.array([b_0,b_1,b_2,b_3,b_4,b_5,b_6,b_7,b_8,b_9,b_10])

plt.title('$J=0$, $\mu=-1.61$, ($k_{on}=1/5k_{off}$)',size=17)
plt.xlabel('k',size=14)
plt.plot([],[],linestyle='',label='$ae^{bn}$')
plt.plot(k,a,label='a',linestyle='--',marker='o')
plt.plot(k,b,label='b',linestyle='--',marker='o')

plt.legend(prop={'size': 15})
