#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 25 11:55:01 2021

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

files0=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/GTraces_Res_J=0.00  _mu= 0.69  .dat'] #All resurrection simulations
data0=[]
#files00=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/GTraces_Res_J=0.00  _mu=0.00  .dat'] #All resurrection simulations
#data00=[]
#files1=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/GTraces_Res_J=1.00  _mu= 0.00  .dat'] #All resurrection simulations
#data1=[]
files2=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/GTraces_Res_J=1.00  _mu=-0.57  .dat'] #All resurrection simulations
data2=[]
#files3=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/GTraces_Res_J=2.00  _mu= 0.00  .dat'] #All resurrection simulations
#data3=[]
files4=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/GTraces_Res_J=2.00  _mu=-1.74  .dat'] #All resurrection simulations
data4=[]
#files5=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/GTraces_Res_J=3.00  _mu= 0.00  .dat'] #All resurrection simulations
#data5=[]
files6=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/GTraces_Res_J=3.00  _mu=-2.81  .dat'] #All resurrection simulations
data6=[]
files8=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/GTraces_Res_J=4.00  _mu=-3.85  .dat'] #All resurrection simulations
data8=[]
files10=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/GTraces_Res_J=5.00  _mu=-4.79  .dat'] #All resurrection simulations
data10=[]

for data_file in files0:
    data0.append(np.loadtxt(data_file))
    
    
for data_file in files2:
    data2.append(np.loadtxt(data_file))
    
    
for data_file in files4:
    data4.append(np.loadtxt(data_file))

for data_file in files6:
    data6.append(np.loadtxt(data_file))
    
for data_file in files8:
    data8.append(np.loadtxt(data_file))
    
for data_file in files10:
    data10.append(np.loadtxt(data_file))

All_Res_J_0_mu_0=data0[0]
#All_Res_J_0_mu_m1=data00[0]
#All_Res_J_1_mu_0=data1[0]
All_Res_J_1_mu_m1=data2[0]
#All_Res_J_2_mu_0=data3[0]
All_Res_J_2_mu_m2=data4[0]
#All_Res_J_3_mu_0=data5[0]
All_Res_J_3_mu_m3=data6[0]
All_Res_J_4_mu_4=data8[0]
All_Res_J_5_mu_5=data10[0]

#%%

sim = 10000
stp = 3000

stp_4=3000
stp_5=3000
nmax=13

k = np.arange(0,11,1)

#%%

dwell_k0_J_0_mu_0 = np.zeros(sim)
dwell_k1_J_0_mu_0 = np.zeros(sim)
dwell_k2_J_0_mu_0 = np.zeros(sim)
dwell_k3_J_0_mu_0 = np.zeros(sim)
dwell_k4_J_0_mu_0 = np.zeros(sim)
dwell_k5_J_0_mu_0 = np.zeros(sim)
dwell_k6_J_0_mu_0 = np.zeros(sim)
dwell_k7_J_0_mu_0 = np.zeros(sim)
dwell_k8_J_0_mu_0 = np.zeros(sim)
dwell_k9_J_0_mu_0 = np.zeros(sim)
dwell_k10_J_0_mu_0 = np.zeros(sim)

start2J_0_0 = np.zeros(sim)
start3J_0_0 = np.zeros(sim)
start4J_0_0 = np.zeros(sim)
start5J_0_0 = np.zeros(sim)
start6J_0_0 = np.zeros(sim)
start7J_0_0 = np.zeros(sim)
start8J_0_0 = np.zeros(sim)
start9J_0_0 = np.zeros(sim)
start10J_0_0 = np.zeros(sim)

stop0J_0_0 = np.zeros(sim)
stop1J_0_0 = np.zeros(sim)
stop2J_0_0 = np.zeros(sim)
stop3J_0_0 = np.zeros(sim)
stop4J_0_0 = np.zeros(sim)
stop5J_0_0 = np.zeros(sim)
stop6J_0_0 = np.zeros(sim)
stop7J_0_0 = np.zeros(sim)
stop8J_0_0 = np.zeros(sim)
stop9J_0_0 = np.zeros(sim)


#%% Dwell time for kon = koff calculation


for i in range(sim):
    for j in range(stp*i,stp*(i+1)-1):
        if All_Res_J_0_mu_0[j,2]*nmax==0.:
            dwell_k0_J_0_mu_0[i] = dwell_k0_J_0_mu_0[i] + 1.
        else:
            stop0J_0_0[i] = j
            break

for i in range(sim):     
    for j in range(int(stop0J_0_0[i]),stp*(i+1)-1):
        if All_Res_J_0_mu_0[j,2]*nmax==1.:
            dwell_k1_J_0_mu_0[i] = dwell_k1_J_0_mu_0[i] + 1.
        else:
            stop1J_0_0[i] = j
            break
        
for i in range(sim):
    for j in range(int(stop1J_0_0[i]),stp*(i+1)-1):
        if All_Res_J_0_mu_0[j,2]*nmax==2.:
            start2J_0_0[i] = j           
            break        
for i in range(sim):
    for j in range(int(start2J_0_0[i]),stp*(i+1)-1):
            if All_Res_J_0_mu_0[j,2]*nmax==2.:
                dwell_k2_J_0_mu_0[i] = dwell_k2_J_0_mu_0[i] + 1.
            else:
                stop2J_0_0[i] = j
                break
            
for i in range(sim):
    for j in range(int(stop2J_0_0[i]),stp*(i+1)-1):
        if All_Res_J_0_mu_0[j,2]*nmax==3.:
            start3J_0_0[i] = j
            break
for i in range(sim):
    for j in range(int(start3J_0_0[i]),stp*(i+1)-1):
        if All_Res_J_0_mu_0[j,2]*nmax==3.:
            dwell_k3_J_0_mu_0[i] = dwell_k3_J_0_mu_0[i] + 1.
        else:
            stop3J_0_0[i] = j
            break

for i in range(sim):
    for j in range(int(stop3J_0_0[i]),stp*(i+1)-1):
        if All_Res_J_0_mu_0[j,2]*nmax==4.:
            start4J_0_0[i] = j
            break
for i in range(sim):
    for j in range(int(start4J_0_0[i]),stp*(i+1)-1):
        if All_Res_J_0_mu_0[j,2]*nmax==4.:
            dwell_k4_J_0_mu_0[i] = dwell_k4_J_0_mu_0[i] + 1.
        else:
            stop4J_0_0[i] = j
            break
        
for i in range(sim):
    for j in range(int(stop4J_0_0[i]),stp*(i+1)-1):
        if All_Res_J_0_mu_0[j,2]*nmax==5.:
            start5J_0_0[i] = j
            break
for i in range(sim):
    for j in range(int(start5J_0_0[i]),stp*(i+1)-1):
        if All_Res_J_0_mu_0[j,2]*nmax==5.:
            dwell_k5_J_0_mu_0[i] = dwell_k5_J_0_mu_0[i] + 1.
        else:
            stop5J_0_0[i] = j
            break
        
for i in range(sim):
    for j in range(int(stop5J_0_0[i]),stp*(i+1)-1):
        if All_Res_J_0_mu_0[j,2]*nmax==6.:
            start6J_0_0[i] = j
            break
for i in range(sim):
    for j in range(int(start6J_0_0[i]),stp*(i+1)-1):
        if All_Res_J_0_mu_0[j,2]*nmax==6.:
            dwell_k6_J_0_mu_0[i] = dwell_k6_J_0_mu_0[i] + 1.
        else:
            stop6J_0_0[i] = j
            break

for i in range(sim):
    for j in range(int(stop6J_0_0[i]),stp*(i+1)-1):
        if All_Res_J_0_mu_0[j,2]*nmax==7.:
            start7J_0_0[i] = j
            break
for i in range(sim):
    for j in range(int(start7J_0_0[i]),stp*(i+1)-1):
        if All_Res_J_0_mu_0[j,2]*nmax==7.:
            dwell_k7_J_0_mu_0[i] = dwell_k7_J_0_mu_0[i] + 1.
        else:
            stop7J_0_0[i] = j
            break
        
for i in range(sim):
    for j in range(int(stop7J_0_0[i]),stp*(i+1)-1):
        if All_Res_J_0_mu_0[j,2]*nmax==8.:
            start8J_0_0[i] = j
            break
for i in range(sim):
    for j in range(int(start8J_0_0[i]),stp*(i+1)-1):
        if All_Res_J_0_mu_0[j,2]*nmax==8.:
            dwell_k8_J_0_mu_0[i] = dwell_k8_J_0_mu_0[i] + 1.
        else:
            stop8J_0_0[i] = j
            break
   
 
for i in range(sim):
    for j in range(int(stop8J_0_0[i]),stp*(i+1)-1):
        if All_Res_J_0_mu_0[j,2]*nmax==9.:
            start9J_0_0[i] = j
            break
for i in range(sim):
    for j in range(int(start9J_0_0[i]),stp*(i+1)-1):
        if All_Res_J_0_mu_0[j,2]*nmax==9.:
            dwell_k9_J_0_mu_0[i] = dwell_k9_J_0_mu_0[i] + 1.
        else:
            stop9J_0_0[i] = j
            break

        
for i in range(sim):
    for j in range(int(stop9J_0_0[i]),stp*(i+1)-1):
        if All_Res_J_0_mu_0[j,2]*nmax==10.:
            start10J_0_0[i] = j
            break
for i in range(sim):
    for j in range(int(start10J_0_0[i]),stp*(i+1)-1):
        if All_Res_J_0_mu_0[j,2]*nmax==10.:
            dwell_k10_J_0_mu_0[i] = dwell_k10_J_0_mu_0[i] + 1.
        else:
            break


#%% J=1, mu = -1

dwell_k0_J_1_mu_m1 = np.zeros(sim)
dwell_k1_J_1_mu_m1 = np.zeros(sim)
dwell_k2_J_1_mu_m1 = np.zeros(sim)
dwell_k3_J_1_mu_m1 = np.zeros(sim)
dwell_k4_J_1_mu_m1 = np.zeros(sim)
dwell_k5_J_1_mu_m1 = np.zeros(sim)
dwell_k6_J_1_mu_m1 = np.zeros(sim)
dwell_k7_J_1_mu_m1 = np.zeros(sim)
dwell_k8_J_1_mu_m1 = np.zeros(sim)
dwell_k9_J_1_mu_m1 = np.zeros(sim)
dwell_k10_J_1_mu_m1 = np.zeros(sim)

start2_m1 = np.zeros(sim)
start3_m1 = np.zeros(sim)
start4_m1 = np.zeros(sim)
start5_m1 = np.zeros(sim)
start6_m1 = np.zeros(sim)
start7_m1 = np.zeros(sim)
start8_m1 = np.zeros(sim)
start9_m1 = np.zeros(sim)
start10_m1 = np.zeros(sim)



#%% J=1, mu=-1


for i in range(sim):
    for j in range(stp*i,stp*(i+1)-1):
        if All_Res_J_1_mu_m1[j,2]*nmax<1.:
            dwell_k0_J_1_mu_m1[i] = dwell_k0_J_1_mu_m1[i] + 1.
        else:
            break

for i in range(sim):     
    for j in range(int(stp*i+dwell_k0_J_1_mu_m1[i]),stp*(i+1)-1):
        if All_Res_J_1_mu_m1[j,2]*nmax==1.:
            dwell_k1_J_1_mu_m1[i] = dwell_k1_J_1_mu_m1[i] + 1.
        else:
            break
        
for i in range(sim):
    for j in range(int(stp*i+dwell_k1_J_1_mu_m1[i]),stp*(i+1)-1):
        if All_Res_J_1_mu_m1[j,2]*nmax==2.:
            start2_m1[i] = j           
            break        
for i in range(sim):
    for j in range(int(start2_m1[i]),stp*(i+1)-1):
            if All_Res_J_1_mu_m1[j,2]*nmax==2.:
                dwell_k2_J_1_mu_m1[i] = dwell_k2_J_1_mu_m1[i] + 1.
            else:
                break
            
for i in range(sim):
    for j in range(int(stp*i+dwell_k2_J_1_mu_m1[i]),stp*(i+1)-1):
        if All_Res_J_1_mu_m1[j,2]*nmax==3.:
            start3_m1[i] = j
            break
for i in range(sim):
    for j in range(int(start3_m1[i]),stp*(i+1)-1):
        if All_Res_J_1_mu_m1[j,2]*nmax==3.:
            dwell_k3_J_1_mu_m1[i] = dwell_k3_J_1_mu_m1[i] + 1.
        else:
            break

for i in range(sim):
    for j in range(int(stp*i+dwell_k3_J_1_mu_m1[i]),stp*(i+1)-1):
        if All_Res_J_1_mu_m1[j,2]*nmax==4.:
            start4_m1[i] = j
            break
for i in range(sim):
    for j in range(int(start4_m1[i]),stp*(i+1)-1):
        if All_Res_J_1_mu_m1[j,2]*nmax==4.:
            dwell_k4_J_1_mu_m1[i] = dwell_k4_J_1_mu_m1[i] + 1.
        else:
            break
        
for i in range(sim):
    for j in range(int(stp*i+dwell_k4_J_1_mu_m1[i]),stp*(i+1)-1):
        if All_Res_J_1_mu_m1[j,2]*nmax==5.:
            start5_m1[i] = j
            break
for i in range(sim):
    for j in range(int(start5_m1[i]),stp*(i+1)-1):
        if All_Res_J_1_mu_m1[j,2]*nmax==5.:
            dwell_k5_J_1_mu_m1[i] = dwell_k5_J_1_mu_m1[i] + 1.
        else:
            break
        
for i in range(sim):
    for j in range(int(stp*i+dwell_k5_J_1_mu_m1[i]),stp*(i+1)-1):
        if All_Res_J_1_mu_m1[j,2]*nmax==6.:
            start6_m1[i] = j
            break
for i in range(sim):
    for j in range(int(start6_m1[i]),stp*(i+1)-1):
        if All_Res_J_1_mu_m1[j,2]*nmax==6.:
            dwell_k6_J_1_mu_m1[i] = dwell_k6_J_1_mu_m1[i] + 1.
        else:
            break
        
for i in range(sim):
    for j in range(int(stp*i+dwell_k6_J_1_mu_m1[i]),stp*(i+1)-1):
        if All_Res_J_1_mu_m1[j,2]*nmax==7.:
            start7_m1[i] = j
            break
for i in range(sim):
    for j in range(int(start7_m1[i]),stp*(i+1)-1):
        if All_Res_J_1_mu_m1[j,2]*nmax==7.:
            dwell_k7_J_1_mu_m1[i] = dwell_k7_J_1_mu_m1[i] + 1.
        else:
            break
        
for i in range(sim):
    for j in range(int(stp*i+dwell_k7_J_1_mu_m1[i]),stp*(i+1)-1):
        if All_Res_J_1_mu_m1[j,2]*nmax==8.:
            start8_m1[i] = j
            break
for i in range(sim):
    for j in range(int(start8_m1[i]),stp*(i+1)-1):
        if All_Res_J_1_mu_m1[j,2]*nmax==8.:
            dwell_k8_J_1_mu_m1[i] = dwell_k8_J_1_mu_m1[i] + 1.
        else:
            break        
        
for i in range(sim):
    for j in range(int(stp*i+dwell_k8_J_1_mu_m1[i]),stp*(i+1)-1):
        if All_Res_J_1_mu_m1[j,2]*nmax==9.:
            start9_m1[i] = j
            break
for i in range(sim):
    for j in range(int(start9_m1[i]),stp*(i+1)-1):
        if All_Res_J_1_mu_m1[j,2]*nmax==9.:
            dwell_k9_J_1_mu_m1[i] = dwell_k9_J_1_mu_m1[i] + 1.
        else:
            break        

for i in range(sim):
    for j in range(int(stp*i+dwell_k9_J_1_mu_m1[i]),stp*(i+1)-1):
        if All_Res_J_1_mu_m1[j,2]*nmax==10.:
            start10_m1[i] = j
            break
for i in range(sim):
    for j in range(int(start10_m1[i]),stp*(i+1)-1):
        if All_Res_J_1_mu_m1[j,2]*nmax==10.:
            dwell_k10_J_1_mu_m1[i] = dwell_k10_J_1_mu_m1[i] + 1.
        else:
            break
        

        
#%% J=2 mu=-2

dwell_k0_J_2_mu_m2 = np.zeros(sim)
dwell_k1_J_2_mu_m2 = np.zeros(sim)
dwell_k2_J_2_mu_m2 = np.zeros(sim)
dwell_k3_J_2_mu_m2 = np.zeros(sim)
dwell_k4_J_2_mu_m2 = np.zeros(sim)
dwell_k5_J_2_mu_m2 = np.zeros(sim)
dwell_k6_J_2_mu_m2 = np.zeros(sim)
dwell_k7_J_2_mu_m2 = np.zeros(sim)
dwell_k8_J_2_mu_m2 = np.zeros(sim)
dwell_k9_J_2_mu_m2 = np.zeros(sim)
dwell_k10_J_2_mu_m2 = np.zeros(sim)

start2_J_2_m2 = np.zeros(sim)
start3_J_2_m2 = np.zeros(sim)
start4_J_2_m2 = np.zeros(sim)
start5_J_2_m2 = np.zeros(sim)
start6_J_2_m2 = np.zeros(sim)
start7_J_2_m2 = np.zeros(sim)
start8_J_2_m2 = np.zeros(sim)
start9_J_2_m2 = np.zeros(sim)
start10_J_2_m2 = np.zeros(sim)



#%% J=2, mu=-2


for i in range(sim):
    for j in range(stp*i,stp*(i+1)-1):
        if All_Res_J_2_mu_m2[j,2]*nmax<1.:
            dwell_k0_J_2_mu_m2[i] = dwell_k0_J_2_mu_m2[i] + 1.
        else:
            break

for i in range(sim):     
    for j in range(int(stp*i+dwell_k0_J_2_mu_m2[i]),stp*(i+1)-1):
        if All_Res_J_2_mu_m2[j,2]*nmax==1.:
            dwell_k1_J_2_mu_m2[i] = dwell_k1_J_2_mu_m2[i] + 1.
        else:
            break
        
for i in range(sim):
    for j in range(int(stp*i+dwell_k1_J_2_mu_m2[i]),stp*(i+1)-1):
        if All_Res_J_2_mu_m2[j,2]*nmax==2.:
            start2_J_2_m2[i] = j           
            break        
for i in range(sim):
    for j in range(int(start2_J_2_m2[i]),stp*(i+1)-1):
            if All_Res_J_2_mu_m2[j,2]*nmax==2.:
                dwell_k2_J_2_mu_m2[i] = dwell_k2_J_2_mu_m2[i] + 1.
            else:
                break
            
for i in range(sim):
    for j in range(int(stp*i+dwell_k2_J_2_mu_m2[i]),stp*(i+1)-1):
        if All_Res_J_2_mu_m2[j,2]*nmax==3.:
            start3_J_2_m2[i] = j
            break
for i in range(sim):
    for j in range(int(start3_J_2_m2[i]),stp*(i+1)-1):
        if All_Res_J_2_mu_m2[j,2]*nmax==3.:
            dwell_k3_J_2_mu_m2[i] = dwell_k3_J_2_mu_m2[i] + 1.
        else:
            break

for i in range(sim):
    for j in range(int(stp*i+dwell_k3_J_2_mu_m2[i]),stp*(i+1)-1):
        if All_Res_J_2_mu_m2[j,2]*nmax==4.:
            start4_J_2_m2[i] = j
            break
for i in range(sim):
    for j in range(int(start4_J_2_m2[i]),stp*(i+1)-1):
        if All_Res_J_2_mu_m2[j,2]*nmax==4.:
            dwell_k4_J_2_mu_m2[i] = dwell_k4_J_2_mu_m2[i] + 1.
        else:
            break
        
for i in range(sim):
    for j in range(int(stp*i+dwell_k4_J_2_mu_m2[i]),stp*(i+1)-1):
        if All_Res_J_2_mu_m2[j,2]*nmax==5.:
            start5_J_2_m2[i] = j
            break
for i in range(sim):
    for j in range(int(start5_J_2_m2[i]),stp*(i+1)-1):
        if All_Res_J_2_mu_m2[j,2]*nmax==5.:
            dwell_k5_J_2_mu_m2[i] = dwell_k5_J_2_mu_m2[i] + 1.
        else:
            break
        
for i in range(sim):
    for j in range(int(stp*i+dwell_k5_J_2_mu_m2[i]),stp*(i+1)-1):
        if All_Res_J_2_mu_m2[j,2]*nmax==6.:
            start6_J_2_m2[i] = j
            break
for i in range(sim):
    for j in range(int(start6_J_2_m2[i]),stp*(i+1)-1):
        if All_Res_J_2_mu_m2[j,2]*nmax==6.:
            dwell_k6_J_2_mu_m2[i] = dwell_k6_J_2_mu_m2[i] + 1.
        else:
            break
        
for i in range(sim):
    for j in range(int(stp*i+dwell_k6_J_2_mu_m2[i]),stp*(i+1)-1):
        if All_Res_J_2_mu_m2[j,2]*nmax==7.:
            start7_J_2_m2[i] = j
            break
for i in range(sim):
    for j in range(int(start7_J_2_m2[i]),stp*(i+1)-1):
        if All_Res_J_2_mu_m2[j,2]*nmax==7.:
            dwell_k7_J_2_mu_m2[i] = dwell_k7_J_2_mu_m2[i] + 1.
        else:
            break
        
for i in range(sim):
    for j in range(int(stp*i+dwell_k7_J_2_mu_m2[i]),stp*(i+1)-1):
        if All_Res_J_2_mu_m2[j,2]*nmax==8.:
            start8_J_2_m2[i] = j
            break
for i in range(sim):
    for j in range(int(start8_J_2_m2[i]),stp*(i+1)-1):
        if All_Res_J_2_mu_m2[j,2]*nmax==8.:
            dwell_k8_J_2_mu_m2[i] = dwell_k8_J_2_mu_m2[i] + 1.
        else:
            break        
        
for i in range(sim):
    for j in range(int(stp*i+dwell_k8_J_2_mu_m2[i]),stp*(i+1)-1):
        if All_Res_J_2_mu_m2[j,2]*nmax==9.:
            start9_J_2_m2[i] = j
            break
for i in range(sim):
    for j in range(int(start9_J_2_m2[i]),stp*(i+1)-1):
        if All_Res_J_2_mu_m2[j,2]*nmax==9.:
            dwell_k9_J_2_mu_m2[i] = dwell_k9_J_2_mu_m2[i] + 1.
        else:
            break        

for i in range(sim):
    for j in range(int(stp*i+dwell_k9_J_2_mu_m2[i]),stp*(i+1)-1):
        if All_Res_J_2_mu_m2[j,2]*nmax==10.:
            start10_J_2_m2[i] = j
            break
for i in range(sim):
    for j in range(int(start10_J_2_m2[i]),stp*(i+1)-1):
        if All_Res_J_2_mu_m2[j,2]*nmax==10.:
            dwell_k10_J_2_mu_m2[i] = dwell_k10_J_2_mu_m2[i] + 1.
        else:
            break



#%%

dwell_k0_J_3_mu_m3 = np.zeros(sim)
dwell_k1_J_3_mu_m3 = np.zeros(sim)
dwell_k2_J_3_mu_m3 = np.zeros(sim)
dwell_k3_J_3_mu_m3 = np.zeros(sim)
dwell_k4_J_3_mu_m3 = np.zeros(sim)
dwell_k5_J_3_mu_m3 = np.zeros(sim)
dwell_k6_J_3_mu_m3 = np.zeros(sim)
dwell_k7_J_3_mu_m3 = np.zeros(sim)
dwell_k8_J_3_mu_m3 = np.zeros(sim)
dwell_k9_J_3_mu_m3 = np.zeros(sim)
dwell_k10_J_3_mu_m3 = np.zeros(sim)

start2_J_3_m3 = np.zeros(sim)
start3_J_3_m3 = np.zeros(sim)
start4_J_3_m3 = np.zeros(sim)
start5_J_3_m3 = np.zeros(sim)
start6_J_3_m3 = np.zeros(sim)
start7_J_3_m3 = np.zeros(sim)
start8_J_3_m3 = np.zeros(sim)
start9_J_3_m3 = np.zeros(sim)
start10_J_3_m3 = np.zeros(sim)



#%% J=2, mu=-2


for i in range(sim):
    for j in range(stp*i,stp*(i+1)-1):
        if All_Res_J_3_mu_m3[j,2]*nmax<1.:
            dwell_k0_J_3_mu_m3[i] = dwell_k0_J_3_mu_m3[i] + 1.
        else:
            break

for i in range(sim):     
    for j in range(int(stp*i+dwell_k0_J_3_mu_m3[i]),stp*(i+1)-1):
        if All_Res_J_3_mu_m3[j,2]*nmax==1.:
            dwell_k1_J_3_mu_m3[i] = dwell_k1_J_3_mu_m3[i] + 1.
        else:
            break
        
for i in range(sim):
    for j in range(int(stp*i+dwell_k1_J_3_mu_m3[i]),stp*(i+1)-1):
        if All_Res_J_3_mu_m3[j,2]*nmax==2.:
            start2_J_3_m3[i] = j           
            break        
for i in range(sim):
    for j in range(int(start2_J_3_m3[i]),stp*(i+1)-1):
            if All_Res_J_3_mu_m3[j,2]*nmax==2.:
                dwell_k2_J_3_mu_m3[i] = dwell_k2_J_3_mu_m3[i] + 1.
            else:
                break
            
for i in range(sim):
    for j in range(int(stp*i+dwell_k2_J_3_mu_m3[i]),stp*(i+1)-1):
        if All_Res_J_3_mu_m3[j,2]*nmax==3.:
            start3_J_3_m3[i] = j
            break
for i in range(sim):
    for j in range(int(start3_J_3_m3[i]),stp*(i+1)-1):
        if All_Res_J_3_mu_m3[j,2]*nmax==3.:
            dwell_k3_J_3_mu_m3[i] = dwell_k3_J_3_mu_m3[i] + 1.
        else:
            break

for i in range(sim):
    for j in range(int(stp*i+dwell_k3_J_3_mu_m3[i]),stp*(i+1)-1):
        if All_Res_J_3_mu_m3[j,2]*nmax==4.:
            start4_J_3_m3[i] = j
            break
for i in range(sim):
    for j in range(int(start4_J_3_m3[i]),stp*(i+1)-1):
        if All_Res_J_3_mu_m3[j,2]*nmax==4.:
            dwell_k4_J_3_mu_m3[i] = dwell_k4_J_3_mu_m3[i] + 1.
        else:
            break
        
for i in range(sim):
    for j in range(int(stp*i+dwell_k4_J_3_mu_m3[i]),stp*(i+1)-1):
        if All_Res_J_3_mu_m3[j,2]*nmax==5.:
            start5_J_3_m3[i] = j
            break
for i in range(sim):
    for j in range(int(start5_J_3_m3[i]),stp*(i+1)-1):
        if All_Res_J_3_mu_m3[j,2]*nmax==5.:
            dwell_k5_J_3_mu_m3[i] = dwell_k5_J_3_mu_m3[i] + 1.
        else:
            break
        
for i in range(sim):
    for j in range(int(stp*i+dwell_k5_J_3_mu_m3[i]),stp*(i+1)-1):
        if All_Res_J_3_mu_m3[j,2]*nmax==6.:
            start6_J_3_m3[i] = j
            break
for i in range(sim):
    for j in range(int(start6_J_3_m3[i]),stp*(i+1)-1):
        if All_Res_J_3_mu_m3[j,2]*nmax==6.:
            dwell_k6_J_3_mu_m3[i] = dwell_k6_J_3_mu_m3[i] + 1.
        else:
            break
        
for i in range(sim):
    for j in range(int(stp*i+dwell_k6_J_3_mu_m3[i]),stp*(i+1)-1):
        if All_Res_J_3_mu_m3[j,2]*nmax==7.:
            start7_J_3_m3[i] = j
            break
for i in range(sim):
    for j in range(int(start7_J_3_m3[i]),stp*(i+1)-1):
        if All_Res_J_3_mu_m3[j,2]*nmax==7.:
            dwell_k7_J_3_mu_m3[i] = dwell_k7_J_3_mu_m3[i] + 1.
        else:
            break
        
for i in range(sim):
    for j in range(int(stp*i+dwell_k7_J_3_mu_m3[i]),stp*(i+1)-1):
        if All_Res_J_3_mu_m3[j,2]*nmax==8.:
            start8_J_3_m3[i] = j
            break
for i in range(sim):
    for j in range(int(start8_J_3_m3[i]),stp*(i+1)-1):
        if All_Res_J_3_mu_m3[j,2]*nmax==8.:
            dwell_k8_J_3_mu_m3[i] = dwell_k8_J_3_mu_m3[i] + 1.
        else:
            break        
        
for i in range(sim):
    for j in range(int(stp*i+dwell_k8_J_3_mu_m3[i]),stp*(i+1)-1):
        if All_Res_J_3_mu_m3[j,2]*nmax==9.:
            start9_J_3_m3[i] = j
            break
for i in range(sim):
    for j in range(int(start9_J_3_m3[i]),stp*(i+1)-1):
        if All_Res_J_3_mu_m3[j,2]*nmax==9.:
            dwell_k9_J_3_mu_m3[i] = dwell_k9_J_3_mu_m3[i] + 1.
        else:
            break        

for i in range(sim):
    for j in range(int(stp*i+dwell_k9_J_3_mu_m3[i]),stp*(i+1)-1):
        if All_Res_J_3_mu_m3[j,2]*nmax==10.:
            start10_J_3_m3[i] = j
            break
for i in range(sim):
    for j in range(int(start10_J_3_m3[i]),stp*(i+1)-1):
        if All_Res_J_3_mu_m3[j,2]*nmax==10.:
            dwell_k10_J_3_mu_m3[i] = dwell_k10_J_3_mu_m3[i] + 1.
        else:
            break
        
#%%

dwell_k0_J_4_mu_4 = np.zeros(sim)
dwell_k1_J_4_mu_4 = np.zeros(sim)
dwell_k2_J_4_mu_4 = np.zeros(sim)
dwell_k3_J_4_mu_4 = np.zeros(sim)
dwell_k4_J_4_mu_4 = np.zeros(sim)
dwell_k5_J_4_mu_4 = np.zeros(sim)
dwell_k6_J_4_mu_4 = np.zeros(sim)
dwell_k7_J_4_mu_4 = np.zeros(sim)
dwell_k8_J_4_mu_4 = np.zeros(sim)
dwell_k9_J_4_mu_4 = np.zeros(sim)
dwell_k10_J_4_mu_4 = np.zeros(sim)

start2J_4_4 = np.zeros(sim)
start3J_4_4 = np.zeros(sim)
start4J_4_4 = np.zeros(sim)
start5J_4_4 = np.zeros(sim)
start6J_4_4 = np.zeros(sim)
start7J_4_4 = np.zeros(sim)
start8J_4_4 = np.zeros(sim)
start9J_4_4 = np.zeros(sim)
start10J_4_4 = np.zeros(sim)

stop0J_4_4 = np.zeros(sim)
stop1J_4_4 = np.zeros(sim)
stop2J_4_4 = np.zeros(sim)
stop3J_4_4 = np.zeros(sim)
stop4J_4_4 = np.zeros(sim)
stop5J_4_4 = np.zeros(sim)
stop6J_4_4 = np.zeros(sim)
stop7J_4_4 = np.zeros(sim)
stop8J_4_4 = np.zeros(sim)
stop9J_4_4 = np.zeros(sim)


#%% Dwell time for kon = koff calculation


for i in range(sim):
    for j in range(stp_4*i,stp_4*(i+1)-1):
        if All_Res_J_4_mu_4[j,2]*nmax==0.:
            dwell_k0_J_4_mu_4[i] = dwell_k0_J_4_mu_4[i] + 1.
        else:
            stop0J_4_4[i] = j
            break

for i in range(sim):     
    for j in range(int(stop0J_4_4[i]),stp_4*(i+1)-1):
        if All_Res_J_4_mu_4[j,2]*nmax==1.:
            dwell_k1_J_4_mu_4[i] = dwell_k1_J_4_mu_4[i] + 1.
        else:
            stop1J_4_4[i] = j
            break
        
for i in range(sim):
    for j in range(int(stop1J_4_4[i]),stp_4*(i+1)-1):
        if All_Res_J_4_mu_4[j,2]*nmax==2.:
            start2J_4_4[i] = j           
            break        
for i in range(sim):
    for j in range(int(start2J_4_4[i]),stp_4*(i+1)-1):
            if All_Res_J_4_mu_4[j,2]*nmax==2.:
                dwell_k2_J_4_mu_4[i] = dwell_k2_J_4_mu_4[i] + 1.
            else:
                stop2J_4_4[i] = j
                break
            
for i in range(sim):
    for j in range(int(stop2J_4_4[i]),stp_4*(i+1)-1):
        if All_Res_J_4_mu_4[j,2]*nmax==3.:
            start3J_4_4[i] = j
            break
for i in range(sim):
    for j in range(int(start3J_4_4[i]),stp_4*(i+1)-1):
        if All_Res_J_4_mu_4[j,2]*nmax==3.:
            dwell_k3_J_4_mu_4[i] = dwell_k3_J_4_mu_4[i] + 1.
        else:
            stop3J_4_4[i] = j
            break

for i in range(sim):
    for j in range(int(stop3J_4_4[i]),stp_4*(i+1)-1):
        if All_Res_J_4_mu_4[j,2]*nmax==4.:
            start4J_4_4[i] = j
            break
for i in range(sim):
    for j in range(int(start4J_4_4[i]),stp_4*(i+1)-1):
        if All_Res_J_4_mu_4[j,2]*nmax==4.:
            dwell_k4_J_4_mu_4[i] = dwell_k4_J_4_mu_4[i] + 1.
        else:
            stop4J_4_4[i] = j
            break
        
for i in range(sim):
    for j in range(int(stop4J_4_4[i]),stp_4*(i+1)-1):
        if All_Res_J_4_mu_4[j,2]*nmax==5.:
            start5J_4_4[i] = j
            break
for i in range(sim):
    for j in range(int(start5J_4_4[i]),stp_4*(i+1)-1):
        if All_Res_J_4_mu_4[j,2]*nmax==5.:
            dwell_k5_J_4_mu_4[i] = dwell_k5_J_4_mu_4[i] + 1.
        else:
            stop5J_4_4[i] = j
            break
        
for i in range(sim):
    for j in range(int(stop5J_4_4[i]),stp_4*(i+1)-1):
        if All_Res_J_4_mu_4[j,2]*nmax==6.:
            start6J_4_4[i] = j
            break
for i in range(sim):
    for j in range(int(start6J_4_4[i]),stp_4*(i+1)-1):
        if All_Res_J_4_mu_4[j,2]*nmax==6.:
            dwell_k6_J_4_mu_4[i] = dwell_k6_J_4_mu_4[i] + 1.
        else:
            stop6J_4_4[i] = j
            break

for i in range(sim):
    for j in range(int(stop6J_4_4[i]),stp_4*(i+1)-1):
        if All_Res_J_4_mu_4[j,2]*nmax==7.:
            start7J_4_4[i] = j
            break
for i in range(sim):
    for j in range(int(start7J_4_4[i]),stp_4*(i+1)-1):
        if All_Res_J_4_mu_4[j,2]*nmax==7.:
            dwell_k7_J_4_mu_4[i] = dwell_k7_J_4_mu_4[i] + 1.
        else:
            stop7J_4_4[i] = j
            break
        
for i in range(sim):
    for j in range(int(stop7J_4_4[i]),stp_4*(i+1)-1):
        if All_Res_J_4_mu_4[j,2]*nmax==8.:
            start8J_4_4[i] = j
            break
for i in range(sim):
    for j in range(int(start8J_4_4[i]),stp_4*(i+1)-1):
        if All_Res_J_4_mu_4[j,2]*nmax==8.:
            dwell_k8_J_4_mu_4[i] = dwell_k8_J_4_mu_4[i] + 1.
        else:
            stop8J_4_4[i] = j
            break
   
 
for i in range(sim):
    for j in range(int(stop8J_4_4[i]),stp_4*(i+1)-1):
        if All_Res_J_4_mu_4[j,2]*nmax==9.:
            start9J_4_4[i] = j
            break
for i in range(sim):
    for j in range(int(start9J_4_4[i]),stp_4*(i+1)-1):
        if All_Res_J_4_mu_4[j,2]*nmax==9.:
            dwell_k9_J_4_mu_4[i] = dwell_k9_J_4_mu_4[i] + 1.
        else:
            stop9J_4_4[i] = j
            break

        
for i in range(sim):
    for j in range(int(stop9J_4_4[i]),stp_4*(i+1)-1):
        if All_Res_J_4_mu_4[j,2]*nmax==10.:
            start10J_4_4[i] = j
            break
for i in range(sim):
    for j in range(int(start10J_4_4[i]),stp_4*(i+1)-1):
        if All_Res_J_4_mu_4[j,2]*nmax==10.:
            dwell_k10_J_4_mu_4[i] = dwell_k10_J_4_mu_4[i] + 1.
        else:
            break

#%%

dwell_k0_J_5_mu_5 = np.zeros(sim)
dwell_k1_J_5_mu_5 = np.zeros(sim)
dwell_k2_J_5_mu_5 = np.zeros(sim)
dwell_k3_J_5_mu_5 = np.zeros(sim)
dwell_k4_J_5_mu_5 = np.zeros(sim)
dwell_k5_J_5_mu_5 = np.zeros(sim)
dwell_k6_J_5_mu_5 = np.zeros(sim)
dwell_k7_J_5_mu_5 = np.zeros(sim)
dwell_k8_J_5_mu_5 = np.zeros(sim)
dwell_k9_J_5_mu_5 = np.zeros(sim)
dwell_k10_J_5_mu_5 = np.zeros(sim)

start2J_5_5 = np.zeros(sim)
start3J_5_5 = np.zeros(sim)
start4J_5_5 = np.zeros(sim)
start5J_5_5 = np.zeros(sim)
start6J_5_5 = np.zeros(sim)
start7J_5_5 = np.zeros(sim)
start8J_5_5 = np.zeros(sim)
start9J_5_5 = np.zeros(sim)
start10J_5_5 = np.zeros(sim)

stop0J_5_5 = np.zeros(sim)
stop1J_5_5 = np.zeros(sim)
stop2J_5_5 = np.zeros(sim)
stop3J_5_5 = np.zeros(sim)
stop4J_5_5 = np.zeros(sim)
stop5J_5_5 = np.zeros(sim)
stop6J_5_5 = np.zeros(sim)
stop7J_5_5 = np.zeros(sim)
stop8J_5_5 = np.zeros(sim)
stop9J_5_5 = np.zeros(sim)


#%% Dwell time for kon = koff calculation


for i in range(sim):
    for j in range(stp_5*i,stp_5*(i+1)-1):
        if All_Res_J_5_mu_5[j,2]*nmax==0.:
            dwell_k0_J_5_mu_5[i] = dwell_k0_J_5_mu_5[i] + 1.
        else:
            stop0J_5_5[i] = j
            break

for i in range(sim):     
    for j in range(int(stop0J_5_5[i]),stp_5*(i+1)-1):
        if All_Res_J_5_mu_5[j,2]*nmax==1.:
            dwell_k1_J_5_mu_5[i] = dwell_k1_J_5_mu_5[i] + 1.
        else:
            stop1J_5_5[i] = j
            break
        
for i in range(sim):
    for j in range(int(stop1J_5_5[i]),stp_5*(i+1)-1):
        if All_Res_J_5_mu_5[j,2]*nmax==2.:
            start2J_5_5[i] = j           
            break        
for i in range(sim):
    for j in range(int(start2J_5_5[i]),stp_5*(i+1)-1):
            if All_Res_J_5_mu_5[j,2]*nmax==2.:
                dwell_k2_J_5_mu_5[i] = dwell_k2_J_5_mu_5[i] + 1.
            else:
                stop2J_5_5[i] = j
                break
            
for i in range(sim):
    for j in range(int(stop2J_5_5[i]),stp_5*(i+1)-1):
        if All_Res_J_5_mu_5[j,2]*nmax==3.:
            start3J_5_5[i] = j
            break
for i in range(sim):
    for j in range(int(start3J_5_5[i]),stp_5*(i+1)-1):
        if All_Res_J_5_mu_5[j,2]*nmax==3.:
            dwell_k3_J_5_mu_5[i] = dwell_k3_J_5_mu_5[i] + 1.
        else:
            stop3J_5_5[i] = j
            break

for i in range(sim):
    for j in range(int(stop3J_5_5[i]),stp_5*(i+1)-1):
        if All_Res_J_5_mu_5[j,2]*nmax==4.:
            start4J_5_5[i] = j
            break
for i in range(sim):
    for j in range(int(start4J_5_5[i]),stp_5*(i+1)-1):
        if All_Res_J_5_mu_5[j,2]*nmax==4.:
            dwell_k4_J_5_mu_5[i] = dwell_k4_J_5_mu_5[i] + 1.
        else:
            stop4J_5_5[i] = j
            break
        
for i in range(sim):
    for j in range(int(stop4J_5_5[i]),stp_5*(i+1)-1):
        if All_Res_J_5_mu_5[j,2]*nmax==5.:
            start5J_5_5[i] = j
            break
for i in range(sim):
    for j in range(int(start5J_5_5[i]),stp_5*(i+1)-1):
        if All_Res_J_5_mu_5[j,2]*nmax==5.:
            dwell_k5_J_5_mu_5[i] = dwell_k5_J_5_mu_5[i] + 1.
        else:
            stop5J_5_5[i] = j
            break
        
for i in range(sim):
    for j in range(int(stop5J_5_5[i]),stp_5*(i+1)-1):
        if All_Res_J_5_mu_5[j,2]*nmax==6.:
            start6J_5_5[i] = j
            break
for i in range(sim):
    for j in range(int(start6J_5_5[i]),stp_5*(i+1)-1):
        if All_Res_J_5_mu_5[j,2]*nmax==6.:
            dwell_k6_J_5_mu_5[i] = dwell_k6_J_5_mu_5[i] + 1.
        else:
            stop6J_5_5[i] = j
            break

for i in range(sim):
    for j in range(int(stop6J_5_5[i]),stp_5*(i+1)-1):
        if All_Res_J_5_mu_5[j,2]*nmax==7.:
            start7J_5_5[i] = j
            break
for i in range(sim):
    for j in range(int(start7J_5_5[i]),stp_5*(i+1)-1):
        if All_Res_J_5_mu_5[j,2]*nmax==7.:
            dwell_k7_J_5_mu_5[i] = dwell_k7_J_5_mu_5[i] + 1.
        else:
            stop7J_5_5[i] = j
            break
        
for i in range(sim):
    for j in range(int(stop7J_5_5[i]),stp_5*(i+1)-1):
        if All_Res_J_5_mu_5[j,2]*nmax==8.:
            start8J_5_5[i] = j
            break
for i in range(sim):
    for j in range(int(start8J_5_5[i]),stp_5*(i+1)-1):
        if All_Res_J_5_mu_5[j,2]*nmax==8.:
            dwell_k8_J_5_mu_5[i] = dwell_k8_J_5_mu_5[i] + 1.
        else:
            stop8J_5_5[i] = j
            break
   
 
for i in range(sim):
    for j in range(int(stop8J_5_5[i]),stp_5*(i+1)-1):
        if All_Res_J_5_mu_5[j,2]*nmax==9.:
            start9J_5_5[i] = j
            break
for i in range(sim):
    for j in range(int(start9J_5_5[i]),stp_5*(i+1)-1):
        if All_Res_J_5_mu_5[j,2]*nmax==9.:
            dwell_k9_J_5_mu_5[i] = dwell_k9_J_5_mu_5[i] + 1.
        else:
            stop9J_5_5[i] = j
            break

        
for i in range(sim):
    for j in range(int(stop9J_5_5[i]),stp_5*(i+1)-1):
        if All_Res_J_5_mu_5[j,2]*nmax==10.:
            start10J_5_5[i] = j
            break
for i in range(sim):
    for j in range(int(start10J_5_5[i]),stp_5*(i+1)-1):
        if All_Res_J_5_mu_5[j,2]*nmax==10.:
            dwell_k10_J_5_mu_5[i] = dwell_k10_J_5_mu_5[i] + 1.
        else:
            break


#%%

dwell_avg_k0_J_0_mu_0 = sum(dwell_k0_J_0_mu_0)/sim
dwell_avg2_k0_J_0_mu_0 = sum(dwell_k0_J_0_mu_0**2.)/sim
dwell_sd_k0_J_0_mu_0 = np.sqrt(dwell_avg2_k0_J_0_mu_0-dwell_avg_k0_J_0_mu_0**2)

dwell_avg_k1_J_0_mu_0 = sum(dwell_k1_J_0_mu_0)/sim
dwell_avg2_k1_J_0_mu_0 = sum(dwell_k1_J_0_mu_0**2.)/sim
dwell_sd_k1_J_0_mu_0 = np.sqrt(dwell_avg2_k1_J_0_mu_0-dwell_avg_k1_J_0_mu_0**2)

dwell_avg_k2_J_0_mu_0 = sum(dwell_k2_J_0_mu_0)/sim
dwell_avg2_k2_J_0_mu_0 = sum(dwell_k2_J_0_mu_0**2.)/sim
dwell_sd_k2_J_0_mu_0 = np.sqrt(dwell_avg2_k2_J_0_mu_0-dwell_avg_k2_J_0_mu_0**2)   

dwell_avg_k3_J_0_mu_0 = sum(dwell_k3_J_0_mu_0)/sim
dwell_avg2_k3_J_0_mu_0 = sum(dwell_k3_J_0_mu_0**2.)/sim
dwell_sd_k3_J_0_mu_0 = np.sqrt(dwell_avg2_k3_J_0_mu_0-dwell_avg_k3_J_0_mu_0**2)    

dwell_avg_k4_J_0_mu_0 = sum(dwell_k4_J_0_mu_0)/sim
dwell_avg2_k4_J_0_mu_0 = sum(dwell_k4_J_0_mu_0**2.)/sim
dwell_sd_k4_J_0_mu_0 = np.sqrt(dwell_avg2_k4_J_0_mu_0-dwell_avg_k4_J_0_mu_0**2) 

dwell_avg_k5_J_0_mu_0 = sum(dwell_k5_J_0_mu_0)/sim
dwell_avg2_k5_J_0_mu_0 = sum(dwell_k5_J_0_mu_0**2.)/sim
dwell_sd_k5_J_0_mu_0 = np.sqrt(dwell_avg2_k5_J_0_mu_0-dwell_avg_k5_J_0_mu_0**2)

dwell_avg_k6_J_0_mu_0 = sum(dwell_k6_J_0_mu_0)/sim
dwell_avg2_k6_J_0_mu_0 = sum(dwell_k6_J_0_mu_0**2.)/sim
dwell_sd_k6_J_0_mu_0 = np.sqrt(dwell_avg2_k6_J_0_mu_0-dwell_avg_k6_J_0_mu_0**2) 

dwell_avg_k7_J_0_mu_0 = sum(dwell_k7_J_0_mu_0)/sim
dwell_avg2_k7_J_0_mu_0 = sum(dwell_k7_J_0_mu_0**2.)/sim
dwell_sd_k7_J_0_mu_0 = np.sqrt(dwell_avg2_k7_J_0_mu_0-dwell_avg_k7_J_0_mu_0**2)

dwell_avg_k8_J_0_mu_0 = sum(dwell_k8_J_0_mu_0)/sim
dwell_avg2_k8_J_0_mu_0 = sum(dwell_k8_J_0_mu_0**2.)/sim
dwell_sd_k8_J_0_mu_0 = np.sqrt(dwell_avg2_k8_J_0_mu_0-dwell_avg_k8_J_0_mu_0**2)  

dwell_avg_k9_J_0_mu_0 = sum(dwell_k9_J_0_mu_0)/sim
dwell_avg2_k9_J_0_mu_0 = sum(dwell_k9_J_0_mu_0**2.)/sim
dwell_sd_k9_J_0_mu_0 = np.sqrt(dwell_avg2_k9_J_0_mu_0-dwell_avg_k9_J_0_mu_0**2)

dwell_avg_k10_J_0_mu_0 = sum(dwell_k10_J_0_mu_0)/sim
dwell_avg2_k10_J_0_mu_0 = sum(dwell_k10_J_0_mu_0**2.)/sim
dwell_sd_k10_J_0_mu_0 = np.sqrt(dwell_avg2_k10_J_0_mu_0-dwell_avg_k10_J_0_mu_0**2)  



#%% Average and sd of kon = 10koff

dwell_avg_k0_J_1_mu_m1 = sum(dwell_k0_J_1_mu_m1)/sim
dwell_avg2_k0_J_1_mu_m1 = sum(dwell_k0_J_1_mu_m1**2.)/sim
dwell_sd_k0_J_1_mu_m1 = np.sqrt(dwell_avg2_k0_J_1_mu_m1-dwell_avg_k0_J_1_mu_m1**2)

dwell_avg_k1_J_1_mu_m1 = sum(dwell_k1_J_1_mu_m1)/sim
dwell_avg2_k1_J_1_mu_m1 = sum(dwell_k1_J_1_mu_m1**2.)/sim
dwell_sd_k1_J_1_mu_m1 = np.sqrt(dwell_avg2_k1_J_1_mu_m1-dwell_avg_k1_J_1_mu_m1**2)

dwell_avg_k2_J_1_mu_m1 = sum(dwell_k2_J_1_mu_m1)/sim
dwell_avg2_k2_J_1_mu_m1 = sum(dwell_k2_J_1_mu_m1**2.)/sim
dwell_sd_k2_J_1_mu_m1 = np.sqrt(dwell_avg2_k2_J_1_mu_m1-dwell_avg_k2_J_1_mu_m1**2)   

dwell_avg_k3_J_1_mu_m1 = sum(dwell_k3_J_1_mu_m1)/sim
dwell_avg2_k3_J_1_mu_m1 = sum(dwell_k3_J_1_mu_m1**2.)/sim
dwell_sd_k3_J_1_mu_m1 = np.sqrt(dwell_avg2_k3_J_1_mu_m1-dwell_avg_k3_J_1_mu_m1**2)    

dwell_avg_k4_J_1_mu_m1 = sum(dwell_k4_J_1_mu_m1)/sim
dwell_avg2_k4_J_1_mu_m1 = sum(dwell_k4_J_1_mu_m1**2.)/sim
dwell_sd_k4_J_1_mu_m1 = np.sqrt(dwell_avg2_k4_J_1_mu_m1-dwell_avg_k4_J_1_mu_m1**2) 

dwell_avg_k5_J_1_mu_m1 = sum(dwell_k5_J_1_mu_m1)/sim
dwell_avg2_k5_J_1_mu_m1 = sum(dwell_k5_J_1_mu_m1**2.)/sim
dwell_sd_k5_J_1_mu_m1 = np.sqrt(dwell_avg2_k5_J_1_mu_m1-dwell_avg_k5_J_1_mu_m1**2)

dwell_avg_k6_J_1_mu_m1 = sum(dwell_k6_J_1_mu_m1)/sim
dwell_avg2_k6_J_1_mu_m1 = sum(dwell_k6_J_1_mu_m1**2.)/sim
dwell_sd_k6_J_1_mu_m1 = np.sqrt(dwell_avg2_k6_J_1_mu_m1-dwell_avg_k6_J_1_mu_m1**2)

dwell_avg_k7_J_1_mu_m1 = sum(dwell_k7_J_1_mu_m1)/sim
dwell_avg2_k7_J_1_mu_m1 = sum(dwell_k7_J_1_mu_m1**2.)/sim
dwell_sd_k7_J_1_mu_m1 = np.sqrt(dwell_avg2_k7_J_1_mu_m1-dwell_avg_k7_J_1_mu_m1**2)

dwell_avg_k8_J_1_mu_m1 = sum(dwell_k8_J_1_mu_m1)/sim
dwell_avg2_k8_J_1_mu_m1 = sum(dwell_k8_J_1_mu_m1**2.)/sim
dwell_sd_k8_J_1_mu_m1 = np.sqrt(dwell_avg2_k8_J_1_mu_m1-dwell_avg_k8_J_1_mu_m1**2)    

dwell_avg_k9_J_1_mu_m1 = sum(dwell_k9_J_1_mu_m1)/sim
dwell_avg2_k9_J_1_mu_m1 = sum(dwell_k9_J_1_mu_m1**2.)/sim
dwell_sd_k9_J_1_mu_m1 = np.sqrt(dwell_avg2_k9_J_1_mu_m1-dwell_avg_k9_J_1_mu_m1**2)

dwell_avg_k10_J_1_mu_m1 = sum(dwell_k10_J_1_mu_m1)/sim
dwell_avg2_k10_J_1_mu_m1 = sum(dwell_k10_J_1_mu_m1**2.)/sim
dwell_sd_k10_J_1_mu_m1 = np.sqrt(dwell_avg2_k10_J_1_mu_m1-dwell_avg_k10_J_1_mu_m1**2)



#%%

dwell_avg_k0_J_2_mu_m2 = sum(dwell_k0_J_2_mu_m2)/sim
dwell_avg2_k0_J_2_mu_m2 = sum(dwell_k0_J_2_mu_m2**2.)/sim
dwell_sd_k0_J_2_mu_m2 = np.sqrt(dwell_avg2_k0_J_2_mu_m2-dwell_avg_k0_J_2_mu_m2**2)

dwell_avg_k1_J_2_mu_m2 = sum(dwell_k1_J_2_mu_m2)/sim
dwell_avg2_k1_J_2_mu_m2 = sum(dwell_k1_J_2_mu_m2**2.)/sim
dwell_sd_k1_J_2_mu_m2 = np.sqrt(dwell_avg2_k1_J_2_mu_m2-dwell_avg_k1_J_2_mu_m2**2)

dwell_avg_k2_J_2_mu_m2 = sum(dwell_k2_J_2_mu_m2)/sim
dwell_avg2_k2_J_2_mu_m2 = sum(dwell_k2_J_2_mu_m2**2.)/sim
dwell_sd_k2_J_2_mu_m2 = np.sqrt(dwell_avg2_k2_J_2_mu_m2-dwell_avg_k2_J_2_mu_m2**2)   

dwell_avg_k3_J_2_mu_m2 = sum(dwell_k3_J_2_mu_m2)/sim
dwell_avg2_k3_J_2_mu_m2 = sum(dwell_k3_J_2_mu_m2**2.)/sim
dwell_sd_k3_J_2_mu_m2 = np.sqrt(dwell_avg2_k3_J_2_mu_m2-dwell_avg_k3_J_2_mu_m2**2)    

dwell_avg_k4_J_2_mu_m2 = sum(dwell_k4_J_2_mu_m2)/sim
dwell_avg2_k4_J_2_mu_m2 = sum(dwell_k4_J_2_mu_m2**2.)/sim
dwell_sd_k4_J_2_mu_m2 = np.sqrt(dwell_avg2_k4_J_2_mu_m2-dwell_avg_k4_J_2_mu_m2**2) 

dwell_avg_k5_J_2_mu_m2 = sum(dwell_k5_J_2_mu_m2)/sim
dwell_avg2_k5_J_2_mu_m2 = sum(dwell_k5_J_2_mu_m2**2.)/sim
dwell_sd_k5_J_2_mu_m2 = np.sqrt(dwell_avg2_k5_J_2_mu_m2-dwell_avg_k5_J_2_mu_m2**2)

dwell_avg_k6_J_2_mu_m2 = sum(dwell_k6_J_2_mu_m2)/sim
dwell_avg2_k6_J_2_mu_m2 = sum(dwell_k6_J_2_mu_m2**2.)/sim
dwell_sd_k6_J_2_mu_m2 = np.sqrt(dwell_avg2_k6_J_2_mu_m2-dwell_avg_k6_J_2_mu_m2**2)

dwell_avg_k7_J_2_mu_m2 = sum(dwell_k7_J_2_mu_m2)/sim
dwell_avg2_k7_J_2_mu_m2 = sum(dwell_k7_J_2_mu_m2**2.)/sim
dwell_sd_k7_J_2_mu_m2 = np.sqrt(dwell_avg2_k7_J_2_mu_m2-dwell_avg_k7_J_2_mu_m2**2)

dwell_avg_k8_J_2_mu_m2 = sum(dwell_k8_J_2_mu_m2)/sim
dwell_avg2_k8_J_2_mu_m2 = sum(dwell_k8_J_2_mu_m2**2.)/sim
dwell_sd_k8_J_2_mu_m2 = np.sqrt(dwell_avg2_k8_J_2_mu_m2-dwell_avg_k8_J_2_mu_m2**2)    

dwell_avg_k9_J_2_mu_m2 = sum(dwell_k9_J_2_mu_m2)/sim
dwell_avg2_k9_J_2_mu_m2 = sum(dwell_k9_J_2_mu_m2**2.)/sim
dwell_sd_k9_J_2_mu_m2 = np.sqrt(dwell_avg2_k9_J_2_mu_m2-dwell_avg_k9_J_2_mu_m2**2)

dwell_avg_k10_J_2_mu_m2 = sum(dwell_k10_J_2_mu_m2)/sim
dwell_avg2_k10_J_2_mu_m2 = sum(dwell_k10_J_2_mu_m2**2.)/sim
dwell_sd_k10_J_2_mu_m2 = np.sqrt(dwell_avg2_k10_J_2_mu_m2-dwell_avg_k10_J_2_mu_m2**2)




#%%

dwell_avg_k0_J_3_mu_m3 = sum(dwell_k0_J_3_mu_m3)/sim
dwell_avg2_k0_J_3_mu_m3 = sum(dwell_k0_J_3_mu_m3**2.)/sim
dwell_sd_k0_J_3_mu_m3 = np.sqrt(dwell_avg2_k0_J_3_mu_m3-dwell_avg_k0_J_3_mu_m3**2)

dwell_avg_k1_J_3_mu_m3 = sum(dwell_k1_J_3_mu_m3)/sim
dwell_avg2_k1_J_3_mu_m3 = sum(dwell_k1_J_3_mu_m3**2.)/sim
dwell_sd_k1_J_3_mu_m3 = np.sqrt(dwell_avg2_k1_J_3_mu_m3-dwell_avg_k1_J_3_mu_m3**2)

dwell_avg_k2_J_3_mu_m3 = sum(dwell_k2_J_3_mu_m3)/sim
dwell_avg2_k2_J_3_mu_m3 = sum(dwell_k2_J_3_mu_m3**2.)/sim
dwell_sd_k2_J_3_mu_m3 = np.sqrt(dwell_avg2_k2_J_3_mu_m3-dwell_avg_k2_J_3_mu_m3**2)   

dwell_avg_k3_J_3_mu_m3 = sum(dwell_k3_J_3_mu_m3)/sim
dwell_avg2_k3_J_3_mu_m3 = sum(dwell_k3_J_3_mu_m3**2.)/sim
dwell_sd_k3_J_3_mu_m3 = np.sqrt(dwell_avg2_k3_J_3_mu_m3-dwell_avg_k3_J_3_mu_m3**2)    

dwell_avg_k4_J_3_mu_m3 = sum(dwell_k4_J_3_mu_m3)/sim
dwell_avg2_k4_J_3_mu_m3 = sum(dwell_k4_J_3_mu_m3**2.)/sim
dwell_sd_k4_J_3_mu_m3 = np.sqrt(dwell_avg2_k4_J_3_mu_m3-dwell_avg_k4_J_3_mu_m3**2) 

dwell_avg_k5_J_3_mu_m3 = sum(dwell_k5_J_3_mu_m3)/sim
dwell_avg2_k5_J_3_mu_m3 = sum(dwell_k5_J_3_mu_m3**2.)/sim
dwell_sd_k5_J_3_mu_m3 = np.sqrt(dwell_avg2_k5_J_3_mu_m3-dwell_avg_k5_J_3_mu_m3**2)

dwell_avg_k6_J_3_mu_m3 = sum(dwell_k6_J_3_mu_m3)/sim
dwell_avg2_k6_J_3_mu_m3 = sum(dwell_k6_J_3_mu_m3**2.)/sim
dwell_sd_k6_J_3_mu_m3 = np.sqrt(dwell_avg2_k6_J_3_mu_m3-dwell_avg_k6_J_3_mu_m3**2)

dwell_avg_k7_J_3_mu_m3 = sum(dwell_k7_J_3_mu_m3)/sim
dwell_avg2_k7_J_3_mu_m3 = sum(dwell_k7_J_3_mu_m3**2.)/sim
dwell_sd_k7_J_3_mu_m3 = np.sqrt(dwell_avg2_k7_J_3_mu_m3-dwell_avg_k7_J_3_mu_m3**2)

dwell_avg_k8_J_3_mu_m3 = sum(dwell_k8_J_3_mu_m3)/sim
dwell_avg2_k8_J_3_mu_m3 = sum(dwell_k8_J_3_mu_m3**2.)/sim
dwell_sd_k8_J_3_mu_m3 = np.sqrt(dwell_avg2_k8_J_3_mu_m3-dwell_avg_k8_J_3_mu_m3**2)    

dwell_avg_k9_J_3_mu_m3 = sum(dwell_k9_J_3_mu_m3)/sim
dwell_avg2_k9_J_3_mu_m3 = sum(dwell_k9_J_3_mu_m3**2.)/sim
dwell_sd_k9_J_3_mu_m3 = np.sqrt(dwell_avg2_k9_J_3_mu_m3-dwell_avg_k9_J_3_mu_m3**2)

dwell_avg_k10_J_3_mu_m3 = sum(dwell_k10_J_3_mu_m3)/sim
dwell_avg2_k10_J_3_mu_m3 = sum(dwell_k10_J_3_mu_m3**2.)/sim
dwell_sd_k10_J_3_mu_m3 = np.sqrt(dwell_avg2_k10_J_3_mu_m3-dwell_avg_k10_J_3_mu_m3**2)

#%%

dwell_avg_k0_J_4_mu_4 = sum(dwell_k0_J_4_mu_4)/sim
dwell_avg2_k0_J_4_mu_4 = sum(dwell_k0_J_4_mu_4**2.)/sim
dwell_sd_k0_J_4_mu_4 = np.sqrt(dwell_avg2_k0_J_4_mu_4-dwell_avg_k0_J_4_mu_4**2)

dwell_avg_k1_J_4_mu_4 = sum(dwell_k1_J_4_mu_4)/sim
dwell_avg2_k1_J_4_mu_4 = sum(dwell_k1_J_4_mu_4**2.)/sim
dwell_sd_k1_J_4_mu_4 = np.sqrt(dwell_avg2_k1_J_4_mu_4-dwell_avg_k1_J_4_mu_4**2)

dwell_avg_k2_J_4_mu_4 = sum(dwell_k2_J_4_mu_4)/sim
dwell_avg2_k2_J_4_mu_4 = sum(dwell_k2_J_4_mu_4**2.)/sim
dwell_sd_k2_J_4_mu_4 = np.sqrt(dwell_avg2_k2_J_4_mu_4-dwell_avg_k2_J_4_mu_4**2)   

dwell_avg_k3_J_4_mu_4 = sum(dwell_k3_J_4_mu_4)/sim
dwell_avg2_k3_J_4_mu_4 = sum(dwell_k3_J_4_mu_4**2.)/sim
dwell_sd_k3_J_4_mu_4 = np.sqrt(dwell_avg2_k3_J_4_mu_4-dwell_avg_k3_J_4_mu_4**2)    

dwell_avg_k4_J_4_mu_4 = sum(dwell_k4_J_4_mu_4)/sim
dwell_avg2_k4_J_4_mu_4 = sum(dwell_k4_J_4_mu_4**2.)/sim
dwell_sd_k4_J_4_mu_4 = np.sqrt(dwell_avg2_k4_J_4_mu_4-dwell_avg_k4_J_4_mu_4**2) 

dwell_avg_k5_J_4_mu_4 = sum(dwell_k5_J_4_mu_4)/sim
dwell_avg2_k5_J_4_mu_4 = sum(dwell_k5_J_4_mu_4**2.)/sim
dwell_sd_k5_J_4_mu_4 = np.sqrt(dwell_avg2_k5_J_4_mu_4-dwell_avg_k5_J_4_mu_4**2)

dwell_avg_k6_J_4_mu_4 = sum(dwell_k6_J_4_mu_4)/sim
dwell_avg2_k6_J_4_mu_4 = sum(dwell_k6_J_4_mu_4**2.)/sim
dwell_sd_k6_J_4_mu_4 = np.sqrt(dwell_avg2_k6_J_4_mu_4-dwell_avg_k6_J_4_mu_4**2) 

dwell_avg_k7_J_4_mu_4 = sum(dwell_k7_J_4_mu_4)/sim
dwell_avg2_k7_J_4_mu_4 = sum(dwell_k7_J_4_mu_4**2.)/sim
dwell_sd_k7_J_4_mu_4 = np.sqrt(dwell_avg2_k7_J_4_mu_4-dwell_avg_k7_J_4_mu_4**2)

dwell_avg_k8_J_4_mu_4 = sum(dwell_k8_J_4_mu_4)/sim
dwell_avg2_k8_J_4_mu_4 = sum(dwell_k8_J_4_mu_4**2.)/sim
dwell_sd_k8_J_4_mu_4 = np.sqrt(dwell_avg2_k8_J_4_mu_4-dwell_avg_k8_J_4_mu_4**2)  

dwell_avg_k9_J_4_mu_4 = sum(dwell_k9_J_4_mu_4)/sim
dwell_avg2_k9_J_4_mu_4 = sum(dwell_k9_J_4_mu_4**2.)/sim
dwell_sd_k9_J_4_mu_4 = np.sqrt(dwell_avg2_k9_J_4_mu_4-dwell_avg_k9_J_4_mu_4**2)

dwell_avg_k10_J_4_mu_4 = sum(dwell_k10_J_4_mu_4)/sim
dwell_avg2_k10_J_4_mu_4 = sum(dwell_k10_J_4_mu_4**2.)/sim
dwell_sd_k10_J_4_mu_4 = np.sqrt(dwell_avg2_k10_J_4_mu_4-dwell_avg_k10_J_4_mu_4**2) 

#%%

dwell_avg_k0_J_5_mu_5 = sum(dwell_k0_J_5_mu_5)/sim
dwell_avg2_k0_J_5_mu_5 = sum(dwell_k0_J_5_mu_5**2.)/sim
dwell_sd_k0_J_5_mu_5 = np.sqrt(dwell_avg2_k0_J_5_mu_5-dwell_avg_k0_J_5_mu_5**2)

dwell_avg_k1_J_5_mu_5 = sum(dwell_k1_J_5_mu_5)/sim
dwell_avg2_k1_J_5_mu_5 = sum(dwell_k1_J_5_mu_5**2.)/sim
dwell_sd_k1_J_5_mu_5 = np.sqrt(dwell_avg2_k1_J_5_mu_5-dwell_avg_k1_J_5_mu_5**2)

dwell_avg_k2_J_5_mu_5 = sum(dwell_k2_J_5_mu_5)/sim
dwell_avg2_k2_J_5_mu_5 = sum(dwell_k2_J_5_mu_5**2.)/sim
dwell_sd_k2_J_5_mu_5 = np.sqrt(dwell_avg2_k2_J_5_mu_5-dwell_avg_k2_J_5_mu_5**2)   

dwell_avg_k3_J_5_mu_5 = sum(dwell_k3_J_5_mu_5)/sim
dwell_avg2_k3_J_5_mu_5 = sum(dwell_k3_J_5_mu_5**2.)/sim
dwell_sd_k3_J_5_mu_5 = np.sqrt(dwell_avg2_k3_J_5_mu_5-dwell_avg_k3_J_5_mu_5**2)    

dwell_avg_k4_J_5_mu_5 = sum(dwell_k4_J_5_mu_5)/sim
dwell_avg2_k4_J_5_mu_5 = sum(dwell_k4_J_5_mu_5**2.)/sim
dwell_sd_k4_J_5_mu_5 = np.sqrt(dwell_avg2_k4_J_5_mu_5-dwell_avg_k4_J_5_mu_5**2) 

dwell_avg_k5_J_5_mu_5 = sum(dwell_k5_J_5_mu_5)/sim
dwell_avg2_k5_J_5_mu_5 = sum(dwell_k5_J_5_mu_5**2.)/sim
dwell_sd_k5_J_5_mu_5 = np.sqrt(dwell_avg2_k5_J_5_mu_5-dwell_avg_k5_J_5_mu_5**2)

dwell_avg_k6_J_5_mu_5 = sum(dwell_k6_J_5_mu_5)/sim
dwell_avg2_k6_J_5_mu_5 = sum(dwell_k6_J_5_mu_5**2.)/sim
dwell_sd_k6_J_5_mu_5 = np.sqrt(dwell_avg2_k6_J_5_mu_5-dwell_avg_k6_J_5_mu_5**2) 

dwell_avg_k7_J_5_mu_5 = sum(dwell_k7_J_5_mu_5)/sim
dwell_avg2_k7_J_5_mu_5 = sum(dwell_k7_J_5_mu_5**2.)/sim
dwell_sd_k7_J_5_mu_5 = np.sqrt(dwell_avg2_k7_J_5_mu_5-dwell_avg_k7_J_5_mu_5**2)

dwell_avg_k8_J_5_mu_5 = sum(dwell_k8_J_5_mu_5)/sim
dwell_avg2_k8_J_5_mu_5 = sum(dwell_k8_J_5_mu_5**2.)/sim
dwell_sd_k8_J_5_mu_5 = np.sqrt(dwell_avg2_k8_J_5_mu_5-dwell_avg_k8_J_5_mu_5**2)  

dwell_avg_k9_J_5_mu_5 = sum(dwell_k9_J_5_mu_5)/sim
dwell_avg2_k9_J_5_mu_5 = sum(dwell_k9_J_5_mu_5**2.)/sim
dwell_sd_k9_J_5_mu_5 = np.sqrt(dwell_avg2_k9_J_5_mu_5-dwell_avg_k9_J_5_mu_5**2)

dwell_avg_k10_J_5_mu_5 = sum(dwell_k10_J_5_mu_5)/sim
dwell_avg2_k10_J_5_mu_5 = sum(dwell_k10_J_5_mu_5**2.)/sim
dwell_sd_k10_J_5_mu_5 = np.sqrt(dwell_avg2_k10_J_5_mu_5-dwell_avg_k10_J_5_mu_5**2) 




#%%

dwell_avg_J_0_mu_0 = np.array([dwell_avg_k0_J_0_mu_0,dwell_avg_k1_J_0_mu_0,dwell_avg_k2_J_0_mu_0,dwell_avg_k3_J_0_mu_0,
                              dwell_avg_k4_J_0_mu_0,dwell_avg_k5_J_0_mu_0,dwell_avg_k6_J_0_mu_0,dwell_avg_k7_J_0_mu_0,
                              dwell_avg_k8_J_0_mu_0,dwell_avg_k9_J_0_mu_0,dwell_avg_k10_J_0_mu_0])

dwell_sd_J_0_mu_0 = np.array([dwell_sd_k0_J_0_mu_0,dwell_sd_k1_J_0_mu_0,dwell_sd_k2_J_0_mu_0,dwell_sd_k3_J_0_mu_0,
                             dwell_sd_k4_J_0_mu_0,dwell_sd_k5_J_0_mu_0,dwell_sd_k6_J_0_mu_0,dwell_sd_k7_J_0_mu_0,
                             dwell_sd_k8_J_0_mu_0,dwell_sd_k9_J_0_mu_0,dwell_sd_k10_J_0_mu_0])


dwell_avg_J_1_mu_m1 = np.array([dwell_avg_k0_J_1_mu_m1,dwell_avg_k1_J_1_mu_m1,dwell_avg_k2_J_1_mu_m1,dwell_avg_k3_J_1_mu_m1,
                              dwell_avg_k4_J_1_mu_m1,dwell_avg_k5_J_1_mu_m1,dwell_avg_k6_J_1_mu_m1,dwell_avg_k7_J_1_mu_m1,
                              dwell_avg_k8_J_1_mu_m1,dwell_avg_k9_J_1_mu_m1,dwell_avg_k10_J_1_mu_m1])

dwell_sd_J_1_mu_m1 = np.array([dwell_sd_k0_J_1_mu_m1,dwell_sd_k1_J_1_mu_m1,dwell_sd_k2_J_1_mu_m1,dwell_sd_k3_J_1_mu_m1,
                             dwell_sd_k4_J_1_mu_m1,dwell_sd_k5_J_1_mu_m1,dwell_sd_k6_J_1_mu_m1,dwell_sd_k7_J_1_mu_m1,
                             dwell_sd_k8_J_1_mu_m1,dwell_sd_k9_J_1_mu_m1,dwell_sd_k10_J_1_mu_m1])


dwell_avg_J_2_mu_m2 = np.array([dwell_avg_k0_J_2_mu_m2,dwell_avg_k1_J_2_mu_m2,dwell_avg_k2_J_2_mu_m2,dwell_avg_k3_J_2_mu_m2,
                              dwell_avg_k4_J_2_mu_m2,dwell_avg_k5_J_2_mu_m2,dwell_avg_k6_J_2_mu_m2,dwell_avg_k7_J_2_mu_m2,
                              dwell_avg_k8_J_2_mu_m2,dwell_avg_k9_J_2_mu_m2,dwell_avg_k10_J_2_mu_m2])

dwell_sd_J_2_mu_m2 = np.array([dwell_sd_k0_J_2_mu_m2,dwell_sd_k1_J_2_mu_m2,dwell_sd_k2_J_2_mu_m2,dwell_sd_k3_J_2_mu_m2,
                             dwell_sd_k4_J_2_mu_m2,dwell_sd_k5_J_2_mu_m2,dwell_sd_k6_J_2_mu_m2,dwell_sd_k7_J_2_mu_m2,
                             dwell_sd_k8_J_2_mu_m2,dwell_sd_k9_J_2_mu_m2,dwell_sd_k10_J_2_mu_m2])


dwell_avg_J_3_mu_m3 = np.array([dwell_avg_k0_J_3_mu_m3,dwell_avg_k1_J_3_mu_m3,dwell_avg_k2_J_3_mu_m3,dwell_avg_k3_J_3_mu_m3,
                              dwell_avg_k4_J_3_mu_m3,dwell_avg_k5_J_3_mu_m3,dwell_avg_k6_J_3_mu_m3,dwell_avg_k7_J_3_mu_m3,
                              dwell_avg_k8_J_3_mu_m3,dwell_avg_k9_J_3_mu_m3,dwell_avg_k10_J_3_mu_m3])

dwell_sd_J_3_mu_m3 = np.array([dwell_sd_k0_J_3_mu_m3,dwell_sd_k1_J_3_mu_m3,dwell_sd_k2_J_3_mu_m3,dwell_sd_k3_J_3_mu_m3,
                             dwell_sd_k4_J_3_mu_m3,dwell_sd_k5_J_3_mu_m3,dwell_sd_k6_J_3_mu_m3,dwell_sd_k7_J_3_mu_m3,
                             dwell_sd_k8_J_3_mu_m3,dwell_sd_k9_J_3_mu_m3,dwell_sd_k10_J_3_mu_m3])


dwell_avg_J_4_mu_4 = np.array([dwell_avg_k0_J_4_mu_4,dwell_avg_k1_J_4_mu_4,dwell_avg_k2_J_4_mu_4,dwell_avg_k3_J_4_mu_4,
                              dwell_avg_k4_J_4_mu_4,dwell_avg_k5_J_4_mu_4,dwell_avg_k6_J_4_mu_4,dwell_avg_k7_J_4_mu_4,
                              dwell_avg_k8_J_4_mu_4,dwell_avg_k9_J_4_mu_4,dwell_avg_k10_J_4_mu_4])

dwell_sd_J_4_mu_4 = np.array([dwell_sd_k0_J_4_mu_4,dwell_sd_k1_J_4_mu_4,dwell_sd_k2_J_4_mu_4,dwell_sd_k3_J_4_mu_4,
                             dwell_sd_k4_J_4_mu_4,dwell_sd_k5_J_4_mu_4,dwell_sd_k6_J_4_mu_4,dwell_sd_k7_J_4_mu_4,
                             dwell_sd_k8_J_4_mu_4,dwell_sd_k9_J_4_mu_4,dwell_sd_k10_J_4_mu_4])

dwell_avg_J_5_mu_5 = np.array([dwell_avg_k0_J_5_mu_5,dwell_avg_k1_J_5_mu_5,dwell_avg_k2_J_5_mu_5,dwell_avg_k3_J_5_mu_5,
                              dwell_avg_k4_J_5_mu_5,dwell_avg_k5_J_5_mu_5,dwell_avg_k6_J_5_mu_5,dwell_avg_k7_J_5_mu_5,
                              dwell_avg_k8_J_5_mu_5,dwell_avg_k9_J_5_mu_5,dwell_avg_k10_J_5_mu_5])

dwell_sd_J_5_mu_5 = np.array([dwell_sd_k0_J_5_mu_5,dwell_sd_k1_J_5_mu_5,dwell_sd_k2_J_5_mu_5,dwell_sd_k3_J_5_mu_5,
                             dwell_sd_k4_J_5_mu_5,dwell_sd_k5_J_5_mu_5,dwell_sd_k6_J_5_mu_5,dwell_sd_k7_J_5_mu_5,
                             dwell_sd_k8_J_5_mu_5,dwell_sd_k9_J_5_mu_5,dwell_sd_k10_J_5_mu_5])


'''
#%%

plt.rcParams["figure.figsize"] = [8.0,6.0]
plt.grid()
plt.title('$\mu = 0.00$',size=17)
plt.xlabel('k',size=14)
plt.ylabel('<n(k)>',size=14)
plt.errorbar(k,dwell_avg_J_0_mu_0,dwell_sd_J_0_mu_0,linestyle='',marker='o',capsize=3,label='J=0')
plt.errorbar(k,dwell_avg_J_1_mu_0,dwell_sd_J_1_mu_0,linestyle='',marker='o',capsize=3,label='J=1')
plt.errorbar(k,dwell_avg_J_2_mu_0,dwell_sd_J_2_mu_0,linestyle='',marker='o',capsize=3,label ='J=2')
plt.errorbar(k,dwell_avg_J_3_mu_0,dwell_sd_J_3_mu_0,linestyle='',marker='o',capsize=3,label ='J=3')

plt.legend()

'''
#%%   

plt.rcParams["figure.figsize"] = [8.0,6.0]
plt.grid()
plt.title('$<\phi>=2/3$',size=17)
plt.xlabel('k',size=14)
plt.ylabel('<n(k)>',size=14)
plt.yticks(np.arange(0, 310, 20))
plt.errorbar(k,dwell_avg_J_0_mu_0,dwell_sd_J_0_mu_0,linestyle='',marker='o',capsize=3,label='J=0, $\mu=0.69$')
plt.errorbar(k,dwell_avg_J_1_mu_m1,dwell_sd_J_1_mu_m1,linestyle='',marker='o',capsize=3,label='J=1, $\mu=-0.57$')
plt.errorbar(k,dwell_avg_J_2_mu_m2,dwell_sd_J_2_mu_m2,linestyle='',marker='o',capsize=3,label='J=2, $\mu=-2.26$')#
plt.errorbar(k,dwell_avg_J_3_mu_m3,dwell_sd_J_3_mu_m3,linestyle='',marker='o',capsize=3,label='J=3, $\mu=-1.74$')
plt.errorbar(k,dwell_avg_J_4_mu_4,dwell_sd_J_4_mu_4,linestyle='',marker='o',capsize=3,label='J=4, $\mu=-3.85$')
plt.errorbar(k,dwell_avg_J_5_mu_5,dwell_sd_J_5_mu_5,linestyle='',marker='o',capsize=3,label='J=5, $\mu=-4.79$')

plt.legend(loc='upper right')

#%% Ratio between k=0 and k=1

J=np.arange(0,6,1)

R_J_0 = dwell_avg_k0_J_0_mu_0/dwell_avg_k1_J_0_mu_0
R_J_1 = dwell_avg_k0_J_1_mu_m1/dwell_avg_k1_J_1_mu_m1
R_J_2 = dwell_avg_k0_J_2_mu_m2/dwell_avg_k1_J_2_mu_m2
R_J_3 = dwell_avg_k0_J_3_mu_m3/dwell_avg_k1_J_3_mu_m3
R_J_4 = dwell_avg_k0_J_4_mu_4/dwell_avg_k1_J_4_mu_4
R_J_5 = dwell_avg_k0_J_5_mu_5/dwell_avg_k1_J_5_mu_5

dR_J_0 = R_J_0*(dwell_sd_k0_J_0_mu_0/dwell_avg_k0_J_0_mu_0 + dwell_sd_k1_J_0_mu_0/dwell_avg_k1_J_0_mu_0)
dR_J_1 = R_J_1*(dwell_sd_k0_J_1_mu_m1/dwell_avg_k0_J_1_mu_m1 + dwell_sd_k1_J_1_mu_m1/dwell_avg_k1_J_1_mu_m1)
dR_J_2 = R_J_2*(dwell_sd_k0_J_2_mu_m2/dwell_avg_k0_J_2_mu_m2 + dwell_sd_k1_J_2_mu_m2/dwell_avg_k1_J_2_mu_m2)
dR_J_3 = R_J_3*(dwell_sd_k0_J_3_mu_m3/dwell_avg_k0_J_3_mu_m3 + dwell_sd_k1_J_3_mu_m3/dwell_avg_k1_J_3_mu_m3)
dR_J_4 = R_J_4*(dwell_sd_k0_J_4_mu_4/dwell_avg_k0_J_4_mu_4 + dwell_sd_k1_J_4_mu_4/dwell_avg_k1_J_4_mu_4)
dR_J_5 = R_J_5*(dwell_sd_k0_J_5_mu_5/dwell_avg_k0_J_5_mu_5 + dwell_sd_k1_J_5_mu_5/dwell_avg_k1_J_5_mu_5)



#%%

R = np.array([R_J_0,R_J_1,R_J_2,R_J_3,R_J_4,R_J_5])
dR = np.array([dR_J_0,dR_J_1,dR_J_2,dR_J_3,dR_J_4,dR_J_5])

np.savetxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Dwell_times/Ratio_phi=0.66',(np.c_[J,R,dR]))

#%%

R_05 = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Dwell_times/Ratio_phi=0.5')
R_03 = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Dwell_times/Ratio_phi=0.33')
R_06 = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Dwell_times/Ratio_phi=0.66')

plt.rcParams["figure.figsize"] = [8.0,6.0]
#plt.yticks(np.arange(0, 27, 5))
plt.grid()
plt.xlabel('J',size=14)
plt.ylabel('n(0)/n(1)',size=14)
plt.errorbar(J,R_05[:,1],R_05[:,2],marker='o',linestyle='',capsize=3,label='<$\phi$>=1/2')
plt.errorbar(J,R_03[:,1],R_03[:,2],marker='o',linestyle='',capsize=3,label='<$\phi$>=1/3')
plt.errorbar(J,R_06[:,1],R_06[:,2],marker='o',linestyle='',capsize=3,label='<$\phi$>=2/3')


plt.legend()

plt.show()

#%% 

def f(t,a,b):
    return a*np.exp(-b*t)

#%%

plt.rcParams["figure.figsize"] = [8.0,6.0]
plt.grid()
plt.title('$J=1$, $\mu=-1$',size=17)
plt.xlabel('k',size=14)
plt.ylabel('<n(k)>/<n(0)>',size=14)

plt.errorbar(k,dwell_avg_J_1_mu_m1/dwell_avg_k0_J_1_mu_m1,dwell_sd_J_1_mu_m1/dwell_sd_k0_J_1_mu_m1,linestyle='',marker='o',capsize=3,label='J=1')
#plt.legend()

#%%

plt.plot(k,dwell_sd_J_0_mu_0,linestyle='',marker='o',label='J=0, $\mu=0.00$')
plt.plot(k,dwell_sd_J_1_mu_m1,linestyle='',marker='o',label='J=1, $\mu=-1.00$')
plt.plot(k,dwell_sd_J_2_mu_m2,linestyle='',marker='o',label='J=2, $\mu=-2.00$')
plt.plot(k,dwell_sd_J_3_mu_m3,linestyle='',marker='o',label='J=3, $\mu=-3.00$')

plt.legend()

#%%

u_dwell_J_1_mu_m1, co_dwell_J_1_mu_m1 = np.unique(dwell_k4_J_1_mu_m1,return_counts=True)

cn = np.linspace(0,np.max(u_dwell_J_1_mu_m1),1000)

po_dwell_J_1_mu_m1, pc_dwell_J_1_mu_m1 = curve_fit(f,u_dwell_J_1_mu_m1-1,co_dwell_J_1_mu_m1/sim)
model_dwell_J_1_mu_m1 = f(cn,*po_dwell_J_1_mu_m1)


plt.title('P(k=4), $J=1$, $\mu = -1$',size=15)
plt.ylabel('Frecuency')
plt.xlabel('n')
plt.bar(u_dwell_J_1_mu_m1-1,co_dwell_J_1_mu_m1/sim,color='blue',alpha=0.5)
plt.plot(cn,model_dwell_J_1_mu_m1,color='red',label='$ae^{-bn}$')
plt.plot([],[],linestyle='',label='a=%.2f'%(po_dwell_J_1_mu_m1[0]))
plt.plot([],[],linestyle='',label='b=%.2f'%(po_dwell_J_1_mu_m1[1]))

plt.legend(prop={'size': 15})