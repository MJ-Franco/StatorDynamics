#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 11:42:24 2021

Program to analyse the dwell times of resurrection

@author: mariajose
"""

import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
#from scipy.optimize import curve_fit

#%%

files1=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/GTraces_Res_J=0.00  .dat'] #All resurrection simulations
data1=[]
files11=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/GTraces_Res_J=0.25  .dat'] #All resurrection simulations
data11=[]
files12=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/GTraces_Res_J=0.50  .dat'] #All resurrection simulations
data12=[]
files13=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/GTraces_Res_J=0.75  .dat'] #All resurrection simulations
data13=[]

files2=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/GTraces_Res_J=1.00  .dat'] #Average of resurrection simulations
data2=[]
files21=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/GTraces_Res_J=1.25  .dat'] #Average of resurrection simulations
data21=[]
files22=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/GTraces_Res_J=1.50  .dat'] #Average of resurrection simulations
data22=[]
files23=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/GTraces_Res_J=1.75  .dat'] #Average of resurrection simulations
data23=[]

files3=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/GTraces_Res_J=2.00  .dat'] #Average of resurrection simulations
data3=[]
files31=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/GTraces_Res_J=2.25  .dat'] #Average of resurrection simulations
data31=[]
files32=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/GTraces_Res_J=2.50  .dat'] #Average of resurrection simulations
data32=[]
files33=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/GTraces_Res_J=2.75  .dat'] #Average of resurrection simulations
data33=[]

files4=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/GTraces_Res_J=3.00  .dat'] #Average of resurrection simulations
data4=[]
files41=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/GTraces_Res_J=3.25  .dat'] #Average of resurrection simulations
data41=[]
files42=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/GTraces_Res_J=3.50  .dat'] #Average of resurrection simulations
data42=[]
files43=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/GTraces_Res_J=3.75  .dat'] #Average of resurrection simulations
data43=[]

files5=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/GTraces_Res_J=4.00  .dat'] #Average of resurrection simulations
data5=[]
files51=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/GTraces_Res_J=4.25  .dat'] #Average of resurrection simulations
data51=[]
files52=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/GTraces_Res_J=4.50  .dat'] #Average of resurrection simulations
data52=[]
files53=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/GTraces_Res_J=4.75  .dat'] #Average of resurrection simulations
data53=[]

files6=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/GTraces_Res_J=5.00  .dat'] #Average of resurrection simulations
data6=[]

for data_file in files1:
    data1.append(np.loadtxt(data_file))
    
for data_file in files11:
    data11.append(np.loadtxt(data_file))

for data_file in files12:
    data12.append(np.loadtxt(data_file))

for data_file in files13:
    data13.append(np.loadtxt(data_file))


for data_file in files2:
    data2.append(np.loadtxt(data_file))
    
for data_file in files21:
    data21.append(np.loadtxt(data_file))
    
for data_file in files22:
    data22.append(np.loadtxt(data_file))

for data_file in files23:
    data23.append(np.loadtxt(data_file))
    
    
for data_file in files3:
    data3.append(np.loadtxt(data_file))
    
for data_file in files31:
    data31.append(np.loadtxt(data_file))
    
for data_file in files32:
    data32.append(np.loadtxt(data_file))

for data_file in files33:
    data33.append(np.loadtxt(data_file))

    
for data_file in files4:
    data4.append(np.loadtxt(data_file))
    
for data_file in files41:
    data41.append(np.loadtxt(data_file))
    
for data_file in files42:
    data42.append(np.loadtxt(data_file))
    
for data_file in files43:
    data43.append(np.loadtxt(data_file))
    
    
for data_file in files5:
    data5.append(np.loadtxt(data_file))
    
for data_file in files51:
    data51.append(np.loadtxt(data_file))

for data_file in files52:
    data52.append(np.loadtxt(data_file))
    
for data_file in files53:
    data53.append(np.loadtxt(data_file))    
    
    
    
for data_file in files6:
    data6.append(np.loadtxt(data_file))
    
    
    
All_Res_J_0=data1[0] #All traces
All_Res_J_025=data11[0]
All_Res_J_050=data12[0]
All_Res_J_075=data13[0]

All_Res_J_1=data2[0]
All_Res_J_125=data21[0]
All_Res_J_150=data22[0]
All_Res_J_175=data23[0]

All_Res_J_2=data3[0]
All_Res_J_225=data31[0]
All_Res_J_250=data32[0]
All_Res_J_275=data33[0]

All_Res_J_3=data4[0]
All_Res_J_325=data41[0]
All_Res_J_350=data42[0]
All_Res_J_375=data43[0]

All_Res_J_4=data5[0]
All_Res_J_425=data51[0]
All_Res_J_450=data52[0]
All_Res_J_475=data53[0]

All_Res_J_5=data6[0]

#%%

sim=10000

S_0=28 #Total number of steps of the simulation
S_025=33
S_050=38
S_075=45

S_1=55
S_125=67
S_150=85
S_175=100

S_2=125
S_225=160
S_250=200
S_275=260

S_3=325
S_325=425
S_350=550
S_375=700

S_4=900
S_425=1100
S_450=1400
S_475=1850

S_5=2300

nmax=13

#%%

dwell_J_0_k0 = np.zeros(sim)
dwell_J_025_k0 = np.zeros(sim)
dwell_J_050_k0 = np.zeros(sim)
dwell_J_075_k0 = np.zeros(sim)

dwell_J_1_k0 = np.zeros(sim)
dwell_J_125_k0 = np.zeros(sim)
dwell_J_150_k0 = np.zeros(sim)
dwell_J_175_k0 = np.zeros(sim)

dwell_J_2_k0 = np.zeros(sim)
dwell_J_225_k0 = np.zeros(sim)
dwell_J_250_k0 = np.zeros(sim)
dwell_J_275_k0 = np.zeros(sim)

dwell_J_3_k0 = np.zeros(sim)
dwell_J_325_k0 = np.zeros(sim)
dwell_J_350_k0 = np.zeros(sim)
dwell_J_375_k0 = np.zeros(sim)

dwell_J_4_k0 = np.zeros(sim)
dwell_J_425_k0 = np.zeros(sim)
dwell_J_450_k0 = np.zeros(sim)
dwell_J_475_k0 = np.zeros(sim)

dwell_J_5_k0 = np.zeros(sim)


dwell_J_0_k1 = np.zeros(sim)
dwell_J_025_k1 = np.zeros(sim)
dwell_J_050_k1 = np.zeros(sim)
dwell_J_075_k1 = np.zeros(sim)

dwell_J_1_k1 = np.zeros(sim)
dwell_J_125_k1 = np.zeros(sim)
dwell_J_150_k1 = np.zeros(sim)
dwell_J_175_k1 = np.zeros(sim)

dwell_J_2_k1 = np.zeros(sim)
dwell_J_225_k1 = np.zeros(sim)
dwell_J_250_k1 = np.zeros(sim)
dwell_J_275_k1 = np.zeros(sim)

dwell_J_3_k1 = np.zeros(sim)
dwell_J_325_k1 = np.zeros(sim)
dwell_J_350_k1 = np.zeros(sim)
dwell_J_375_k1 = np.zeros(sim)

dwell_J_4_k1 = np.zeros(sim)
dwell_J_425_k1 = np.zeros(sim)
dwell_J_450_k1 = np.zeros(sim)
dwell_J_475_k1 = np.zeros(sim)

dwell_J_5_k1 = np.zeros(sim)

#%% Dwell time J=0s for k=0,k=1

for i in range(sim):
    for j in range(S_0*i,S_0*(i+1)-1):
        if All_Res_J_0[j,2]*nmax==1.:
            dwell_J_0_k0[i] = j-S_0*i
            break

for i in range(sim):     
    for j in range(int(S_0*i+dwell_J_0_k0[i]),S_0*(i+1)-1):
        if All_Res_J_0[j,2]*nmax!=1.:
            
            dwell_J_0_k1[i] = j-S_0*i
            break
        
for i in range(sim):
    for j in range(S_025*i,S_025*(i+1)-1):
        if All_Res_J_025[j,2]*nmax==1.:
            dwell_J_025_k0[i] = j-S_025*i
            break

for i in range(sim):     
    for j in range(int(S_025*i+dwell_J_025_k0[i]),S_025*(i+1)-1):
        if All_Res_J_025[j,2]*nmax!=1.:
            
            dwell_J_025_k1[i] = j-S_025*i
            break
        
for i in range(sim):
    for j in range(S_050*i,S_050*(i+1)-1):
        if All_Res_J_050[j,2]*nmax==1.:
            dwell_J_050_k0[i] = j-S_050*i
            break

for i in range(sim):     
    for j in range(int(S_050*i+dwell_J_050_k0[i]),S_050*(i+1)-1):
        if All_Res_J_050[j,2]*nmax!=1.:
            
            dwell_J_050_k1[i] = j-S_050*i
            break
        
for i in range(sim):
    for j in range(S_075*i,S_075*(i+1)-1):
        if All_Res_J_075[j,2]*nmax==1.:
            dwell_J_075_k0[i] = j-S_075*i
            break

for i in range(sim):     
    for j in range(int(S_075*i+dwell_J_075_k0[i]),S_075*(i+1)-1):
        if All_Res_J_075[j,2]*nmax!=1.:
            
            dwell_J_075_k1[i] = j-S_075*i
            break
 
#%% Dwell time J=1 for k=0, k=1

for i in range(sim):
    for j in range(S_1*i,S_1*(i+1)-1):
        if All_Res_J_1[j,2]*nmax==1.:
            dwell_J_1_k0[i] = j-S_1*i
            break
        
for i in range(sim):     
    for j in range(int(S_1*i+dwell_J_1_k0[i]),S_1*(i+1)-1):
        if All_Res_J_1[j,2]*nmax!=1.:
            
            dwell_J_1_k1[i] = j-S_1*i
            break
        
for i in range(sim):
    for j in range(S_125*i,S_125*(i+1)-1):
        if All_Res_J_125[j,2]*nmax==1.:
            dwell_J_125_k0[i] = j-S_125*i
            break
        
for i in range(sim):     
    for j in range(int(S_125*i+dwell_J_125_k0[i]),S_125*(i+1)-1):
        if All_Res_J_125[j,2]*nmax!=1.:
            
            dwell_J_125_k1[i] = j-S_125*i
            break

for i in range(sim):
    for j in range(S_150*i,S_150*(i+1)-1):
        if All_Res_J_150[j,2]*nmax==1.:
            dwell_J_150_k0[i] = j-S_150*i
            break

for i in range(sim):     
    for j in range(int(S_150*i+dwell_J_150_k0[i]),S_150*(i+1)-1):
        if All_Res_J_150[j,2]*nmax!=1.:
            
            dwell_J_150_k1[i] = j-S_150*i
            break

for i in range(sim):
    for j in range(S_175*i,S_175*(i+1)-1):
        if All_Res_J_175[j,2]*nmax==1.:
            dwell_J_175_k0[i] = j-S_175*i
            break

for i in range(sim):     
    for j in range(int(S_175*i+dwell_J_175_k0[i]),S_175*(i+1)-1):
        if All_Res_J_175[j,2]*nmax!=1.:
            
            dwell_J_175_k1[i] = j-S_175*i
            break

#%% Dwell time J=2 for k=0,k=1
        
for i in range(sim):
    for j in range(S_2*i,S_2*(i+1)-1):
        if All_Res_J_2[j,2]*nmax==1.:
            dwell_J_2_k0[i] = j-S_2*i
            break
        
for i in range(sim):     
    for j in range(int(S_2*i+dwell_J_2_k0[i]),S_2*(i+1)-1):
        if All_Res_J_2[j,2]*nmax!=1.:
            
            dwell_J_2_k1[i] = j-S_2*i
            break
        
for i in range(sim):
    for j in range(S_225*i,S_225*(i+1)-1):
        if All_Res_J_225[j,2]*nmax==1.:
            dwell_J_225_k0[i] = j-S_225*i
            break
        
for i in range(sim):     
    for j in range(int(S_225*i+dwell_J_225_k0[i]),S_225*(i+1)-1):
        if All_Res_J_225[j,2]*nmax!=1.:
            
            dwell_J_225_k1[i] = j-S_225*i
            break
        
for i in range(sim):
    for j in range(S_250*i,S_250*(i+1)-1):
        if All_Res_J_250[j,2]*nmax==1.:
            dwell_J_250_k0[i] = j-S_250*i
            break
        
for i in range(sim):     
    for j in range(int(S_250*i+dwell_J_250_k0[i]),S_250*(i+1)-1):
        if All_Res_J_250[j,2]*nmax!=1.:
            
            dwell_J_250_k1[i] = j-S_250*i
            break

for i in range(sim):
    for j in range(S_275*i,S_275*(i+1)-1):
        if All_Res_J_275[j,2]*nmax==1.:
            dwell_J_275_k0[i] = j-S_275*i
            break
        
for i in range(sim):     
    for j in range(int(S_275*i+dwell_J_275_k0[i]),S_275*(i+1)-1):
        if All_Res_J_275[j,2]*nmax!=1.:
            
            dwell_J_275_k1[i] = j-S_275*i
            break

#%% Dwell time J=3 for k=0,k=1

for i in range(sim):
    for j in range(S_3*i,S_3*(i+1)-1):
        if All_Res_J_3[j,2]*nmax==1.:
            dwell_J_3_k0[i] = j-S_3*i
            break
        
for i in range(sim):     
    for j in range(int(S_3*i+dwell_J_3_k0[i]),S_3*(i+1)-1):
        if All_Res_J_3[j,2]*nmax!=1.:
            
            dwell_J_3_k1[i] = j-S_3*i
            break
        
for i in range(sim):
    for j in range(S_325*i,S_325*(i+1)-1):
        if All_Res_J_325[j,2]*nmax==1.:
            dwell_J_325_k0[i] = j-S_325*i
            break
        
for i in range(sim):     
    for j in range(int(S_325*i+dwell_J_325_k0[i]),S_325*(i+1)-1):
        if All_Res_J_325[j,2]*nmax!=1.:
            
            dwell_J_325_k1[i] = j-S_325*i
            break
        
for i in range(sim):
    for j in range(S_350*i,S_350*(i+1)-1):
        if All_Res_J_350[j,2]*nmax==1.:
            dwell_J_350_k0[i] = j-S_350*i
            break
        
for i in range(sim):     
    for j in range(int(S_350*i+dwell_J_350_k0[i]),S_350*(i+1)-1):
        if All_Res_J_350[j,2]*nmax!=1.:
            
            dwell_J_350_k1[i] = j-S_350*i
            break
        
for i in range(sim):
    for j in range(S_375*i,S_375*(i+1)-1):
        if All_Res_J_375[j,2]*nmax==1.:
            dwell_J_375_k0[i] = j-S_375*i
            break
        
for i in range(sim):     
    for j in range(int(S_375*i+dwell_J_375_k0[i]),S_375*(i+1)-1):
        if All_Res_J_375[j,2]*nmax!=1.:
            
            dwell_J_375_k1[i] = j-S_375*i
            break
        
#%% Dwell time J=4 for k=0,k=1
        
for i in range(sim):
    for j in range(S_4*i,S_4*(i+1)-1):
        if All_Res_J_4[j,2]*nmax==1.:
            dwell_J_4_k0[i] = j-S_4*i
            break
        
for i in range(sim):     
    for j in range(int(S_4*i+dwell_J_4_k0[i]),S_4*(i+1)-1):
        if All_Res_J_4[j,2]*nmax!=1.:
            
            dwell_J_4_k1[i] = j-S_4*i
            break
        
for i in range(sim):
    for j in range(S_425*i,S_425*(i+1)-1):
        if All_Res_J_425[j,2]*nmax==1.:
            dwell_J_425_k0[i] = j-S_425*i
            break
        
for i in range(sim):     
    for j in range(int(S_425*i+dwell_J_425_k0[i]),S_425*(i+1)-1):
        if All_Res_J_425[j,2]*nmax!=1.:
            
            dwell_J_425_k1[i] = j-S_425*i
            break
        
for i in range(sim):
    for j in range(S_450*i,S_450*(i+1)-1):
        if All_Res_J_450[j,2]*nmax==1.:
            dwell_J_450_k0[i] = j-S_450*i
            break
        
for i in range(sim):     
    for j in range(int(S_450*i+dwell_J_450_k0[i]),S_450*(i+1)-1):
        if All_Res_J_450[j,2]*nmax!=1.:
            
            dwell_J_450_k1[i] = j-S_450*i
            break
        
for i in range(sim):
    for j in range(S_475*i,S_475*(i+1)-1):
        if All_Res_J_475[j,2]*nmax==1.:
            dwell_J_475_k0[i] = j-S_475*i
            break
        
for i in range(sim):     
    for j in range(int(S_475*i+dwell_J_475_k0[i]),S_475*(i+1)-1):
        if All_Res_J_475[j,2]*nmax!=1.:
            
            dwell_J_475_k1[i] = j-S_475*i
            break
        
#%% Dwell time J=5 for k=0,k=1

for i in range(sim):
    for j in range(S_5*i,S_5*(i+1)-1):
        if All_Res_J_5[j,2]*nmax==1.:
            dwell_J_5_k0[i] = j-S_5*i
            break
        
for i in range(sim):     
    for j in range(int(S_5*i+dwell_J_5_k0[i]),S_5*(i+1)-1):
        if All_Res_J_5[j,2]*nmax!=1.:
            
            dwell_J_5_k1[i] = j-S_5*i
            break

#%% Average and sd dwell times J=0

dwell_avg_J_0_k0 = sum(dwell_J_0_k0)/sim
dwell_avg2_J_0_k0 = sum(dwell_J_0_k0**2)/sim
dwell_sd_J_0_k0 = np.sqrt(dwell_avg2_J_0_k0 - dwell_avg_J_0_k0**2)

dwell_avg_J_0_k1 = sum(dwell_J_0_k1)/sim
dwell_avg2_J_0_k1 = sum(dwell_J_0_k1**2)/sim
dwell_sd_J_0_k1 = np.sqrt(dwell_avg2_J_0_k1 - dwell_avg_J_0_k1**2)


dwell_avg_J_025_k0 = sum(dwell_J_025_k0)/sim
dwell_avg2_J_025_k0 = sum(dwell_J_025_k0**2)/sim
dwell_sd_J_025_k0 = np.sqrt(dwell_avg2_J_025_k0 - dwell_avg_J_025_k0**2)

dwell_avg_J_025_k1 = sum(dwell_J_025_k1)/sim
dwell_avg2_J_025_k1 = sum(dwell_J_025_k1**2)/sim
dwell_sd_J_025_k1 = np.sqrt(dwell_avg2_J_025_k1 - dwell_avg_J_025_k1**2)


dwell_avg_J_050_k0 = sum(dwell_J_050_k0)/sim
dwell_avg2_J_050_k0 = sum(dwell_J_050_k0**2)/sim
dwell_sd_J_050_k0 = np.sqrt(dwell_avg2_J_050_k0 - dwell_avg_J_050_k0**2)

dwell_avg_J_050_k1 = sum(dwell_J_050_k1)/sim
dwell_avg2_J_050_k1 = sum(dwell_J_050_k1**2)/sim
dwell_sd_J_050_k1 = np.sqrt(dwell_avg2_J_050_k1 - dwell_avg_J_050_k1**2)


dwell_avg_J_075_k0 = sum(dwell_J_075_k0)/sim
dwell_avg2_J_075_k0 = sum(dwell_J_075_k0**2)/sim
dwell_sd_J_075_k0 = np.sqrt(dwell_avg2_J_075_k0 - dwell_avg_J_075_k0**2)

dwell_avg_J_075_k1 = sum(dwell_J_075_k1)/sim
dwell_avg2_J_075_k1 = sum(dwell_J_075_k1**2)/sim
dwell_sd_J_075_k1 = np.sqrt(dwell_avg2_J_075_k1 - dwell_avg_J_075_k1**2)

#%% Average and sd dwell times J=1

dwell_avg_J_1_k0 = sum(dwell_J_1_k0)/sim
dwell_avg2_J_1_k0 = sum(dwell_J_1_k0**2)/sim
dwell_sd_J_1_k0 = np.sqrt(dwell_avg2_J_1_k0 - dwell_avg_J_1_k0**2)

dwell_avg_J_1_k1 = sum(dwell_J_1_k1)/sim
dwell_avg2_J_1_k1 = sum(dwell_J_1_k1**2)/sim
dwell_sd_J_1_k1 = np.sqrt(dwell_avg2_J_1_k1 - dwell_avg_J_1_k1**2)


dwell_avg_J_125_k0 = sum(dwell_J_125_k0)/sim
dwell_avg2_J_125_k0 = sum(dwell_J_125_k0**2)/sim
dwell_sd_J_125_k0 = np.sqrt(dwell_avg2_J_125_k0 - dwell_avg_J_125_k0**2)

dwell_avg_J_125_k1 = sum(dwell_J_125_k1)/sim
dwell_avg2_J_125_k1 = sum(dwell_J_125_k1**2)/sim
dwell_sd_J_125_k1 = np.sqrt(dwell_avg2_J_125_k1 - dwell_avg_J_125_k1**2)


dwell_avg_J_150_k0 = sum(dwell_J_150_k0)/sim
dwell_avg2_J_150_k0 = sum(dwell_J_150_k0**2)/sim
dwell_sd_J_150_k0 = np.sqrt(dwell_avg2_J_150_k0 - dwell_avg_J_150_k0**2)

dwell_avg_J_150_k1 = sum(dwell_J_150_k1)/sim
dwell_avg2_J_150_k1 = sum(dwell_J_150_k1**2)/sim
dwell_sd_J_150_k1 = np.sqrt(dwell_avg2_J_150_k1 - dwell_avg_J_150_k1**2)


dwell_avg_J_175_k0 = sum(dwell_J_175_k0)/sim
dwell_avg2_J_175_k0 = sum(dwell_J_175_k0**2)/sim
dwell_sd_J_175_k0 = np.sqrt(dwell_avg2_J_175_k0 - dwell_avg_J_175_k0**2)

dwell_avg_J_175_k1 = sum(dwell_J_175_k1)/sim
dwell_avg2_J_175_k1 = sum(dwell_J_175_k1**2)/sim
dwell_sd_J_175_k1 = np.sqrt(dwell_avg2_J_175_k1 - dwell_avg_J_175_k1**2)

#%% Average and sd dwell times J=2

dwell_avg_J_2_k0 = sum(dwell_J_2_k0)/sim
dwell_avg2_J_2_k0 = sum(dwell_J_2_k0**2)/sim
dwell_sd_J_2_k0 = np.sqrt(dwell_avg2_J_2_k0 - dwell_avg_J_2_k0**2)

dwell_avg_J_2_k1 = sum(dwell_J_2_k1)/sim
dwell_avg2_J_2_k1 = sum(dwell_J_2_k1**2)/sim
dwell_sd_J_2_k1 = np.sqrt(dwell_avg2_J_2_k1 - dwell_avg_J_2_k1**2)


dwell_avg_J_225_k0 = sum(dwell_J_225_k0)/sim
dwell_avg2_J_225_k0 = sum(dwell_J_225_k0**2)/sim
dwell_sd_J_225_k0 = np.sqrt(dwell_avg2_J_225_k0 - dwell_avg_J_225_k0**2)

dwell_avg_J_225_k1 = sum(dwell_J_225_k1)/sim
dwell_avg2_J_225_k1 = sum(dwell_J_225_k1**2)/sim
dwell_sd_J_225_k1 = np.sqrt(dwell_avg2_J_225_k1 - dwell_avg_J_225_k1**2)


dwell_avg_J_250_k0 = sum(dwell_J_250_k0)/sim
dwell_avg2_J_250_k0 = sum(dwell_J_250_k0**2)/sim
dwell_sd_J_250_k0 = np.sqrt(dwell_avg2_J_250_k0 - dwell_avg_J_250_k0**2)

dwell_avg_J_250_k1 = sum(dwell_J_250_k1)/sim
dwell_avg2_J_250_k1 = sum(dwell_J_250_k1**2)/sim
dwell_sd_J_250_k1 = np.sqrt(dwell_avg2_J_250_k1 - dwell_avg_J_250_k1**2)


dwell_avg_J_275_k0 = sum(dwell_J_275_k0)/sim
dwell_avg2_J_275_k0 = sum(dwell_J_275_k0**2)/sim
dwell_sd_J_275_k0 = np.sqrt(dwell_avg2_J_275_k0 - dwell_avg_J_275_k0**2)

dwell_avg_J_275_k1 = sum(dwell_J_275_k1)/sim
dwell_avg2_J_275_k1 = sum(dwell_J_275_k1**2)/sim
dwell_sd_J_275_k1 = np.sqrt(dwell_avg2_J_275_k1 - dwell_avg_J_275_k1**2)


#%% Average and sd dwell times J=3

dwell_avg_J_3_k0 = sum(dwell_J_3_k0)/sim
dwell_avg2_J_3_k0 = sum(dwell_J_3_k0**2)/sim
dwell_sd_J_3_k0 = np.sqrt(dwell_avg2_J_3_k0 - dwell_avg_J_3_k0**2)

dwell_avg_J_3_k1 = sum(dwell_J_3_k1)/sim
dwell_avg2_J_3_k1 = sum(dwell_J_3_k1**2)/sim
dwell_sd_J_3_k1 = np.sqrt(dwell_avg2_J_3_k1 - dwell_avg_J_3_k1**2)


dwell_avg_J_325_k0 = sum(dwell_J_325_k0)/sim
dwell_avg2_J_325_k0 = sum(dwell_J_325_k0**2)/sim
dwell_sd_J_325_k0 = np.sqrt(dwell_avg2_J_325_k0 - dwell_avg_J_325_k0**2)

dwell_avg_J_325_k1 = sum(dwell_J_325_k1)/sim
dwell_avg2_J_325_k1 = sum(dwell_J_325_k1**2)/sim
dwell_sd_J_325_k1 = np.sqrt(dwell_avg2_J_325_k1 - dwell_avg_J_325_k1**2)


dwell_avg_J_350_k0 = sum(dwell_J_350_k0)/sim
dwell_avg2_J_350_k0 = sum(dwell_J_350_k0**2)/sim
dwell_sd_J_350_k0 = np.sqrt(dwell_avg2_J_350_k0 - dwell_avg_J_350_k0**2)

dwell_avg_J_350_k1 = sum(dwell_J_350_k1)/sim
dwell_avg2_J_350_k1 = sum(dwell_J_350_k1**2)/sim
dwell_sd_J_350_k1 = np.sqrt(dwell_avg2_J_350_k1 - dwell_avg_J_350_k1**2)


dwell_avg_J_375_k0 = sum(dwell_J_375_k0)/sim
dwell_avg2_J_375_k0 = sum(dwell_J_375_k0**2)/sim
dwell_sd_J_375_k0 = np.sqrt(dwell_avg2_J_375_k0 - dwell_avg_J_375_k0**2)

dwell_avg_J_375_k1 = sum(dwell_J_375_k1)/sim
dwell_avg2_J_375_k1 = sum(dwell_J_375_k1**2)/sim
dwell_sd_J_375_k1 = np.sqrt(dwell_avg2_J_375_k1 - dwell_avg_J_375_k1**2)

#%% Aaverage and sd dwell times J=4

dwell_avg_J_4_k0 = sum(dwell_J_4_k0)/sim
dwell_avg2_J_4_k0 = sum(dwell_J_4_k0**2)/sim
dwell_sd_J_4_k0 = np.sqrt(dwell_avg2_J_4_k0 - dwell_avg_J_4_k0**2)

dwell_avg_J_4_k1 = sum(dwell_J_4_k1)/sim
dwell_avg2_J_4_k1 = sum(dwell_J_4_k1**2)/sim
dwell_sd_J_4_k1 = np.sqrt(dwell_avg2_J_4_k1 - dwell_avg_J_4_k1**2)


dwell_avg_J_425_k0 = sum(dwell_J_425_k0)/sim
dwell_avg2_J_425_k0 = sum(dwell_J_425_k0**2)/sim
dwell_sd_J_425_k0 = np.sqrt(dwell_avg2_J_425_k0 - dwell_avg_J_425_k0**2)

dwell_avg_J_425_k1 = sum(dwell_J_425_k1)/sim
dwell_avg2_J_425_k1 = sum(dwell_J_425_k1**2)/sim
dwell_sd_J_425_k1 = np.sqrt(dwell_avg2_J_425_k1 - dwell_avg_J_425_k1**2)


dwell_avg_J_450_k0 = sum(dwell_J_450_k0)/sim
dwell_avg2_J_450_k0 = sum(dwell_J_450_k0**2)/sim
dwell_sd_J_450_k0 = np.sqrt(dwell_avg2_J_450_k0 - dwell_avg_J_450_k0**2)

dwell_avg_J_450_k1 = sum(dwell_J_450_k1)/sim
dwell_avg2_J_450_k1 = sum(dwell_J_450_k1**2)/sim
dwell_sd_J_450_k1 = np.sqrt(dwell_avg2_J_450_k1 - dwell_avg_J_450_k1**2)


dwell_avg_J_475_k0 = sum(dwell_J_475_k0)/sim
dwell_avg2_J_475_k0 = sum(dwell_J_475_k0**2)/sim
dwell_sd_J_475_k0 = np.sqrt(dwell_avg2_J_475_k0 - dwell_avg_J_475_k0**2)

dwell_avg_J_475_k1 = sum(dwell_J_475_k1)/sim
dwell_avg2_J_475_k1 = sum(dwell_J_475_k1**2)/sim
dwell_sd_J_475_k1 = np.sqrt(dwell_avg2_J_475_k1 - dwell_avg_J_475_k1**2)

#%% Average and sd dwell times J=5

dwell_avg_J_5_k0 = sum(dwell_J_5_k0)/sim
dwell_avg2_J_5_k0 = sum(dwell_J_5_k0**2)/sim
dwell_sd_J_5_k0 = np.sqrt(dwell_avg2_J_5_k0 - dwell_avg_J_5_k0**2)

dwell_avg_J_5_k1 = sum(dwell_J_5_k1)/sim
dwell_avg2_J_5_k1 = sum(dwell_J_5_k1**2)/sim
dwell_sd_J_5_k1 = np.sqrt(dwell_avg2_J_5_k1 - dwell_avg_J_5_k1**2)

#%%

J=np.arange(0,5.25,0.25)

dwell_avg_k0 = np.array([dwell_avg_J_0_k0,dwell_avg_J_025_k0,dwell_avg_J_050_k0,dwell_avg_J_075_k0,dwell_avg_J_1_k0,
                         dwell_avg_J_125_k0,dwell_avg_J_150_k0,dwell_avg_J_175_k0,dwell_avg_J_2_k0,dwell_avg_J_225_k0,
                         dwell_avg_J_250_k0,dwell_avg_J_275_k0,dwell_avg_J_3_k0,dwell_avg_J_325_k0,dwell_avg_J_350_k0,
                         dwell_avg_J_375_k0,dwell_avg_J_4_k0,dwell_avg_J_425_k0,dwell_avg_J_450_k0,dwell_avg_J_475_k0,
                         dwell_avg_J_5_k0])

dwell_sd_k0 = np.array([dwell_sd_J_0_k0,dwell_sd_J_025_k0,dwell_sd_J_050_k0,dwell_sd_J_075_k0,dwell_sd_J_1_k0,
                        dwell_sd_J_125_k0,dwell_sd_J_150_k0,dwell_sd_J_175_k0,dwell_sd_J_2_k0,dwell_sd_J_225_k0,
                        dwell_sd_J_250_k0,dwell_sd_J_275_k0,dwell_sd_J_3_k0,dwell_sd_J_325_k0,dwell_sd_J_350_k0,
                        dwell_sd_J_375_k0,dwell_sd_J_4_k0,dwell_sd_J_425_k0,dwell_sd_J_450_k0,dwell_sd_J_475_k0,
                        dwell_sd_J_5_k0])


dwell_avg_k1 = np.array([dwell_avg_J_0_k1,dwell_avg_J_025_k1,dwell_avg_J_050_k1,dwell_avg_J_075_k1,dwell_avg_J_1_k1,
                         dwell_avg_J_125_k1,dwell_avg_J_150_k1,dwell_avg_J_175_k1,dwell_avg_J_2_k1,dwell_avg_J_225_k1,
                         dwell_avg_J_250_k1,dwell_avg_J_275_k1,dwell_avg_J_3_k1,dwell_avg_J_325_k1,dwell_avg_J_350_k1,
                         dwell_avg_J_375_k1,dwell_avg_J_4_k1,dwell_avg_J_425_k1,dwell_avg_J_450_k1,dwell_avg_J_475_k1,
                         dwell_avg_J_5_k1])

dwell_sd_k1 = np.array([dwell_sd_J_0_k1,dwell_sd_J_025_k1,dwell_sd_J_050_k1,dwell_sd_J_075_k1,dwell_sd_J_1_k1,
                        dwell_sd_J_125_k1,dwell_sd_J_150_k1,dwell_sd_J_175_k1,dwell_sd_J_2_k1,dwell_sd_J_225_k1,
                        dwell_sd_J_250_k1,dwell_sd_J_275_k1,dwell_sd_J_3_k1,dwell_sd_J_325_k1,dwell_sd_J_350_k1,
                        dwell_sd_J_375_k1,dwell_sd_J_4_k1,dwell_sd_J_425_k1,dwell_sd_J_450_k1,dwell_sd_J_475_k1,
                        dwell_sd_J_5_k1])

#%%

plt.rcParams["figure.figsize"] = [8.0,6.0]

plt.title('$<\phi>=1/2$',fontsize=17)
plt.xlabel('J ($k_B T$)',fontsize=14)
plt.ylabel('Dwell time',fontsize=14)
plt.errorbar(J,dwell_avg_k0,dwell_sd_k0,marker='o',linestyle='',capsize=3,label='k=0')
plt.errorbar(J,dwell_avg_k1,dwell_sd_k1,marker='o',capsize=3,linestyle='',label='k=1')

plt.grid()

plt.legend()

#%%

np.savetxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Dwell_times/dwell_times_mu-1_k0.dat',
           (np.c_[J,dwell_avg_k0,dwell_sd_k0]))

np.savetxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Dwell_times/dwell_times_mu-1_k1.dat',
           (np.c_[J,dwell_avg_k1,dwell_sd_k1]))

#%%

DW_mu_1_k0 = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Dwell_times/dwell_times_mu-1_k0.dat')
DW_mu_1_k1 = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Dwell_times/dwell_times_mu-1_k1.dat')

DW_mu_2_k0 = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Dwell_times/dwell_times_mu-2_k0.dat')
DW_mu_2_k1 = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Dwell_times/dwell_times_mu-2_k1.dat')

DW_mu_3_k0 = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Dwell_times/dwell_times_mu-3_k0.dat')
DW_mu_3_k1 = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Dwell_times/dwell_times_mu-3_k1.dat')

#%%

plt.xlabel('J ($k_BT$)', fontsize=14)
plt.ylabel('Dwell time',fontsize=14)

plt.errorbar(DW_mu_1_k0[:,0], DW_mu_1_k0[:,1], DW_mu_1_k0[:,2], linestyle='',capsize=3,marker='o',color='blue')
plt.errorbar(DW_mu_1_k1[:,0], DW_mu_1_k1[:,1], DW_mu_1_k1[:,2], linestyle='',capsize=3,marker='^',color='blue')

plt.errorbar(DW_mu_2_k0[:,0], DW_mu_2_k0[:,1], DW_mu_2_k0[:,2], linestyle='',capsize=3,marker='o',color='orange')
plt.errorbar(DW_mu_2_k1[:,0], DW_mu_2_k1[:,1], DW_mu_2_k1[:,2], linestyle='',capsize=3,marker='^',color='orange')

plt.errorbar(DW_mu_3_k0[:,0], DW_mu_3_k0[:,1], DW_mu_3_k0[:,2], linestyle='',capsize=3,marker='o',color='green')
plt.errorbar(DW_mu_3_k1[:,0], DW_mu_3_k1[:,1], DW_mu_3_k1[:,2], linestyle='',capsize=3,marker='^',color='green')


plt.grid()


k0 = mlines.Line2D([],[],linestyle='', color='black', marker='o', label='Dwell time k=0')
k1 = mlines.Line2D([],[],linestyle='', color='black', marker='^', label='Dwell time k=1')
blue_patch = mpatches.Patch(color='blue', label='$\mu=-1$')
orange_patch = mpatches.Patch(color='orange', label='$\mu = -2$')
green_patch = mpatches.Patch(color='green', label='$\mu = -3$')

plt.legend(loc='upper right', prop={'size': 15},handles=[k0,k1,blue_patch,orange_patch,green_patch])

#%%


k=np.arange(0,5,1)

k_J0 = np.arange(0,14,1)

dwell_J_0_k=np.zeros([len(k_J0),sim])
dwell_J_1_k=np.zeros([len(k),sim])
dwell_J_2_k=np.zeros([len(k),sim])
dwell_J_3_k=np.zeros([len(k),sim])
dwell_J_4_k=np.zeros([len(k),sim])
dwell_J_5_k=np.zeros([len(k),sim])

for n in range(len(k)):
    for i in range(sim):
#        for j in range(S_0*i,S_0*(i+1)-1):
#            if All_Res_J_0[j,2]*nmax==n:
#                dwell_J_0_k[n,i] = dwell_J_0_k[n,i] + 1.
                
        for j in range(S_1*i,S_1*(i+1)-1):
            if All_Res_J_1[j,2]*nmax==n:
                dwell_J_1_k[n,i] = dwell_J_1_k[n,i] + 1.
                
        for j in range(S_2*i,S_2*(i+1)-1):
            if All_Res_J_2[j,2]*nmax==n:
                dwell_J_2_k[n,i] = dwell_J_2_k[n,i] + 1.
                
        for j in range(S_3*i,S_3*(i+1)-1):
            if All_Res_J_3[j,2]*nmax==n:
                dwell_J_3_k[n,i] = dwell_J_3_k[n,i] + 1.
        
        for j in range(S_4*i,S_4*(i+1)-1):
            if All_Res_J_4[j,2]*nmax==n:
                dwell_J_4_k[n,i] = dwell_J_4_k[n,i] + 1.
                
        for j in range(S_5*i,S_5*(i+1)-1):
            if All_Res_J_5[j,2]*nmax==n:
                dwell_J_5_k[n,i] = dwell_J_5_k[n,i] + 1.

 

for n in range(len(k_J0)):
    for i in range(sim):
        for j in range(S_0*i,S_0*(i+1)-1):
            if All_Res_J_0[j,2]*nmax==n:
                dwell_J_0_k[n,i] = dwell_J_0_k[n,i] + 1.
       
#%%

dwell_avg_J_0_k0 = sum(dwell_J_0_k[0,:])/sim
dwell_avg2_J_0_k0 = sum(dwell_J_0_k[0,:]**2)/sim
dwell_sd_J_0_k0 = np.sqrt(dwell_avg2_J_0_k0-dwell_avg_J_0_k0**2)

dwell_avg_J_0_k1 = sum(dwell_J_0_k[1,:])/sim
dwell_avg2_J_0_k1 = sum(dwell_J_0_k[1,:]**2)/sim
dwell_sd_J_0_k1 = np.sqrt(dwell_avg2_J_0_k1-dwell_avg_J_0_k1**2)

dwell_avg_J_0_k2 = sum(dwell_J_0_k[2,:])/sim
dwell_avg2_J_0_k2 = sum(dwell_J_0_k[2,:]**2)/sim
dwell_sd_J_0_k2 = np.sqrt(dwell_avg2_J_0_k2-dwell_avg_J_0_k2**2)

dwell_avg_J_0_k3 = sum(dwell_J_0_k[3,:])/sim
dwell_avg2_J_0_k3 = sum(dwell_J_0_k[3,:]**2)/sim
dwell_sd_J_0_k3 = np.sqrt(dwell_avg2_J_0_k3-dwell_avg_J_0_k3**2)

dwell_avg_J_0_k4 = sum(dwell_J_0_k[4,:])/sim
dwell_avg2_J_0_k4 = sum(dwell_J_0_k[4,:]**2)/sim
dwell_sd_J_0_k4 = np.sqrt(dwell_avg2_J_0_k4-dwell_avg_J_0_k4**2)

dwell_avg_J_0_k5 = sum(dwell_J_0_k[5,:])/sim
dwell_avg2_J_0_k5 = sum(dwell_J_0_k[5,:]**2)/sim
dwell_sd_J_0_k5 = np.sqrt(dwell_avg2_J_0_k5-dwell_avg_J_0_k5**2)

dwell_avg_J_0_k6 = sum(dwell_J_0_k[6,:])/sim
dwell_avg2_J_0_k6 = sum(dwell_J_0_k[6,:]**2)/sim
dwell_sd_J_0_k6 = np.sqrt(dwell_avg2_J_0_k6-dwell_avg_J_0_k6**2)

dwell_avg_J_0_k7 = sum(dwell_J_0_k[7,:])/sim
dwell_avg2_J_0_k7 = sum(dwell_J_0_k[7,:]**2)/sim
dwell_sd_J_0_k7 = np.sqrt(dwell_avg2_J_0_k7-dwell_avg_J_0_k7**2)

dwell_avg_J_0_k8 = sum(dwell_J_0_k[8,:])/sim
dwell_avg2_J_0_k8 = sum(dwell_J_0_k[8,:]**2)/sim
dwell_sd_J_0_k8 = np.sqrt(dwell_avg2_J_0_k8-dwell_avg_J_0_k8**2)

dwell_avg_J_0_k9 = sum(dwell_J_0_k[9,:])/sim
dwell_avg2_J_0_k9 = sum(dwell_J_0_k[9,:]**2)/sim
dwell_sd_J_0_k9 = np.sqrt(dwell_avg2_J_0_k9-dwell_avg_J_0_k9**2)

dwell_avg_J_0_k10 = sum(dwell_J_0_k[10,:])/sim
dwell_avg2_J_0_k10 = sum(dwell_J_0_k[10,:]**2)/sim
dwell_sd_J_0_k10 = np.sqrt(dwell_avg2_J_0_k10-dwell_avg_J_0_k10**2)

dwell_avg_J_0_k11 = sum(dwell_J_0_k[11,:])/sim
dwell_avg2_J_0_k11 = sum(dwell_J_0_k[11,:]**2)/sim
dwell_sd_J_0_k11 = np.sqrt(dwell_avg2_J_0_k11-dwell_avg_J_0_k11**2)

dwell_avg_J_0_k12 = sum(dwell_J_0_k[12,:])/sim
dwell_avg2_J_0_k12 = sum(dwell_J_0_k[12,:]**2)/sim
dwell_sd_J_0_k12 = np.sqrt(dwell_avg2_J_0_k12-dwell_avg_J_0_k12**2)

dwell_avg_J_0_k13 = sum(dwell_J_0_k[13,:])/sim
dwell_avg2_J_0_k13 = sum(dwell_J_0_k[13,:]**2)/sim
dwell_sd_J_0_k13 = np.sqrt(dwell_avg2_J_0_k13-dwell_avg_J_0_k13**2)



dwell_avg_J_1_k0 = sum(dwell_J_1_k[0,:])/sim
dwell_avg2_J_1_k0 = sum(dwell_J_1_k[0,:]**2)/sim
dwell_sd_J_1_k0 = np.sqrt(dwell_avg2_J_1_k0-dwell_avg_J_1_k0**2)

dwell_avg_J_1_k1 = sum(dwell_J_1_k[1,:])/sim
dwell_avg2_J_1_k1 = sum(dwell_J_1_k[1,:]**2)/sim
dwell_sd_J_1_k1 = np.sqrt(dwell_avg2_J_1_k1-dwell_avg_J_1_k1**2)

dwell_avg_J_1_k2 = sum(dwell_J_1_k[2,:])/sim
dwell_avg2_J_1_k2 = sum(dwell_J_1_k[2,:]**2)/sim
dwell_sd_J_1_k2 = np.sqrt(dwell_avg2_J_1_k2-dwell_avg_J_1_k2**2)

dwell_avg_J_1_k3 = sum(dwell_J_1_k[3,:])/sim
dwell_avg2_J_1_k3 = sum(dwell_J_1_k[3,:]**2)/sim
dwell_sd_J_1_k3 = np.sqrt(dwell_avg2_J_1_k3-dwell_avg_J_1_k3**2)

dwell_avg_J_1_k4 = sum(dwell_J_1_k[4,:])/sim
dwell_avg2_J_1_k4 = sum(dwell_J_1_k[4,:]**2)/sim
dwell_sd_J_1_k4 = np.sqrt(dwell_avg2_J_1_k4-dwell_avg_J_1_k4**2)


dwell_avg_J_2_k0 = sum(dwell_J_2_k[0,:])/sim
dwell_avg2_J_2_k0 = sum(dwell_J_2_k[0,:]**2)/sim
dwell_sd_J_2_k0 = np.sqrt(dwell_avg2_J_2_k0-dwell_avg_J_2_k0**2)

dwell_avg_J_2_k1 = sum(dwell_J_2_k[1,:])/sim
dwell_avg2_J_2_k1 = sum(dwell_J_2_k[1,:]**2)/sim
dwell_sd_J_2_k1 = np.sqrt(dwell_avg2_J_2_k1-dwell_avg_J_2_k1**2)

dwell_avg_J_2_k2 = sum(dwell_J_2_k[2,:])/sim
dwell_avg2_J_2_k2 = sum(dwell_J_2_k[2,:]**2)/sim
dwell_sd_J_2_k2 = np.sqrt(dwell_avg2_J_2_k2-dwell_avg_J_2_k2**2)

dwell_avg_J_2_k3 = sum(dwell_J_2_k[3,:])/sim
dwell_avg2_J_2_k3 = sum(dwell_J_2_k[3,:]**2)/sim
dwell_sd_J_2_k3 = np.sqrt(dwell_avg2_J_2_k3-dwell_avg_J_2_k3**2)

dwell_avg_J_2_k4 = sum(dwell_J_2_k[4,:])/sim
dwell_avg2_J_2_k4 = sum(dwell_J_2_k[4,:]**2)/sim
dwell_sd_J_2_k4 = np.sqrt(dwell_avg2_J_2_k4-dwell_avg_J_2_k4**2)


dwell_avg_J_3_k0 = sum(dwell_J_3_k[0,:])/sim
dwell_avg2_J_3_k0 = sum(dwell_J_3_k[0,:]**2)/sim
dwell_sd_J_3_k0 = np.sqrt(dwell_avg2_J_3_k0-dwell_avg_J_3_k0**2)

dwell_avg_J_3_k1 = sum(dwell_J_3_k[1,:])/sim
dwell_avg2_J_3_k1 = sum(dwell_J_3_k[1,:]**2)/sim
dwell_sd_J_3_k1 = np.sqrt(dwell_avg2_J_3_k1-dwell_avg_J_3_k1**2)

dwell_avg_J_3_k2 = sum(dwell_J_3_k[2,:])/sim
dwell_avg2_J_3_k2 = sum(dwell_J_3_k[2,:]**2)/sim
dwell_sd_J_3_k2 = np.sqrt(dwell_avg2_J_3_k2-dwell_avg_J_3_k2**2)

dwell_avg_J_3_k3 = sum(dwell_J_3_k[3,:])/sim
dwell_avg2_J_3_k3 = sum(dwell_J_3_k[3,:]**2)/sim
dwell_sd_J_3_k3 = np.sqrt(dwell_avg2_J_3_k3-dwell_avg_J_3_k3**2)

dwell_avg_J_3_k4 = sum(dwell_J_3_k[4,:])/sim
dwell_avg2_J_3_k4 = sum(dwell_J_3_k[4,:]**2)/sim
dwell_sd_J_3_k4 = np.sqrt(dwell_avg2_J_3_k4-dwell_avg_J_3_k4**2)


dwell_avg_J_4_k0 = sum(dwell_J_4_k[0,:])/sim

dwell_avg_J_4_k1 = sum(dwell_J_4_k[1,:])/sim

dwell_avg_J_4_k2 = sum(dwell_J_4_k[2,:])/sim

dwell_avg_J_4_k3 = sum(dwell_J_4_k[3,:])/sim

dwell_avg_J_4_k4 = sum(dwell_J_4_k[4,:])/sim


dwell_avg_J_5_k0 = sum(dwell_J_0_k[0,:])/sim

dwell_avg_J_5_k1 = sum(dwell_J_0_k[1,:])/sim

dwell_avg_J_5_k2 = sum(dwell_J_0_k[2,:])/sim

dwell_avg_J_5_k3 = sum(dwell_J_0_k[3,:])/sim

dwell_avg_J_5_k4 = sum(dwell_J_0_k[4,:])/sim



#dwell_avg_J_0_k = np.array([dwell_avg_J_0_k0,dwell_avg_J_0_k1,dwell_avg_J_0_k2,dwell_avg_J_0_k3,dwell_avg_J_0_k4])
#dwell_sd_J_0_k = np.array([dwell_sd_J_0_k0,dwell_sd_J_0_k1,dwell_sd_J_0_k2,dwell_sd_J_0_k3,dwell_sd_J_0_k4])

dwell_avg_J_0_k = np.array([dwell_avg_J_0_k0,dwell_avg_J_0_k1,dwell_avg_J_0_k2,dwell_avg_J_0_k3,dwell_avg_J_0_k4,
                            dwell_avg_J_0_k5,dwell_avg_J_0_k6,dwell_avg_J_0_k7,dwell_avg_J_0_k8,dwell_avg_J_0_k9,
                            dwell_avg_J_0_k10,dwell_avg_J_0_k11,dwell_avg_J_0_k12,dwell_avg_J_0_k13])
dwell_sd_J_0_k = np.array([dwell_sd_J_0_k0,dwell_sd_J_0_k1,dwell_sd_J_0_k2,dwell_sd_J_0_k3,dwell_sd_J_0_k4,
                           dwell_sd_J_0_k5,dwell_sd_J_0_k6,dwell_sd_J_0_k7,dwell_sd_J_0_k8,dwell_sd_J_0_k9,
                           dwell_sd_J_0_k10,dwell_sd_J_0_k11,dwell_sd_J_0_k12,dwell_sd_J_0_k13])

dwell_avg_J_1_k = np.array([dwell_avg_J_1_k0,dwell_avg_J_1_k1,dwell_avg_J_1_k2,dwell_avg_J_1_k3,dwell_avg_J_1_k4])
dwell_sd_J_1_k = np.array([dwell_sd_J_1_k0,dwell_sd_J_1_k1,dwell_sd_J_1_k2,dwell_sd_J_1_k3,dwell_sd_J_1_k4])

dwell_avg_J_2_k = np.array([dwell_avg_J_2_k0,dwell_avg_J_2_k1,dwell_avg_J_2_k2,dwell_avg_J_2_k3,dwell_avg_J_2_k4])
dwell_sd_J_2_k = np.array([dwell_sd_J_2_k0,dwell_sd_J_2_k1,dwell_sd_J_2_k2,dwell_sd_J_2_k3,dwell_sd_J_2_k4])

dwell_avg_J_3_k = np.array([dwell_avg_J_3_k0,dwell_avg_J_3_k1,dwell_avg_J_3_k2,dwell_avg_J_3_k3,dwell_avg_J_3_k4])
dwell_sd_J_3_k = np.array([dwell_sd_J_3_k0,dwell_sd_J_3_k1,dwell_sd_J_3_k2,dwell_sd_J_3_k3,dwell_sd_J_3_k4])

dwell_avg_J_4_k = np.array([dwell_avg_J_4_k0,dwell_avg_J_4_k1,dwell_avg_J_4_k2,dwell_avg_J_4_k3,dwell_avg_J_4_k4])

dwell_avg_J_5_k = np.array([dwell_avg_J_5_k0,dwell_avg_J_5_k1,dwell_avg_J_5_k2,dwell_avg_J_5_k3,dwell_avg_J_5_k4])

#%%


plt.rcParams["figure.figsize"] = [8.0,6.0]
plt.title('$<\phi>=8/13$',size=17)
plt.xlabel('k',size=14)
plt.ylabel('Dwell time',size=14)


plt.errorbar(k_J0,dwell_avg_J_0_k,dwell_sd_J_0_k,marker='o',linestyle='',capsize=3,label='J = 0 $k_BT$')
#plt.errorbar(k,dwell_avg_J_1_k,dwell_sd_J_1_k,marker='o',linestyle='',capsize=3,label='J = 1 $k_BT$')
#plt.errorbar(k,dwell_avg_J_2_k,dwell_sd_J_2_k,marker='o',linestyle='',capsize=3,label='J = 2 $k_BT$')
#plt.errorbar(k,dwell_avg_J_3_k,dwell_sd_J_3_k,marker='o',linestyle='',capsize=3,label='J = 3 $k_BT$')
#plt.plot(k,dwell_avg_J_4_k,marker='o',linestyle='')
#plt.plot(k,dwell_avg_J_5_k,marker='o',linestyle='')

plt.legend(loc='upper right', prop={'size': 15})


#%%k

u_dwell, co_dwell = np.unique(dwell_J_0,return_counts=True)


plt.bar(u_dwell,co_dwell/len(dwell),width=0.5)
plt.vlines(dwell_avg, 0, 0.01, linestyles='dashed')

avg_line = mlines.Line2D([],[],linestyle='--', color='black', label='Average ={:1.2f}'.format(dwell_avg))

plt.legend(loc='upper right',handles=[avg_line])

#%%


