#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 15:06:05 2021

@author: mariajose
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#%% Useful variables and functions

n=100 #Number of trajectories we want to plot
 #This numbers limits the number of steps we represent
#m=S-l #number used for representation
nmax=12

#Function used for fits
def N(x,tau,nss):
    return nss + (n0-nss)*np.exp(-x/tau) 

#%% Files for resurrection

files1=['/home/mariajose/Escritorio/Simulations/Without depletion/Relaxation times/Traces_Res_Avg_J=0.dat'] 
data1=[]
files2=['/home/mariajose/Escritorio/Simulations/Without depletion/Relaxation times/Traces_Res_Avg_J=1.dat'] 
data2=[]
files3=['/home/mariajose/Escritorio/Simulations/Without depletion/Relaxation times/Traces_Res_Avg_J=2.dat'] 
data3=[]
files4=['/home/mariajose/Escritorio/Simulations/Without depletion/Relaxation times/Traces_Res_Avg_J=3.dat'] 
data4=[]
files5=['/home/mariajose/Escritorio/Simulations/Without depletion/Relaxation times/Traces_Res_Avg_J=4.dat'] 
data5=[]

for data_file in files1:
    data1.append(np.genfromtxt(data_file))
    
for data_file in files2:
    data2.append(np.genfromtxt(data_file))
    
for data_file in files3:
    data3.append(np.genfromtxt(data_file))
    
for data_file in files4:
    data4.append(np.genfromtxt(data_file))
    
for data_file in files5:
    data5.append(np.genfromtxt(data_file))
    
Res_J0 = data1[0]
Res_J1 = data2[0]
Res_J2 = data3[0]
Res_J3 = data4[0]
Res_J4 = data5[0]

#plt.plot(Res_J0[:,0],Res_J0[:,19]*nmax)

#%%

n0=0

t_J0 = np.zeros(n)
ns_J0 = np.zeros(n)
t_err_J0 = np.zeros(n)
ns_err_J0 = np.zeros(n)

t_J1 = np.zeros(n)
ns_J1 = np.zeros(n)
t_err_J1 = np.zeros(n)
ns_err_J1 = np.zeros(n)
    
t_J2 = np.zeros(n)
ns_J2 = np.zeros(n)
t_err_J2 = np.zeros(n)
ns_err_J2 = np.zeros(n)

t_J3 = np.zeros(n)
ns_J3 = np.zeros(n)
t_err_J3 = np.zeros(n)
ns_err_J3 = np.zeros(n)

t_J4 = np.zeros(n)
ns_J4 = np.zeros(n)
t_err_J4 = np.zeros(n)
ns_err_J4 = np.zeros(n)

for i in range(n-1):
    
    
    poRes_J0, pcRes_J0 = curve_fit(N,Res_J0[0:1000,0],Res_J0[0:1000,2*i+1]*nmax)  
    t_J0[i] = poRes_J0[0]
    ns_J0[i] = poRes_J0[1]
    t_err_J0[i] = pcRes_J0[0,0]
    ns_err_J0[i] = pcRes_J0[1,1]
    
    poRes_J1, pcRes_J1 = curve_fit(N,Res_J1[0:2000,0],Res_J1[0:2000,2*i+1]*nmax)  
    t_J1[i] = poRes_J1[0]
    ns_J1[i] = poRes_J1[1]
    t_err_J1[i] = pcRes_J1[0,0]
    ns_err_J1[i] = pcRes_J1[1,1]


    poRes_J2, pcRes_J2 = curve_fit(N,Res_J2[0:2500,0],Res_J2[0:2500,2*i+1]*nmax)  
    t_J2[i] = poRes_J2[0]
    ns_J2[i] = poRes_J2[1]
    t_err_J2[i] = pcRes_J2[0,0]
    ns_err_J2[i] = pcRes_J2[1,1]

    poRes_J3, pcRes_J3 = curve_fit(N,Res_J3[0:5000,0],Res_J3[0:5000,2*i+1]*nmax)  
    t_J3[i] = poRes_J3[0]
    ns_J3[i] = poRes_J3[1]
    t_err_J3[i] = pcRes_J3[0,0]
    ns_err_J3[i] = pcRes_J3[1,1]

    poRes_J4, pcRes_J4 = curve_fit(N,Res_J4[0:9000,0],Res_J4[0:9000,2*i+1]*nmax)  
    t_J4[i] = poRes_J4[0]
    ns_J4[i] = poRes_J4[1]          
    t_err_J4[i] = pcRes_J4[0,0]
    ns_err_J4[i] = pcRes_J4[1,1]
    
        
plt.xlabel('Steady state')
plt.ylabel('Relaxation time')    
plt.errorbar(ns_J0[0:n-1],t_J0[0:n-1],t_err_J0[0:n-1],ns_err_J0[0:n-1],label='J=0',marker='.')
plt.errorbar(ns_J1[0:n-1],t_J1[0:n-1],t_err_J1[0:n-1],ns_err_J1[0:n-1],label='J=1',marker='.')
plt.errorbar(ns_J2[0:n-1],t_J2[0:n-1],t_err_J2[0:n-1],ns_err_J2[0:n-1],label='J=2',marker='.')
#plt.errorbar(ns_J3[0:n-1],t_J3[0:n-1],t_err_J3[0:n-1],ns_err_J3[0:n-1],label='J=3',marker='.')
#plt.errorbar(ns_J4[0:n-1],t_J4[0:n-1],t_err_J4[0:n-1],ns_err_J4[0:n-1],label='J=4',marker='.')

plt.legend()

#np.savetxt('/home/mariajose/Escritorio/Simulations/Without depletion/Relaxation times/ResTime_J=0.dat',(np.c_[ns_J0[0:n-1],t_J0[0:n-1]]))
#np.savetxt('/home/mariajose/Escritorio/Simulations/Without depletion/Relaxation times/ResTime_J=1.dat',(np.c_[ns_J1[0:n-1],t_J1[0:n-1]]))
#np.savetxt('/home/mariajose/Escritorio/Simulations/Without depletion/Relaxation times/ResTime_J=2.dat',(np.c_[ns_J2[0:n-1],t_J2[0:n-1]]))
#np.savetxt('/home/mariajose/Escritorio/Simulations/Without depletion/Relaxation times/ResTime_J=3.dat',(np.c_[ns_J3[0:n-1],t_J3[0:n-1]]))
#np.savetxt('/home/mariajose/Escritorio/Simulations/Without depletion/Relaxation times/ResTime_J=4.dat',(np.c_[ns_J4[0:n-1],t_J4[0:n-1]]))

#%%

files1=['/home/mariajose/Escritorio/Simulations/Without depletion/Relaxation times/Traces_Stall_Avg_J=0.dat'] 
data1=[]
files2=['/home/mariajose/Escritorio/Simulations/Without depletion/Relaxation times/Traces_Stall_Avg_J=1.dat'] 
data2=[]
files3=['/home/mariajose/Escritorio/Simulations/Without depletion/Relaxation times/Traces_Stall_Avg_J=2.dat'] 
data3=[]
files4=['/home/mariajose/Escritorio/Simulations/Without depletion/Relaxation times/Traces_Stall_Avg_J=3.dat'] 
data4=[]
files5=['/home/mariajose/Escritorio/Simulations/Without depletion/Relaxation times/Traces_Stall_Avg_J=3.dat'] 
data5=[]

for data_file in files1:
    data1.append(np.genfromtxt(data_file))
    
for data_file in files2:
    data2.append(np.genfromtxt(data_file))
    
for data_file in files3:
    data3.append(np.genfromtxt(data_file))
    
for data_file in files4:
    data4.append(np.genfromtxt(data_file))
    
for data_file in files5:
    data5.append(np.genfromtxt(data_file))
    
Stall_J0 = data1[0]
Stall_J1 = data2[0]
Stall_J2 = data3[0]
Stall_J3 = data4[0]
Stall_J4 = data5[0]

#plt.plot(Stall_J0[:,0],Stall_J0[:,27]*nmax)


#%%

n0=12

t_J0=np.zeros(n)
ns_J0 = np.zeros(n)
t_err_J0 = np.zeros(n)
ns_err_J0 = np.zeros(n)

t_J1 = np.zeros(n)
ns_J1 = np.zeros(n)
t_err_J1 = np.zeros(n)
ns_err_J1 = np.zeros(n)
    
t_J2 = np.zeros(n)
ns_J2 = np.zeros(n)
t_err_J2 = np.zeros(n)
ns_err_J2 = np.zeros(n)

t_J3 = np.zeros(n)
ns_J3 = np.zeros(n)
t_err_J3 = np.zeros(n)
ns_err_J3 = np.zeros(n)

t_J4 = np.zeros(n)
ns_J4 = np.zeros(n)
t_err_J4 = np.zeros(n)
ns_err_J4 = np.zeros(n)

for i in range(n-1):
    poStall_J0, pcStall_J0 = curve_fit(N,Stall_J0[:,0],Stall_J0[:,2*i+1]*nmax)  
    t_J0[i] = poStall_J0[0]
    ns_J0[i] = poStall_J0[1]
    t_err_J0[i] = pcStall_J0[0,0]
    ns_err_J0[i] = pcStall_J0[1,1]    
    
    poStall_J1, pcStall_J1 = curve_fit(N,Stall_J1[:,0],Stall_J1[:,2*i+1]*nmax)  
    t_J1[i] = poStall_J1[0]
    ns_J1[i] = poStall_J1[1]
    t_err_J1[i] = pcStall_J1[0,0]
    ns_err_J1[i] = pcStall_J1[1,1]

    poStall_J2, pcStall_J2 = curve_fit(N,Stall_J2[:,0],Stall_J2[:,2*i+1]*nmax)  
    t_J2[i] = poStall_J2[0]
    ns_J2[i] = poStall_J2[1] 
    t_err_J2[i] = pcStall_J2[0,0]
    ns_err_J2[i] = pcStall_J2[1,1]

    poStall_J3, pcStall_J3 = curve_fit(N,Stall_J3[:,0],Stall_J3[:,2*i+1]*nmax)  
    t_J3[i] = poStall_J3[0]
    ns_J3[i] = poStall_J3[1]
    t_err_J3[i] = pcStall_J3[0,0]
    ns_err_J3[i] = pcStall_J3[1,1]

    poStall_J4, pcStall_J4 = curve_fit(N,Stall_J4[:,0],Stall_J4[:,2*i+1]*nmax)  
    t_J4[i] = poStall_J4[0]
    ns_J4[i] = poStall_J4[1]         
    
        
plt.errorbar(ns_J0[0:n-1],t_J0[0:n-1],t_err_J0[0:n-1],ns_err_J0[0:n-1],label='J=0',marker='.')
plt.errorbar(ns_J1[0:n-1],t_J1[0:n-1],t_err_J1[0:n-1],ns_err_J1[0:n-1],label='J=1',marker='.')
plt.errorbar(ns_J2[0:n-1],t_J2[0:n-1],t_err_J2[0:n-1],ns_err_J2[0:n-1],label='J=2',marker='.')
plt.errorbar(ns_J3[0:n-1],t_J3[0:n-1],t_err_J3[0:n-1],ns_err_J3[0:n-1],label='J=3',marker='.')
plt.errorbar(ns_J4[0:n-1],t_J4[0:n-1],t_err_J4[0:n-1],ns_err_J4[0:n-1],label='J=',marker='.')

np.savetxt('/home/mariajose/Escritorio/Simulations/Without depletion/Relaxation times/StallTime_J=0.dat',(np.c_[ns_J0[0:n-1],t_J0[0:n-1]]))
np.savetxt('/home/mariajose/Escritorio/Simulations/Without depletion/Relaxation times/StallTime_J=1.dat',(np.c_[ns_J1[0:n-1],t_J1[0:n-1]]))
np.savetxt('/home/mariajose/Escritorio/Simulations/Without depletion/Relaxation times/StallTime_J=2.dat',(np.c_[ns_J2[0:n-1],t_J2[0:n-1]]))
np.savetxt('/home/mariajose/Escritorio/Simulations/Without depletion/Relaxation times/StallTime_J=3.dat',(np.c_[ns_J3[0:n-1],t_J3[0:n-1]]))
np.savetxt('/home/mariajose/Escritorio/Simulations/Without depletion/Relaxation times/StallTime_J=4.dat',(np.c_[ns_J4[0:n-1],t_J4[0:n-1]])) 



