#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 12:08:10 2021

@author: mariajose
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.lines as mlines


#%% Useful variables and functions

T=10 #Number of trajectories we want to plot
S=7000 #Total number of steps of the simulation
l=20000 #This numbers limits the number of steps we represen
m=S-l #number used for representation
nmax=13

#Function used for fits
def N(x,tau,nss):
    return nss + (n0-nss)*np.exp(-x/tau) 


#%% Files

files1=['/media/luke/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg.dat']
data1=[]

files2=['/media/luke/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg.dat']
data2=[]

for data_file in files1:
    data1.append(np.loadtxt(data_file))

for data_file in files2:
    data2.append(np.loadtxt(data_file))

Res=data1[0]
Stall=data2[0]


Res[0:l,1]


#%%

n0=0

#p0 = [10000,2.5]

poRes, pcRes = curve_fit(N,Res[0:l,0],Res[0:l,1])
modelRes=N(Res[0:l,0],*poRes) 

#%%

n0=6.25

#p0 = [200,2.5]

poStall, pcStall = curve_fit(N,Stall[0:l,0],Stall[0:l,1])
modelStall=N(Stall[0:l,0],*poStall)


#%%

plt.title(r'300 nm',fontsize=14)
plt.xlabel('Time',fontsize=16)
plt.ylabel('Stator number',fontsize=16)
plt.rcParams["figure.figsize"] = [8.0,6.0]
plt.grid()

plt.plot(Res[0:l,0],Res[0:l,1],label=r'$\tau$={:1.2f}'.format(poRes[0]))
plt.plot(Stall[0:l,0],Stall[0:l,1],label=r'$\tau$={:1.2f}'.format(poStall[0]))

plt.legend()

#%%For Langmuir

def Nl(t,kon,koff,n0):
    nss = nmax/(1. + koff/kon)
    return nss + (n0-nss)*np.exp(-t*(kon+koff)) 

t = np.linspace(0,1800,10000)

kon_rel = 2.18e-3
koff_rel = 3.96e-3

kon_res = 0.0189
koff_res = 0.0451

plt.title(r'300 nm',fontsize=14)
plt.xlabel('Time',fontsize=16)
plt.ylabel('Stator number',fontsize=16)
plt.rcParams["figure.figsize"] = [8.0,6.0]
plt.grid()

plt.plot(t, Nl(t,kon_rel,koff_rel,0.),color='blue')
plt.plot(t, Nl(t,kon_rel,koff_rel,6.2),color='orange')
plt.plot(Res[0:l,0],Res[0:l,1],label=r'$\tau$={:1.2f}'.format(poRes[0]),color='blue',linestyle='--')
plt.plot(Stall[0:l,0],Stall[0:l,1],label=r'$\tau$={:1.2f}'.format(poStall[0]),color='orange',linestyle='--')

black_line1 = mlines.Line2D([],[],linestyle='-', label='Langmuir')
black_line2 = mlines.Line2D([],[],linestyle='--', label='Simulation')

plt.legend(handles=[black_line1,black_line2])


#%%For cooperativity

files3=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Traces_Res_J=4.00  _mu=-4.00  .dat']
data3=[]

files4=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Traces_Stall_J=4.00  _mu=-4.00  _.dat']
data4=[]

for data_file in files3:
    data3.append(np.loadtxt(data_file))

for data_file in files4:
    data4.append(np.loadtxt(data_file))


Res_Glauber=data3[0]
Stall_Glauber=data4[0]

plt.title(r'J=4, $\mu=-4$',fontsize=14)
plt.xlabel('Time',fontsize=16)
plt.ylabel('Stator number',fontsize=16)
plt.rcParams["figure.figsize"] = [8.0,6.0]
plt.grid()

m=2700

plt.plot(Res_Glauber[0:m,0]/13,Res_Glauber[0:m,1]*nmax,label='Glauber',linestyle='--',color='black')
plt.plot(Stall_Glauber[0:m,0]/13,Stall_Glauber[0:m,1]*nmax,linestyle='--',color='black')
plt.plot(Res[0:l,0],Res[0:l,1],label='KMC',color='red',alpha=0.5)
plt.plot(Stall[0:l,0],Stall[0:l,1],color='red',alpha=0.5)

plt.legend()
