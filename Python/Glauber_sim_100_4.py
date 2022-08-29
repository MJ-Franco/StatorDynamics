#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 16:10:59 2022

@author: mfranco
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
plt.rcParams["figure.figsize"] = [10.0,8.0]



#%% DATA READING

FileR_J0=['/home/procyon/mfranco/Glauber_simulations/Traces_100_Res_J=0.00_N=4.dat']
DataR_J0=[]

FileS_J0=['/home/procyon/mfranco/Glauber_simulations/Traces_100_Stall_J=0.00_N=4.dat']
DataS_J0=[]

FileS_J0_fixed=['/home/procyon/mfranco/Glauber_simulations/Traces_100_Stall_J=0.00_N=4_fixed.dat']
DataS_J0_fixed=[]


for data_file in FileR_J0:
    DataR_J0.append(np.loadtxt(data_file))
    
for data_file in FileS_J0:
    DataS_J0.append(np.loadtxt(data_file))

for data_file in FileS_J0_fixed:
    DataS_J0_fixed.append(np.loadtxt(data_file))
    
Res_J0=DataR_J0[0]
len(Res_J0)

Rel_J0=DataS_J0[0]

Rel_J0_fixed=DataS_J0_fixed[0]




#%%

def N(x,tau,nss):
    return nss + (n0-nss)*np.exp(-x/tau)

nmax=13

nstep_J0=1000

nrep=100

t_J0 = np.zeros((nrep,nstep_J0))

N_Res_J0 = np.zeros((nrep,nstep_J0))
N_Rel_J0 = np.zeros((nrep,nstep_J0))
N_Rel_J0_fixed = np.zeros((nrep,nstep_J0))

tRes_J0 = np.zeros(nrep)
errtRes_J0 = np.zeros(nrep)

tRel_J0= np.zeros(nrep)
errtRel_J0 = np.zeros(nrep)

tRel_J0_fixed=np.zeros(nrep)
errtRel_J0_fixed = np.zeros(nrep)

Dt_J0=np.zeros(nrep)
errDt_J0=np.zeros(nrep)

Dt_J0_fixed=np.zeros(nrep)
errDt_J0_fixed=np.zeros(nrep)


#%% SEPARATE DATA IN ARRAYS OF 100 ROWS (EACH ROW IS ONE SIMULATION)

for i in range(0,nrep):
    for j in range(nstep_J0*i,nstep_J0*(i+1)):
        
        t_J0[i,j-nstep_J0*i] = Res_J0[j,1]
        
        N_Res_J0[i,j-nstep_J0*i] = Res_J0[j,2]
        
        N_Rel_J0[i,j-nstep_J0*i] = Rel_J0[j,2]
        
        N_Rel_J0_fixed[i,j-nstep_J0*i] = Rel_J0_fixed[j,2]
        
    


#%% FIT RESURRECTION

n0=0

for i in range(nrep):
    poRes_J0, pcRes_J0 = curve_fit(N,t_J0[i,:],N_Res_J0[i,:]*nmax)
    
    tRes_J0[i] = poRes_J0[0]
    errRes_J0 = np.sqrt(np.diag(pcRes_J0))
    errtRes_J0[i] = errRes_J0[0]
    

#%% FIT STALL

n0=6

for i in range(nrep):
    
    poRel_J0, pcRel_J0 = curve_fit(N,t_J0[i,:],N_Rel_J0[i,:]*nmax)
    poRel_J0_fixed, pcRel_J0_fixed = curve_fit(N,t_J0[i,:],N_Rel_J0_fixed[i,:]*nmax)
    
    tRel_J0[i] = poRel_J0[0]
    errRel_J0 = np.sqrt(np.diag(pcRel_J0))
    errtRel_J0[i] = errRel_J0[0]
    
    tRel_J0_fixed[i] = poRel_J0_fixed[0]
    errRel_J0_fixed = np.sqrt(np.diag(pcRel_J0_fixed))
    errtRel_J0_fixed[i] = errRel_J0_fixed[0]

    
    
#%% RELAX TIME RELATIVE DIFFERENCE

for i in range(nrep):
    
    Dt_J0[i]=(tRel_J0[i]-tRes_J0[i])/(tRel_J0[i]+tRes_J0[i])
    errDt_J0[i]=2/(tRel_J0[i]+tRes_J0[i])**2*(tRes_J0[i]*errtRel_J0[i] + tRel_J0[i]*errtRes_J0[i])
    
    Dt_J0_fixed[i]=(tRel_J0_fixed[i]-tRes_J0[i])/(tRel_J0_fixed[i]+tRes_J0[i])
    errDt_J0_fixed[i]=2/(tRel_J0_fixed[i]+tRes_J0[i])**2*(tRes_J0[i]*errtRel_J0_fixed[i] + tRel_J0_fixed[i]*errtRes_J0[i])


#%%

J=np.linspace(0,5,6)

t100 = np.linspace(0,100,100)


plt.title(r'J=0, <$\phi$>=4/13, <$n_0$>=6',fontsize=16)
plt.xlabel('# Simulation',fontsize=14)
plt.ylabel(r'$\Delta$',fontsize=14)
plt.grid()
plt.errorbar(t100,Dt_J0,errDt_J0,ls='',marker='o',capsize=3)

plt.show()

plt.title(r'J=0, <$\phi$>=4/13, $n_0$=6',fontsize=16)
plt.xlabel('# Simulation',fontsize=14)
plt.ylabel(r'$\Delta$',fontsize=14)
plt.grid()
plt.ylim(-0.01,0.01)
plt.errorbar(t100,Dt_J0_fixed,errDt_J0_fixed,ls='',marker='^',capsize=3,color='darkorange')

plt.show()


#%%
