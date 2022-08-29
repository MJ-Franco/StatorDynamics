#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 11:38:59 2022

@author: mfranco
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 12:50:18 2022

@author: mfranco
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
plt.rcParams["figure.figsize"] = [10.0,8.0]



#%% DATA READING

FileR_J0=['/home/procyon/mfranco/Glauber_simulations/Traces_100_Res_J=0.00.dat']
DataR_J0=[]

FileS_J0=['/home/procyon/mfranco/Glauber_simulations/Traces_100_Stall_J=0.00.dat']
DataS_J0=[]

FileS_J0_fixed=['/home/procyon/mfranco/Glauber_simulations/Traces_100_Stall_J=0.00_fixed.dat']
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



FileR_J3=['/home/procyon/mfranco/Glauber_simulations/Traces_100_Res_J=3.00.dat']
DataR_J3=[]

FileS_J3=['/home/procyon/mfranco/Glauber_simulations/Traces_100_Stall_J=3.00.dat']
DataS_J3=[]

FileS_J3_fixed=['/home/procyon/mfranco/Glauber_simulations/Traces_100_Stall_J=3.00_fixed.dat']
DataS_J3_fixed=[]


for data_file in FileR_J3:
    DataR_J3.append(np.loadtxt(data_file))
    
for data_file in FileS_J3:
    DataS_J3.append(np.loadtxt(data_file))

for data_file in FileS_J3_fixed:
    DataS_J3_fixed.append(np.loadtxt(data_file))
    
Res_J3=DataR_J3[0]

Rel_J3=DataS_J3[0]

Rel_J3_fixed=DataS_J3_fixed[0]



FileR_J5=['/home/procyon/mfranco/Glauber_simulations/Traces_100_Res_J=5.00.dat']
DataR_J5=[]

FileS_J5=['/home/procyon/mfranco/Glauber_simulations/Traces_100_Stall_J=5.00.dat']
DataS_J5=[]

FileS_J5_fixed=['/home/procyon/mfranco/Glauber_simulations/Traces_100_Stall_J=5.00_fixed.dat']
DataS_J5_fixed=[]

for data_file in FileR_J5:
    DataR_J5.append(np.loadtxt(data_file))
    
for data_file in FileS_J5:
    DataS_J5.append(np.loadtxt(data_file))
    
for data_file in FileS_J5_fixed:
    DataS_J5_fixed.append(np.loadtxt(data_file))

Res_J5=DataR_J5[0]

Rel_J5=DataS_J5[0]

Rel_J5_fixed=DataS_J5_fixed[0]


#%%

def N(x,tau,nss):
    return nss + (n0-nss)*np.exp(-x/tau)

nmax=13

nstep_J0=1000
nstep_J3=5000
nstep_J5=10000
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

t_J3 = np.zeros((nrep,nstep_J3))

N_Res_J3 = np.zeros((nrep,nstep_J3))
N_Rel_J3 = np.zeros((nrep,nstep_J3))
N_Rel_J3_fixed = np.zeros((nrep,nstep_J3))

tRes_J3 = np.zeros(nrep)
errtRes_J3 = np.zeros(nrep)

tRel_J3= np.zeros(nrep)
errtRel_J3 = np.zeros(nrep)

tRel_J3_fixed=np.zeros(nrep)
errtRel_J3_fixed = np.zeros(nrep)

Dt_J3=np.zeros(nrep)
errDt_J3=np.zeros(nrep)

Dt_J3_fixed=np.zeros(nrep)
errDt_J3_fixed=np.zeros(nrep)


t_J5 = np.zeros((nrep,nstep_J5))


N_Res_J5 = np.zeros((nrep,nstep_J5))
N_Rel_J5 = np.zeros((nrep,nstep_J5))
N_Rel_J5_fixed = np.zeros((nrep,nstep_J5))

tRes_J5 = np.zeros(nrep)
errtRes_J5 = np.zeros(nrep)

tRel_J5= np.zeros(nrep)
errtRel_J5 = np.zeros(nrep)

Dt_J5=np.zeros(nrep)
errDt_J5=np.zeros(nrep)

tRel_J5_fixed= np.zeros(nrep)
errtRel_J5_fixed = np.zeros(nrep)

Dt_J5_fixed=np.zeros(nrep)
errDt_J5_fixed=np.zeros(nrep)

#%% SEPARATE DATA IN ARRAYS OF 100 ROWS (EACH ROW IS ONE SIMULATION)

for i in range(0,nrep):
    for j in range(nstep_J0*i,nstep_J0*(i+1)):
        
        t_J0[i,j-nstep_J0*i] = Res_J0[j,1]
        
        N_Res_J0[i,j-nstep_J0*i] = Res_J0[j,2]
        
        N_Rel_J0[i,j-nstep_J0*i] = Rel_J0[j,2]
        
        N_Rel_J0_fixed[i,j-nstep_J0*i] = Rel_J0_fixed[j,2]
        
    
for i in range(0,nrep):
    for j in range(nstep_J3*i,nstep_J3*(i+1)):
        
        t_J3[i,j-nstep_J3*i] = Res_J3[j,1]
        
        N_Res_J3[i,j-nstep_J3*i] = Res_J3[j,2]
        
        N_Rel_J3[i,j-nstep_J3*i] = Rel_J3[j,2]
        
        N_Rel_J3_fixed[i,j-nstep_J3*i] = Rel_J3_fixed[j,2]
        
        
for i in range(0,nrep):
    for j in range(nstep_J5*i,nstep_J5*(i+1)):
        
        t_J5[i,j-nstep_J5*i] = Res_J5[j,1]
        
        N_Res_J5[i,j-nstep_J5*i] = Res_J5[j,2]
        
        N_Rel_J5[i,j-nstep_J5*i] = Rel_J5[j,2]
        
        N_Rel_J5_fixed[i,j-nstep_J5*i] = Rel_J5_fixed[j,2]

#%% FIT RESURRECTION

n0=0

for i in range(nrep):
    poRes_J0, pcRes_J0 = curve_fit(N,t_J0[i,:],N_Res_J0[i,:]*nmax)
    
    tRes_J0[i] = poRes_J0[0]
    errRes_J0 = np.sqrt(np.diag(pcRes_J0))
    errtRes_J0[i] = errRes_J0[0]
    
    poRes_J3, pcRes_J3 = curve_fit(N,t_J3[i,:],N_Res_J3[i,:]*nmax)
    
    tRes_J3[i] = poRes_J3[0]
    errRes_J3 = np.sqrt(np.diag(pcRes_J3))
    errtRes_J3[i] = errRes_J3[0]
    
    poRes_J5, pcRes_J5 = curve_fit(N,t_J5[i,:],N_Res_J5[i,:]*nmax)
    
    tRes_J5[i] = poRes_J5[0]
    errRes_J5 = np.sqrt(np.diag(pcRes_J5))
    errtRes_J5[i] = errRes_J5[0]
    
#%% FIT STALL

n0=9

for i in range(nrep):
    
    poRel_J0, pcRel_J0 = curve_fit(N,t_J0[i,:],N_Rel_J0[i,:]*nmax)
    poRel_J0_fixed, pcRel_J0_fixed = curve_fit(N,t_J0[i,:],N_Rel_J0_fixed[i,:]*nmax)
    
    tRel_J0[i] = poRel_J0[0]
    errRel_J0 = np.sqrt(np.diag(pcRel_J0))
    errtRel_J0[i] = errRel_J0[0]
    
    tRel_J0_fixed[i] = poRel_J0_fixed[0]
    errRel_J0_fixed = np.sqrt(np.diag(pcRel_J0_fixed))
    errtRel_J0_fixed[i] = errRel_J0_fixed[0]
    
    
    poRel_J3, pcRel_J3 = curve_fit(N,t_J3[i,:],N_Rel_J3[i,:]*nmax)
    poRel_J3_fixed, pcRel_J3_fixed = curve_fit(N,t_J3[i,:],N_Rel_J3_fixed[i,:]*nmax)
    
    tRel_J3[i] = poRel_J3[0]
    errRel_J3 = np.sqrt(np.diag(pcRel_J3))
    errtRel_J3[i] = errRel_J3[0]
    
    tRel_J3_fixed[i] = poRel_J3_fixed[0]
    errRel_J3_fixed = np.sqrt(np.diag(pcRel_J3_fixed))
    errtRel_J3_fixed[i] = errRel_J3_fixed[0]
    
    
    poRel_J5, pcRel_J5 = curve_fit(N,t_J5[i,:],N_Rel_J5[i,:]*nmax)
    poRel_J5_fixed, pcRel_J5_fixed = curve_fit(N,t_J5[i,:],N_Rel_J5_fixed[i,:]*nmax)
    
    tRel_J5[i] = poRel_J5[0]
    errRel_J5 = np.sqrt(np.diag(pcRel_J5))
    errtRel_J5[i] = errRel_J5[0]
    
    tRel_J5_fixed[i] = poRel_J5_fixed[0]
    errRel_J5_fixed = np.sqrt(np.diag(pcRel_J5_fixed))
    errtRel_J5_fixed[i] = errRel_J5_fixed[0]
    
    
#%% RELAX TIME RELATIVE DIFFERENCE

for i in range(nrep):
    
    Dt_J0[i]=(tRel_J0[i]-tRes_J0[i])/(tRel_J0[i]+tRes_J0[i])
    errDt_J0[i]=2/(tRel_J0[i]+tRes_J0[i])**2*(tRes_J0[i]*errtRel_J0[i] + tRel_J0[i]*errtRes_J0[i])
    
    Dt_J0_fixed[i]=(tRel_J0_fixed[i]-tRes_J0[i])/(tRel_J0_fixed[i]+tRes_J0[i])
    errDt_J0_fixed[i]=2/(tRel_J0_fixed[i]+tRes_J0[i])**2*(tRes_J0[i]*errtRel_J0_fixed[i] + tRel_J0_fixed[i]*errtRes_J0[i])
    
    
    Dt_J3[i]=(tRel_J3[i]-tRes_J3[i])/(tRel_J3[i]+tRes_J3[i])
    errDt_J3[i]=2/(tRel_J3[i]+tRes_J3[i])**2*(tRes_J3[i]*errtRel_J3[i] + tRel_J3[i]*errtRes_J3[i])
    
    Dt_J3_fixed[i]=(tRel_J3_fixed[i]-tRes_J3[i])/(tRel_J3_fixed[i]+tRes_J3[i])
    errDt_J3_fixed[i]=2/(tRel_J3_fixed[i]+tRes_J3[i])**2*(tRes_J3[i]*errtRel_J3_fixed[i] + tRel_J3_fixed[i]*errtRes_J3[i])
    
    
    Dt_J5[i]=(tRel_J5[i]-tRes_J5[i])/(tRel_J5[i]+tRes_J5[i])
    errDt_J5[i]=2/(tRel_J5[i]+tRes_J5[i])**2*(tRes_J5[i]*errtRel_J5[i] + tRel_J5[i]*errtRes_J5[i])
    
    Dt_J5_fixed[i]=(tRel_J5_fixed[i]-tRes_J5[i])/(tRel_J5_fixed[i]+tRes_J5[i])
    errDt_J5_fixed[i]=2/(tRel_J5_fixed[i]+tRes_J5[i])**2*(tRes_J5[i]*errtRel_J5_fixed[i] + tRel_J5_fixed[i]*errtRes_J5[i])

#%%

J=np.linspace(0,5,6)

t100 = np.linspace(0,100,100)


plt.title(r'J=0, <$\phi$>=1/2, <$n_0$>=9',fontsize=16)
plt.xlabel('# Simulation',fontsize=14)
plt.ylabel(r'$\Delta$',fontsize=14)
plt.grid()
plt.errorbar(t100,Dt_J0,errDt_J0,ls='',marker='o',capsize=3)

plt.show()

plt.title(r'J=0, <$\phi$>=1/2, $n_0$=9',fontsize=16)
plt.xlabel('# Simulation',fontsize=14)
plt.ylabel(r'$\Delta$',fontsize=14)
plt.grid()
plt.ylim(0,0.015)
plt.errorbar(t100,Dt_J0_fixed,errDt_J0_fixed,ls='',marker='^',capsize=3,color='darkorange')

plt.show()

plt.title(r'J=3, <$\phi$>=1/2, <$n_0$>=9',fontsize=16)
plt.xlabel('# Simulation',fontsize=14)
plt.ylabel(r'$\Delta$',fontsize=14)
plt.grid()
plt.errorbar(t100,Dt_J3,errDt_J3,ls='',marker='o',capsize=3)

plt.show()


plt.title(r'J=3, <$\phi$>=1/2, $n_0$=9',fontsize=16)
plt.xlabel('# Simulation',fontsize=14)
plt.ylabel(r'$\Delta$',fontsize=14)
plt.grid()
plt.ylim(-0.01,-0.0050)
plt.errorbar(t100,Dt_J3_fixed,errDt_J3_fixed,ls='',marker='^',capsize=3,color='darkorange')

plt.show()


plt.title(r'J=5, <$\phi$>=1/2, <$n_0$>=9',fontsize=16)
plt.xlabel('# Simulation',fontsize=14)
plt.ylabel(r'$\Delta$',fontsize=14)
plt.grid()
plt.errorbar(t100,Dt_J5,errDt_J5,ls='',marker='o',capsize=3)


plt.show()

plt.title(r'J=5, <$\phi$>=1/2, $n_0$=9',fontsize=16)
plt.xlabel('# Simulation',fontsize=14)
plt.ylabel(r'$\Delta$',fontsize=14)
plt.grid()
plt.errorbar(t100,Dt_J5_fixed,errDt_J5_fixed,ls='',marker='o',capsize=3)


plt.show()

#%%

t0=np.linspace(0,0,100)
t3=np.linspace(3,3,100)
t5=np.linspace(5,5,100)


plt.title('100 simulations',fontsize=16)
plt.ylabel(r'$\Delta$',fontsize=14)
plt.xlabel('J',fontsize=14)
plt.grid()
plt.errorbar(t0,Dt_J0,errDt_J0,ls='',marker='o',capsize=3,color='hotpink',label='<$n_0$>=9')
plt.errorbar(t0,Dt_J0_fixed,errDt_J0_fixed,ls='',marker='^',capsize=3,color='dodgerblue',label='$n_0$=9')
plt.errorbar(t3,Dt_J3,errDt_J3,ls='',marker='o',capsize=3,color='hotpink')
plt.errorbar(t3,Dt_J3_fixed,errDt_J3_fixed,ls='',marker='^',capsize=3,color='dodgerblue')
plt.errorbar(t5,Dt_J5,errDt_J5,ls='',marker='o',capsize=3,color='hotpink')
plt.errorbar(t5,Dt_J5_fixed,errDt_J5_fixed,ls='',marker='^',capsize=3,color='dodgerblue')

plt.legend(prop={'size': 12})

#%%

avDt=np.zeros(3)
erravDt=np.zeros(3)

#%%

avDt_J0 = sum(Dt_J0)/nrep
avDt_J3 = sum(Dt_J3)/nrep
avDt_J5 = sum(Dt_J5)/nrep

erravDt_J0=np.sqrt(sum((Dt_J0-avDt_J0)**2)/nrep)
erravDt_J3=np.sqrt(sum((Dt_J3-avDt_J3)**2)/nrep)
erravDt_J5=np.sqrt(sum((Dt_J3-avDt_J3)**2)/nrep)

avDt[0] = avDt_J0
avDt[1] = avDt_J3
avDt[2] = avDt_J5

erravDt[0] = erravDt_J0
erravDt[1] = erravDt_J3
erravDt[2] = erravDt_J5

J=(0,3,5)


plt.title('Average 100 simulations',fontsize=16)
plt.ylabel(r'$\Delta$',fontsize=14)
plt.xlabel('J',fontsize=14)
plt.grid()
plt.errorbar(J,avDt,erravDt, ls='',marker='o',capsize=3,color='hotpink',label='<$n_0$>=9')
plt.errorbar(t0,Dt_J0_fixed,errDt_J0_fixed,ls='',marker='^',capsize=3,color='dodgerblue',label='$n_0$=9')
plt.errorbar(t3,Dt_J3_fixed,errDt_J3_fixed,ls='',marker='^',capsize=3,color='dodgerblue')
plt.errorbar(t5,Dt_J5_fixed,errDt_J5_fixed,ls='',marker='^',capsize=3,color='dodgerblue')

plt.legend(prop={'size': 12})