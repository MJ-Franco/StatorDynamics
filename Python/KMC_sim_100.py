#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 12:12:36 2022

@author: mfranco
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
plt.rcParams["figure.figsize"] = [10.0,8.0]



#%% DATA READING

FileRHF_J0=['/home/procyon/mfranco/KMC_simulations/Traces_100_Res_J=0.00_N=6.50.dat']
DataRHF_J0=[]

FileSHF_J0=['/home/procyon/mfranco/KMC_simulations/Traces_100_Stall_J=0.00_N=6.50.dat']
DataSHF_J0=[]

FileSHF_J0_fixed=['/home/procyon/mfranco/KMC_simulations/Traces_100_Stall_J=0.00_N=6.50_fixed.dat']
DataSHF_J0_fixed=[]


for data_file in FileRHF_J0:
    DataRHF_J0.append(np.loadtxt(data_file))
    
for data_file in FileSHF_J0:
    DataSHF_J0.append(np.loadtxt(data_file))

for data_file in FileSHF_J0_fixed:
    DataSHF_J0_fixed.append(np.loadtxt(data_file))
    
ResHF_J0=DataRHF_J0[0]

RelHF_J0=DataSHF_J0[0]

RelHF_J0_fixed=DataSHF_J0_fixed[0]



FileR4_J0=['/home/procyon/mfranco/KMC_simulations/Traces_100_Res_J=0.00_N=4.00.dat']
DataR4_J0=[]

FileS4_J0=['/home/procyon/mfranco/KMC_simulations/Traces_100_Stall_J=0.00_N=4.00.dat']
DataS4_J0=[]

FileS4_J0_fixed=['/home/procyon/mfranco/KMC_simulations/Traces_100_Stall_J=0.00_N=4.00_fixed.dat']
DataS4_J0_fixed=[]


for data_file in FileR4_J0:
    DataR4_J0.append(np.loadtxt(data_file))
    
for data_file in FileS4_J0:
    DataS4_J0.append(np.loadtxt(data_file))

for data_file in FileS4_J0_fixed:
    DataS4_J0_fixed.append(np.loadtxt(data_file))
    
Res4_J0=DataR4_J0[0]

Rel4_J0=DataS4_J0[0]

Rel4_J0_fixed=DataS4_J0_fixed[0]



FileR8_J0=['/home/procyon/mfranco/KMC_simulations/Traces_100_Res_J=0.00_N=8.00.dat']
DataR8_J0=[]

FileS8_J0=['/home/procyon/mfranco/KMC_simulations/Traces_100_Stall_J=0.00_N=8.00.dat']
DataS8_J0=[]

FileS8_J0_fixed=['/home/procyon/mfranco/KMC_simulations/Traces_100_Stall_J=0.00_N=8.00_fixed.dat']
DataS8_J0_fixed=[]


for data_file in FileR8_J0:
    DataR8_J0.append(np.loadtxt(data_file))
    
for data_file in FileS8_J0:
    DataS8_J0.append(np.loadtxt(data_file))

for data_file in FileS8_J0_fixed:
    DataS8_J0_fixed.append(np.loadtxt(data_file))
    
Res8_J0=DataR8_J0[0]

Rel8_J0=DataS8_J0[0]

Rel8_J0_fixed=DataS8_J0_fixed[0]

#%%

nmax=13

def N(x,tau,nss):
    return nss + (n0-nss)*np.exp(-x/tau)


nrep=100

tResHF_J0=np.zeros(nrep)
errtResHF_J0=np.zeros(nrep)
tRelHF_J0=np.zeros(nrep)
errtRelHF_J0=np.zeros(nrep)
tRelHF_J0_fixed=np.zeros(nrep)
errtRelHF_J0_fixed=np.zeros(nrep)
DtHF_J0=np.zeros(nrep)
errDtHF_J0=np.zeros(nrep)
DtHF_J0_fixed=np.zeros(nrep)
errDtHF_J0_fixed=np.zeros(nrep)

tRes4_J0=np.zeros(nrep)
errtRes4_J0=np.zeros(nrep)
tRel4_J0=np.zeros(nrep)
errtRel4_J0=np.zeros(nrep)
tRel4_J0_fixed=np.zeros(nrep)
errtRel4_J0_fixed=np.zeros(nrep)
Dt4_J0=np.zeros(nrep)
errDt4_J0=np.zeros(nrep)
Dt4_J0_fixed=np.zeros(nrep)
errDt4_J0_fixed=np.zeros(nrep)

tRes8_J0=np.zeros(nrep)
errtRes8_J0=np.zeros(nrep)
tRel8_J0=np.zeros(nrep)
errtRel8_J0=np.zeros(nrep)
tRel8_J0_fixed=np.zeros(nrep)
errtRel8_J0_fixed=np.zeros(nrep)
Dt8_J0=np.zeros(nrep)
errDt8_J0=np.zeros(nrep)
Dt8_J0_fixed=np.zeros(nrep)
errDt8_J0_fixed=np.zeros(nrep)

#%%

#ARRAYS FOR HALF-FILLING
        
j=0
len_ResHF_J0=0
while(ResHF_J0[j,0]==1.):
        
    len_ResHF_J0 = len_ResHF_J0 + 1
        
    j = j + 1 


tHF_J0 = np.zeros((nrep,len_ResHF_J0))
N_ResHF_J0 = np.zeros((nrep,len_ResHF_J0))
N_RelHF_J0 = np.zeros((nrep,len_ResHF_J0))
N_RelHF_J0_fixed = np.zeros((nrep,len_ResHF_J0))


# ARRAYS FOR N=4

j=0
len_Res4_J0=0
while(Res4_J0[j,0]==1.):
        
    len_Res4_J0 = len_Res4_J0 + 1
        
    j = j + 1 
    
t4Res_J0 = np.zeros((nrep,len_Res4_J0))
N_Res4_J0 = np.zeros((nrep,len_Res4_J0))

j=0
len_Rel4_J0=0
while(Rel4_J0[j,0]==1.):
        
    len_Rel4_J0 = len_Rel4_J0 + 1
        
    j = j + 1 

t4Rel_J0 = np.zeros((nrep,len_Rel4_J0))
N_Rel4_J0 = np.zeros((nrep,len_Rel4_J0))
N_Rel4_J0_fixed = np.zeros((nrep,len_Rel4_J0))


# ARRAYS FOR N=8

j=0
len_Res8_J0=0
while(Res8_J0[j,0]==1.):
        
    len_Res8_J0 = len_Res8_J0 + 1
        
    j = j + 1 
    
t8Res_J0 = np.zeros((nrep,len_Res8_J0))
N_Res8_J0 = np.zeros((nrep,len_Res8_J0))

j=0
len_Rel8_J0=0
while(Rel8_J0[j,0]==1.):
        
    len_Rel8_J0 = len_Rel8_J0 + 1
        
    j = j + 1 


t8Rel_J0 = np.zeros((nrep,len_Rel8_J0))
N_Rel8_J0 = np.zeros((nrep,len_Rel8_J0))

j=0
len_Rel8_J0_fixed=0
while(Rel8_J0_fixed[j,0]==1.):
        
    len_Rel8_J0_fixed = len_Rel8_J0_fixed + 1
        
    j = j + 1 

t8Rel_J0_fixed = np.zeros((nrep,len_Rel8_J0_fixed))
N_Rel8_J0_fixed = np.zeros((nrep,len_Rel8_J0_fixed))




#%%

# HALF-FILLING 

j=0        
for i in range(nrep):
    
    flag = float(i+1)
    
    kResHF_J0=0
    while(ResHF_J0[j,0]<flag):
        
        tHF_J0[i,kResHF_J0] = ResHF_J0[j,1]
        
        N_ResHF_J0[i,kResHF_J0] = ResHF_J0[j,2]
        
        kResHF_J0 = kResHF_J0 + 1
        
        j = j + 1
  

j=0        
for i in range(nrep):
    
    flag = float(i+1)
      
    kRelHF_J0=0
    while(RelHF_J0[j,0]<flag):
        
        N_RelHF_J0[i,kRelHF_J0] = RelHF_J0[j,2]
        
        N_RelHF_J0_fixed[i,kRelHF_J0] = RelHF_J0_fixed[j,2]
        
        kRelHF_J0 = kRelHF_J0 + 1
        
        j = j + 1


# N=4

j=0        
for i in range(nrep):
    
    flag = float(i+1)

    
    kRes4_J0=0
    while(Res4_J0[j,0]<flag):
        
        if i==0:
            print(i)
        
        t4Res_J0[i,kRes4_J0] = Res4_J0[j,1]
        
        N_Res4_J0[i,kRes4_J0] = Res4_J0[j,2]
        
        kRes4_J0 = kRes4_J0 + 1
        
        j = j + 1

j=0        
for i in range(nrep):
    
    flag = float(i+1)
       
    kRel4_J0=0
    while(Rel4_J0[j,0]<flag):
        t4Rel_J0[i,kRel4_J0] = Rel4_J0[j,1]
        
        N_Rel4_J0[i,kRel4_J0] = Rel4_J0[j,2]
        
        N_Rel4_J0_fixed[i,kRel4_J0] = Rel4_J0_fixed[j,2]
        
        kRel4_J0 = kRel4_J0 + 1
        
        j = j + 1
        

# N=8 
        
j=0        
for i in range(nrep):
    
    flag = float(i+1)

    
    kRes8_J0=0
    while(Res8_J0[j,0]<flag):
        
        if i==0:
            print(flag)
        
        t8Res_J0[i,kRes8_J0] = Res8_J0[j,1]
        
        N_Res8_J0[i,kRes8_J0] = Res8_J0[j,2]
        
        kRes8_J0 = kRes8_J0 + 1
        
        j = j + 1


        
for i in range(nrep):
    
    flag = float(i+1)
    
    kRel8_J0=0
    for k in range(len(Rel8_J0[:,0])):
        
        if ((Rel8_J0[k,0]==flag) and (kRel8_J0<len_Rel8_J0)):
            
            t8Rel_J0[i,kRel8_J0] = Rel8_J0[k,1]
            
            N_Rel8_J0[i,kRel8_J0] = Rel8_J0[k,2]
            
            kRel8_J0 = kRel8_J0 + 1
            
        if Rel8_J0[k,0]>flag:
            
            break

    kRel8_J0_fixed=0
    for k in range(len(Rel8_J0_fixed[:,0])):
        
        if ((Rel8_J0_fixed[k,0]==flag) and (kRel8_J0_fixed<len_Rel8_J0_fixed)):
            
            t8Rel_J0_fixed[i,kRel8_J0_fixed] = Rel8_J0_fixed[k,1]
        
            N_Rel8_J0_fixed[i,kRel8_J0_fixed] = Rel8_J0_fixed[k,2]
            
            kRel8_J0_fixed = kRel8_J0_fixed + 1
            
        if Rel8_J0_fixed[k,0]>flag:
             
            break
       
      
#%% FIT RESURRECTION

n0=0.0

for i in range(nrep):
    poResHF_J0, pcResHF_J0 = curve_fit(N,tHF_J0[i,:],N_ResHF_J0[i,:])
    
    tResHF_J0[i] = poResHF_J0[0]
    errResHF_J0 = np.sqrt(np.diag(pcResHF_J0))
    errtResHF_J0[i] = errResHF_J0[0]
    
    
    poRes4_J0, pcRes4_J0 = curve_fit(N,t4Res_J0[i,:],N_Res4_J0[i,:])
    
    tRes4_J0[i] = poRes4_J0[0]
    errRes4_J0 = np.sqrt(np.diag(pcRes4_J0))
    errtRes4_J0[i] = errRes4_J0[0]
    
    
    poRes8_J0, pcRes8_J0 = curve_fit(N,t8Res_J0[i,:],N_Res8_J0[i,:])
    
    tRes8_J0[i] = poRes8_J0[0]
    errRes8_J0 = np.sqrt(np.diag(pcRes8_J0))
    errtRes8_J0[i] = errRes8_J0[0]




#%%

n0=9.0

for i in range(nrep):
    
    poRelHF_J0, pcRelHF_J0 = curve_fit(N,tHF_J0[i,:],N_RelHF_J0[i,:])
    poRelHF_J0_fixed, pcRelHF_J0_fixed = curve_fit(N,tHF_J0[i,:],N_RelHF_J0_fixed[i,:])
    
    tRelHF_J0[i] = poRelHF_J0[0]
    errRelHF_J0 = np.sqrt(np.diag(pcRelHF_J0))
    errtRelHF_J0[i] = errRelHF_J0[0]
    
    tRelHF_J0_fixed[i] = poRelHF_J0_fixed[0]
    errRelHF_J0_fixed = np.sqrt(np.diag(pcRelHF_J0_fixed))
    errtRelHF_J0_fixed[i] = errRelHF_J0_fixed[0]
    
    
    poRel8_J0, pcRel8_J0 = curve_fit(N,t8Rel_J0[i,:],N_Rel8_J0[i,:])
    poRel8_J0_fixed, pcRel8_J0_fixed = curve_fit(N,t8Rel_J0_fixed[i,:],N_Rel8_J0_fixed[i,:])
    
    tRel8_J0[i] = poRel8_J0[0]
    errRel8_J0 = np.sqrt(np.diag(pcRel8_J0))
    errtRel8_J0[i] = errRel8_J0[0]
    
    tRel8_J0_fixed[i] = poRel8_J0_fixed[0]
    errRel8_J0_fixed = np.sqrt(np.diag(pcRel8_J0_fixed))
    errtRel8_J0_fixed[i] = errRel8_J0_fixed[0]

#%%

n0=6.0

for i in range(nrep):
    
    poRel4_J0, pcRel4_J0 = curve_fit(N,t4Rel_J0[i,:],N_Rel4_J0[i,:])
    poRel4_J0_fixed, pcRel4_J0_fixed = curve_fit(N,t4Rel_J0[i,:],N_Rel4_J0_fixed[i,:])
    
    tRel4_J0[i] = poRel4_J0[0]
    errRel4_J0 = np.sqrt(np.diag(pcRel4_J0))
    errtRel4_J0[i] = errRel4_J0[0]
    
    tRel4_J0_fixed[i] = poRel4_J0_fixed[0]
    errRel4_J0_fixed = np.sqrt(np.diag(pcRel4_J0_fixed))
    errtRel4_J0_fixed[i] = errRel4_J0_fixed[0]

#%%

for i in range(nrep):
    
    DtHF_J0[i]=(tRelHF_J0[i]-tResHF_J0[i])/(tRelHF_J0[i]+tResHF_J0[i])
    errDtHF_J0[i]=2./(tRelHF_J0[i]+tResHF_J0[i])**2*(tResHF_J0[i]*errtRelHF_J0[i] + tRelHF_J0[i]*errtResHF_J0[i])
    
    DtHF_J0_fixed[i]=(tRelHF_J0_fixed[i]-tResHF_J0[i])/(tRelHF_J0_fixed[i]+tResHF_J0[i])
    errDtHF_J0_fixed[i]=2./(tRelHF_J0_fixed[i]+tResHF_J0[i])**2*(tResHF_J0[i]*errtRelHF_J0_fixed[i] + tRelHF_J0_fixed[i]*errtResHF_J0[i])
    
    
    Dt4_J0[i]=(tRel4_J0[i]-tRes4_J0[i])/(tRel4_J0[i]+tRes4_J0[i])
    errDt4_J0[i]=2./(tRel4_J0[i]+tRes4_J0[i])**2*(tRes4_J0[i]*errtRel4_J0[i] + tRel4_J0[i]*errtRes4_J0[i])
    
    Dt4_J0_fixed[i]=(tRel4_J0_fixed[i]-tRes4_J0[i])/(tRel4_J0_fixed[i]+tRes4_J0[i])
    errDt4_J0_fixed[i]=2./(tRel4_J0_fixed[i]+tRes4_J0[i])**2*(tRes4_J0[i]*errtRel4_J0_fixed[i] + tRel4_J0_fixed[i]*errtRes4_J0[i])


    Dt8_J0[i]=(tRel8_J0[i]-tRes8_J0[i])/(tRel8_J0[i]+tRes8_J0[i])
    errDt8_J0[i]=2./(tRel8_J0[i]+tRes8_J0[i])**2*(tRes8_J0[i]*errtRel8_J0[i] + tRel8_J0[i]*errtRes8_J0[i])
    
    Dt8_J0_fixed[i]=(tRel8_J0_fixed[i]-tRes8_J0[i])/(tRel8_J0_fixed[i]+tRes8_J0[i])
    errDt8_J0_fixed[i]=2./(tRel8_J0_fixed[i]+tRes8_J0[i])**2*(tRes8_J0[i]*errtRel8_J0_fixed[i] + tRel8_J0_fixed[i]*errtRes8_J0[i])


#%%

plt.plot(t8Res_J0[80,:],N_Res8_J0[80,:])
plt.plot(t8Rel_J0[1,:],N_Rel8_J0[1,:])
plt.show()

#%%

J=np.linspace(0,5,6)

t100 = np.linspace(0,100,100)


plt.title(r'J=0, <$\phi$>=1/2, <$n_0$>=$n_0$=9',fontsize=16)
plt.xlabel('# KMC simulation',fontsize=14)
plt.ylabel(r'$\Delta$',fontsize=14)
plt.grid()
plt.errorbar(t100,DtHF_J0,errDtHF_J0,ls='',marker='o',capsize=3,label='IC Distribution')
plt.errorbar(t100,DtHF_J0_fixed,errDtHF_J0_fixed,ls='',marker='^',capsize=3, label='IC Fixed')

plt.legend(prop={'size': 12})

plt.show()


plt.title(r'J=0, <$\phi$>=4/13, <$n_0$>=$n_0$=6',fontsize=16)
plt.xlabel('# KMC simulation',fontsize=14)
plt.ylabel(r'$\Delta$',fontsize=14)
plt.grid()
plt.errorbar(t100,Dt4_J0,errDt4_J0,ls='',marker='o',capsize=3,label='IC Distribution')
plt.errorbar(t100,Dt4_J0_fixed,errDt4_J0_fixed,ls='',marker='^',capsize=3,label='IC Fixed')

plt.legend(prop={'size': 12})

plt.show()


plt.title(r'J=0, <$\phi$>=8/13, <$n_0$>=$n_0$=9',fontsize=16)
plt.xlabel('# KMC Simulation',fontsize=14)
plt.ylabel(r'$\Delta$',fontsize=14)
plt.grid()
plt.errorbar(t100[1:100],Dt8_J0[1:100],errDt8_J0[1:100],ls='',marker='o',capsize=3,label='IC Distribution')
plt.errorbar(t100[1:100],Dt8_J0_fixed[1:100],errDt8_J0_fixed[1:100],ls='',marker='^',capsize=3,label='IC Fixed')

plt.legend(prop={'size': 12})

plt.show()


#%%
                               ########### SAME FOR J=1 ##############           

#%%

FileRHF_J1=['/home/procyon/mfranco/KMC_simulations/Traces_100_Res_J=1.00_N=6.50.dat']
DataRHF_J1=[]

FileSHF_J1=['/home/procyon/mfranco/KMC_simulations/Traces_100_Stall_J=1.00_N=6.50.dat']
DataSHF_J1=[]

FileSHF_J1_fixed=['/home/procyon/mfranco/KMC_simulations/Traces_100_Stall_J=1.00_N=6.50_fixed.dat']
DataSHF_J1_fixed=[]


for data_file in FileRHF_J1:
    DataRHF_J1.append(np.loadtxt(data_file))
    
for data_file in FileSHF_J1:
    DataSHF_J1.append(np.loadtxt(data_file))

for data_file in FileSHF_J1_fixed:
    DataSHF_J1_fixed.append(np.loadtxt(data_file))
    
ResHF_J1=DataRHF_J1[0]

RelHF_J1=DataSHF_J1[0]

RelHF_J1_fixed=DataSHF_J1_fixed[0]


#%%

nmax=13

def N(x,tau,nss):
    return nss + (n0-nss)*np.exp(-x/tau)


nrep=100

tResHF_J1=np.zeros(nrep)
errtResHF_J1=np.zeros(nrep)
tRelHF_J1=np.zeros(nrep)
errtRelHF_J1=np.zeros(nrep)
tRelHF_J1_fixed=np.zeros(nrep)
errtRelHF_J1_fixed=np.zeros(nrep)
DtHF_J1=np.zeros(nrep)
errDtHF_J1=np.zeros(nrep)
DtHF_J1_fixed=np.zeros(nrep)
errDtHF_J1_fixed=np.zeros(nrep)

#%%

#ARRAYS FOR HALF-FILLING
        
j=0
len_ResHF_J1=0
while(ResHF_J1[j,0]==1.):
        
    len_ResHF_J1 = len_ResHF_J1 + 1
        
    j = j + 1 

tHFRes_J1 = np.zeros((nrep,len_ResHF_J1))
N_ResHF_J1 = np.zeros((nrep,len_ResHF_J1))



j=0
len_RelHF_J1=0
while(RelHF_J1[j,0]==1.):
        
    len_RelHF_J1 = len_RelHF_J1 + 1
        
    j = j + 1 

tHFRel_J1 = np.zeros((nrep,len_RelHF_J1))
N_RelHF_J1 = np.zeros((nrep,len_RelHF_J1))

j=0
len_RelHF_J1_fixed=0
while(RelHF_J1_fixed[j,0]==1.):
        
    len_RelHF_J1_fixed = len_RelHF_J1_fixed + 1
        
    j = j + 1 

tHFRel_J1_fixed = np.zeros((nrep,len_RelHF_J1_fixed))
N_RelHF_J1_fixed = np.zeros((nrep,len_RelHF_J1_fixed))



#%%

# HALF-FILLING 

for i in range(nrep):
    
    flag = float(i+1)
    
    kResHF_J1=0
    for k in range(len(ResHF_J1[:,0])):
        
        if ((ResHF_J1[k,0]==flag) and (kResHF_J1<len_ResHF_J1)):
            
            tHFRes_J1[i,kResHF_J1] = ResHF_J1[k,1]
            
            N_ResHF_J1[i,kResHF_J1] = ResHF_J1[k,2]
            
            kResHF_J1 = kResHF_J1 + 1
            
        if ResHF_J1[k,0]>flag:
            
            break
    
    kRelHF_J1=0
    for k in range(len(RelHF_J1[:,0])):
        
        if ((RelHF_J1[k,0]==flag) and (kRelHF_J1<len_RelHF_J1)):
            
            tHFRel_J1[i,kRelHF_J1] = RelHF_J1[k,1]
            
            N_RelHF_J1[i,kRelHF_J1] = RelHF_J1[k,2]
            
            kRelHF_J1 = kRelHF_J1 + 1
            
        if RelHF_J1[k,0]>flag:
            
            break
        
    kRelHF_J1_fixed=0
    for k in range(len(RelHF_J1_fixed[:,0])):
        
        if ((RelHF_J1_fixed[k,0]==flag) and (kRelHF_J1_fixed<len_RelHF_J1_fixed)):
            
            tHFRel_J1_fixed[i,kRelHF_J1_fixed] = RelHF_J1_fixed[k,1]
            
            N_RelHF_J1_fixed[i,kRelHF_J1_fixed] = RelHF_J1_fixed[k,2]
            
            kRelHF_J1_fixed = kRelHF_J1_fixed + 1
            
        if  RelHF_J1_fixed[k,0]>flag:
            
            break
        
#%% FIT RESURRECTION

n0=0.0

for i in range(nrep):
    poResHF_J1, pcResHF_J1 = curve_fit(N,tHFRes_J1[i,:],N_ResHF_J1[i,:])
    
    tResHF_J1[i] = poResHF_J1[0]
    errResHF_J1 = np.sqrt(np.diag(pcResHF_J1))
    errtResHF_J1[i] = errResHF_J1[0]


#%%

n0=9.0

for i in range(nrep):
    
    poRelHF_J1, pcRelHF_J1 = curve_fit(N,tHFRel_J1[i,:],N_RelHF_J1[i,:])
    poRelHF_J1_fixed, pcRelHF_J1_fixed = curve_fit(N,tHFRel_J1_fixed[i,:],N_RelHF_J1_fixed[i,:])
    
    tRelHF_J1[i] = poRelHF_J1[0]
    errRelHF_J1 = np.sqrt(np.diag(pcRelHF_J1))
    errtRelHF_J1[i] = errRelHF_J1[0]
    
    tRelHF_J1_fixed[i] = poRelHF_J1_fixed[0]
    errRelHF_J1_fixed = np.sqrt(np.diag(pcRelHF_J1_fixed))
    errtRelHF_J1_fixed[i] = errRelHF_J1_fixed[0]


#%%

for i in range(nrep):
    
    DtHF_J1[i]=(tRelHF_J1[i]-tResHF_J1[i])/(tRelHF_J1[i]+tResHF_J1[i])
    errDtHF_J1[i]=2./(tRelHF_J1[i]+tResHF_J1[i])**2*(tResHF_J1[i]*errtRelHF_J1[i] + tRelHF_J1[i]*errtResHF_J1[i])
    
    DtHF_J1_fixed[i]=(tRelHF_J1_fixed[i]-tResHF_J1[i])/(tRelHF_J1_fixed[i]+tResHF_J1[i])
    errDtHF_J1_fixed[i]=2./(tRelHF_J1_fixed[i]+tResHF_J1[i])**2*(tResHF_J1[i]*errtRelHF_J1_fixed[i] + tRelHF_J1_fixed[i]*errtResHF_J1[i])


#%%

plt.plot(tHFRes_J1[60,:],N_ResHF_J1[60,:])
plt.plot(tHFRel_J1[70,:],N_RelHF_J1[70,:])
plt.show()

#%%

J=np.linspace(0,5,6)

t100 = np.linspace(0,100,100)


plt.title(r'J=1, <$\phi$>=1/2, <$n_0$>=$n_0$=9',fontsize=16)
plt.xlabel('# KMC simulation',fontsize=14)
plt.ylabel(r'$\Delta$',fontsize=14)
plt.grid()
plt.errorbar(t100,DtHF_J1,errDtHF_J1,ls='',marker='o',capsize=3,label='IC Distribution')
plt.errorbar(t100,DtHF_J1_fixed,errDtHF_J1_fixed,ls='',marker='^',capsize=3, label='IC Fixed')

plt.legend(prop={'size': 12})

plt.show()



#%%
                                     ######### SAME FOR J=2 #########
 
#%%
                                     
FileRHF_J2=['/home/procyon/mfranco/KMC_simulations/Traces_100_Res_J=2.00_N=6.50.dat']
DataRHF_J2=[]

FileSHF_J2=['/home/procyon/mfranco/KMC_simulations/Traces_100_Stall_J=2.00_N=6.50.dat']
DataSHF_J2=[]

FileSHF_J2_fixed=['/home/procyon/mfranco/KMC_simulations/Traces_100_Stall_J=2.00_N=6.50_fixed.dat']
DataSHF_J2_fixed=[]


for data_file in FileRHF_J2:
    DataRHF_J2.append(np.loadtxt(data_file))
    
for data_file in FileSHF_J2:
    DataSHF_J2.append(np.loadtxt(data_file))

for data_file in FileSHF_J2_fixed:
    DataSHF_J2_fixed.append(np.loadtxt(data_file))
    
ResHF_J2=DataRHF_J2[0]

RelHF_J2=DataSHF_J2[0]

RelHF_J2_fixed=DataSHF_J2_fixed[0]


#%%

nmax=13

def N(x,tau,nss):
    return nss + (n0-nss)*np.exp(-x/tau)


nrep=100

tResHF_J2=np.zeros(nrep)
errtResHF_J2=np.zeros(nrep)
tRelHF_J2=np.zeros(nrep)
errtRelHF_J2=np.zeros(nrep)
tRelHF_J2_fixed=np.zeros(nrep)
errtRelHF_J2_fixed=np.zeros(nrep)
DtHF_J2=np.zeros(nrep)
errDtHF_J2=np.zeros(nrep)
DtHF_J2_fixed=np.zeros(nrep)
errDtHF_J2_fixed=np.zeros(nrep)

#%%

#ARRAYS FOR HALF-FILLING
        
j=0
len_ResHF_J2=0
while(ResHF_J2[j,0]==1.):
        
    len_ResHF_J2 = len_ResHF_J2 + 1
        
    j = j + 1 

tHF_J2 = np.zeros((nrep,len_ResHF_J2))
N_ResHF_J2 = np.zeros((nrep,len_ResHF_J2))



j=0
len_RelHF_J2=0
while(RelHF_J2[j,0]==1.):
        
    len_RelHF_J2 = len_RelHF_J2 + 1
        
    j = j + 1 

tHFRel_J2 = np.zeros((nrep,len_RelHF_J2))
N_RelHF_J2 = np.zeros((nrep,len_RelHF_J2))

j=0
len_RelHF_J2_fixed=0
while(RelHF_J2_fixed[j,0]==1.):
        
    len_RelHF_J2_fixed = len_RelHF_J2_fixed + 1
        
    j = j + 1 

tHFRel_J2_fixed = np.zeros((nrep,len_RelHF_J2_fixed))
N_RelHF_J2_fixed = np.zeros((nrep,len_RelHF_J2_fixed))



#%%

# HALF-FILLING 

for i in range(nrep):
    
    flag = float(i+1)
    
    kResHF_J2=0
    for k in range(len(ResHF_J2[:,0])):
        
        if ((ResHF_J2[k,0]==flag) and (kResHF_J2<len_ResHF_J2)):
            
            tHF_J2[i,kResHF_J2] = ResHF_J2[k,1]
            
            N_ResHF_J2[i,kResHF_J2] = ResHF_J2[k,2]
            
            kResHF_J2 = kResHF_J2 + 1
            
        if ResHF_J2[k,0]>flag:
            
            break
    
    kRelHF_J2=0
    for k in range(len(RelHF_J2[:,0])):
        
        if ((RelHF_J2[k,0]==flag) and (kRelHF_J2<len_RelHF_J2)):
            
            tHFRel_J2[i,kRelHF_J2] = RelHF_J2[k,1]
            
            N_RelHF_J2[i,kRelHF_J2] = RelHF_J2[k,2]
            
            kRelHF_J2 = kRelHF_J2 + 1
            
        if RelHF_J2[k,0]>flag:
            
            break
        
    kRelHF_J2_fixed=0
    for k in range(len(RelHF_J2_fixed[:,0])):
        
        if ((RelHF_J2_fixed[k,0]==flag) and (kRelHF_J2_fixed<len_RelHF_J2_fixed)):
            
            tHFRel_J2_fixed[i,kRelHF_J2_fixed] = RelHF_J2_fixed[k,1]
            
            N_RelHF_J2_fixed[i,kRelHF_J2_fixed] = RelHF_J2_fixed[k,2]
            
            kRelHF_J2_fixed = kRelHF_J2_fixed + 1
            
        if  RelHF_J2_fixed[k,0]>flag:
            
            break
        
#%% FIT RESURRECTION

n0=0.0

for i in range(nrep):
    poResHF_J2, pcResHF_J2 = curve_fit(N,tHF_J2[i,:],N_ResHF_J2[i,:])
    
    tResHF_J2[i] = poResHF_J2[0]
    errResHF_J2 = np.sqrt(np.diag(pcResHF_J2))
    errtResHF_J2[i] = errResHF_J2[0]


#%%

n0=9.0

for i in range(nrep):
    
    poRelHF_J2, pcRelHF_J2 = curve_fit(N,tHFRel_J2[i,:],N_RelHF_J2[i,:])
    poRelHF_J2_fixed, pcRelHF_J2_fixed = curve_fit(N,tHFRel_J2_fixed[i,:],N_RelHF_J2_fixed[i,:])
    
    tRelHF_J2[i] = poRelHF_J2[0]
    errRelHF_J2 = np.sqrt(np.diag(pcRelHF_J2))
    errtRelHF_J2[i] = errRelHF_J2[0]
    
    tRelHF_J2_fixed[i] = poRelHF_J2_fixed[0]
    errRelHF_J2_fixed = np.sqrt(np.diag(pcRelHF_J2_fixed))
    errtRelHF_J2_fixed[i] = errRelHF_J2_fixed[0]


#%%

for i in range(nrep):
    
    DtHF_J2[i]=(tRelHF_J2[i]-tResHF_J2[i])/(tRelHF_J2[i]+tResHF_J2[i])
    errDtHF_J2[i]=2./(tRelHF_J2[i]+tResHF_J2[i])**2*(tResHF_J2[i]*errtRelHF_J2[i] + tRelHF_J2[i]*errtResHF_J2[i])
    
    DtHF_J2_fixed[i]=(tRelHF_J2_fixed[i]-tResHF_J2[i])/(tRelHF_J2_fixed[i]+tResHF_J2[i])
    errDtHF_J2_fixed[i]=2./(tRelHF_J2_fixed[i]+tResHF_J2[i])**2*(tResHF_J2[i]*errtRelHF_J2_fixed[i] + tRelHF_J2_fixed[i]*errtResHF_J2[i])


#%%

J=np.linspace(0,5,6)

t100 = np.linspace(0,100,100)


plt.title(r'J=2, <$\phi$>=1/2, <$n_0$>=$n_0$=9',fontsize=16)
plt.xlabel('# KMC simulation',fontsize=14)
plt.ylabel(r'$\Delta$',fontsize=14)
plt.grid()
plt.errorbar(t100,DtHF_J2,errDtHF_J2,ls='',marker='o',capsize=3,label='IC Distribution')
plt.errorbar(t100,DtHF_J2_fixed,errDtHF_J2_fixed,ls='',marker='^',capsize=3, label='IC Fixed')

plt.legend(prop={'size': 12})

plt.show()


#%%
                                 ### SAME CODE FOR J=3 ####

#%% DATA READING

FileRHF_J3=['/home/procyon/mfranco/KMC_simulations/Traces_100_Res_J=3.00_N=6.50.dat']
DataRHF_J3=[]

FileSHF_J3=['/home/procyon/mfranco/KMC_simulations/Traces_100_Stall_J=3.00_N=6.50.dat']
DataSHF_J3=[]

FileSHF_J3_fixed=['/home/procyon/mfranco/KMC_simulations/Traces_100_Stall_J=3.00_N=6.50_fixed.dat']
DataSHF_J3_fixed=[]


for data_file in FileRHF_J3:
    DataRHF_J3.append(np.loadtxt(data_file))
    
for data_file in FileSHF_J3:
    DataSHF_J3.append(np.loadtxt(data_file))

for data_file in FileSHF_J3_fixed:
    DataSHF_J3_fixed.append(np.loadtxt(data_file))
    
ResHF_J3=DataRHF_J3[0]

RelHF_J3=DataSHF_J3[0]

RelHF_J3_fixed=DataSHF_J3_fixed[0]



FileR4_J3=['/home/procyon/mfranco/KMC_simulations/Traces_100_Res_J=3.00_N=4.00.dat']
DataR4_J3=[]

FileS4_J3=['/home/procyon/mfranco/KMC_simulations/Traces_100_Stall_J=3.00_N=4.00.dat']
DataS4_J3=[]

FileS4_J3_fixed=['/home/procyon/mfranco/KMC_simulations/Traces_100_Stall_J=3.00_N=4.00_fixed.dat']
DataS4_J3_fixed=[]


for data_file in FileR4_J3:
    DataR4_J3.append(np.loadtxt(data_file))
    
for data_file in FileS4_J3:
    DataS4_J3.append(np.loadtxt(data_file))

for data_file in FileS4_J3_fixed:
    DataS4_J3_fixed.append(np.loadtxt(data_file))
    
Res4_J3=DataR4_J3[0]

Rel4_J3=DataS4_J3[0]

Rel4_J3_fixed=DataS4_J3_fixed[0]



FileR8_J3=['/home/procyon/mfranco/KMC_simulations/Traces_100_Res_J=3.00_N=8.00.dat']
DataR8_J3=[]

FileS8_J3=['/home/procyon/mfranco/KMC_simulations/Traces_100_Stall_J=3.00_N=8.00.dat']
DataS8_J3=[]

FileS8_J3_fixed=['/home/procyon/mfranco/KMC_simulations/Traces_100_Stall_J=3.00_N=8.00_fixed.dat']
DataS8_J3_fixed=[]


for data_file in FileR8_J3:
    DataR8_J3.append(np.loadtxt(data_file))
    
for data_file in FileS8_J3:
    DataS8_J3.append(np.loadtxt(data_file))

for data_file in FileS8_J3_fixed:
    DataS8_J3_fixed.append(np.loadtxt(data_file))
    
Res8_J3=DataR8_J3[0]

Rel8_J3=DataS8_J3[0]

Rel8_J3_fixed=DataS8_J3_fixed[0]


#%%

nmax=13

def N(x,tau,nss):
    return nss + (n0-nss)*np.exp(-x/tau)


nrep=100

tResHF_J3=np.zeros(nrep)
errtResHF_J3=np.zeros(nrep)
tRelHF_J3=np.zeros(nrep)
errtRelHF_J3=np.zeros(nrep)
tRelHF_J3_fixed=np.zeros(nrep)
errtRelHF_J3_fixed=np.zeros(nrep)
DtHF_J3=np.zeros(nrep)
errDtHF_J3=np.zeros(nrep)
DtHF_J3_fixed=np.zeros(nrep)
errDtHF_J3_fixed=np.zeros(nrep)

tRes4_J3=np.zeros(nrep)
errtRes4_J3=np.zeros(nrep)
tRel4_J3=np.zeros(nrep)
errtRel4_J3=np.zeros(nrep)
tRel4_J3_fixed=np.zeros(nrep)
errtRel4_J3_fixed=np.zeros(nrep)
Dt4_J3=np.zeros(nrep)
errDt4_J3=np.zeros(nrep)
Dt4_J3_fixed=np.zeros(nrep)
errDt4_J3_fixed=np.zeros(nrep)

tRes8_J3=np.zeros(nrep)
errtRes8_J3=np.zeros(nrep)
tRel8_J3=np.zeros(nrep)
errtRel8_J3=np.zeros(nrep)
tRel8_J3_fixed=np.zeros(nrep)
errtRel8_J3_fixed=np.zeros(nrep)
Dt8_J3=np.zeros(nrep)
errDt8_J3=np.zeros(nrep)
Dt8_J3_fixed=np.zeros(nrep)
errDt8_J3_fixed=np.zeros(nrep)

#%%

#ARRAYS FOR HALF-FILLING
        
j=0
len_ResHF_J3=0
while(ResHF_J3[j,0]==1.):
        
    len_ResHF_J3 = len_ResHF_J3 + 1
        
    j = j + 1 

tHF_J3 = np.zeros((nrep,len_ResHF_J3))
N_ResHF_J3 = np.zeros((nrep,len_ResHF_J3))



j=0
len_RelHF_J3=0
while(RelHF_J3[j,0]==1.):
        
    len_RelHF_J3 = len_RelHF_J3 + 1
        
    j = j + 1 

tHFRel_J3 = np.zeros((nrep,len_RelHF_J3))
N_RelHF_J3 = np.zeros((nrep,len_RelHF_J3))

j=0
len_RelHF_J3_fixed=0
while(RelHF_J3_fixed[j,0]==1.):
        
    len_RelHF_J3_fixed = len_RelHF_J3_fixed + 1
        
    j = j + 1 

tHFRel_J3_fixed = np.zeros((nrep,len_RelHF_J3_fixed))
N_RelHF_J3_fixed = np.zeros((nrep,len_RelHF_J3_fixed))



# ARRAYS FOR N=4

j=0
len_Res4_J3=0
while(Res4_J3[j,0]==1.):
        
    len_Res4_J3 = len_Res4_J3 + 1
        
    j = j + 1 
    
t4Res_J3 = np.zeros((nrep,len_Res4_J3))
N_Res4_J3 = np.zeros((nrep,len_Res4_J3))

j=0
len_Rel4_J3=0
while(Rel4_J3[j,0]==1.):
        
    len_Rel4_J3 = len_Rel4_J3 + 1
        
    j = j + 1 

t4Rel_J3 = np.zeros((nrep,len_Rel4_J3))
N_Rel4_J3 = np.zeros((nrep,len_Rel4_J3))

j=0
len_Rel4_J3_fixed=0
while(Rel4_J3_fixed[j,0]==1.):
        
    len_Rel4_J3_fixed = len_Rel4_J3_fixed + 1
        
    j = j + 1 

t4Rel_J3_fixed = np.zeros((nrep,len_Rel4_J3_fixed))
N_Rel4_J3_fixed = np.zeros((nrep,len_Rel4_J3_fixed))


# ARRAYS FOR N=8

j=0
len_Res8_J3=0
while(Res8_J3[j,0]==1.):
        
    len_Res8_J3 = len_Res8_J3 + 1
        
    j = j + 1 
    
t8Res_J3 = np.zeros((nrep,len_Res8_J3))
N_Res8_J3 = np.zeros((nrep,len_Res8_J3))

j=0
len_Rel8_J3=0
while(Rel8_J3[j,0]==1.):
        
    len_Rel8_J3 = len_Rel8_J3 + 1
        
    j = j + 1 

len_Rel8_J3

t8Rel_J3 = np.zeros((nrep,len_Rel8_J3))
N_Rel8_J3 = np.zeros((nrep,len_Rel8_J3))


j=0
len_Rel8_J3_fixed=0
while(Rel8_J3_fixed[j,0]==1.):
        
    len_Rel8_J3_fixed = len_Rel8_J3_fixed + 1
        
    j = j + 1 

t8Rel_J3_fixed = np.zeros((nrep,len_Rel8_J3_fixed))
N_Rel8_J3_fixed = np.zeros((nrep,len_Rel8_J3_fixed))


#%%

# HALF-FILLING 

for i in range(nrep):
    
    flag = float(i+1)
    
    kResHF_J3=0
    for k in range(len(ResHF_J3[:,0])):
        
        if ((ResHF_J3[k,0]==flag) and (kResHF_J3<len_ResHF_J3)):
            
            tHF_J3[i,kResHF_J3] = ResHF_J3[k,1]
            
            N_ResHF_J3[i,kResHF_J3] = ResHF_J3[k,2]
            
            kResHF_J3 = kResHF_J3 + 1
            
        if ResHF_J3[k,0]>flag:
            
            break
    
    kRelHF_J3=0
    for k in range(len(RelHF_J3[:,0])):
        
        if ((RelHF_J3[k,0]==flag) and (kRelHF_J3<len_RelHF_J3)):
            
            tHFRel_J3[i,kRelHF_J3] = RelHF_J3[k,1]
            
            N_RelHF_J3[i,kRelHF_J3] = RelHF_J3[k,2]
            
            kRelHF_J3 = kRelHF_J3 + 1
            
        if RelHF_J3[k,0]>flag:
            
            break
        
    kRelHF_J3_fixed=0
    for k in range(len(RelHF_J3_fixed[:,0])):
        
        if ((RelHF_J3_fixed[k,0]==flag) and (kRelHF_J3_fixed<len_RelHF_J3_fixed)):
            
            tHFRel_J3_fixed[i,kRelHF_J3_fixed] = RelHF_J3_fixed[k,1]
            
            N_RelHF_J3_fixed[i,kRelHF_J3_fixed] = RelHF_J3_fixed[k,2]
            
            kRelHF_J3_fixed = kRelHF_J3_fixed + 1
            
        if  RelHF_J3_fixed[k,0]>flag:
            
            break
   
#%%

# N=4

for i in range(nrep):
    
    flag = float(i+1)
    
    kRes4_J3=0
    for k in range(len(Res4_J3[:,0])):
        
        if ((Res4_J3[k,0]==flag) and (kRes4_J3<len_Res4_J3)):
            
            t4Res_J3[i,kRes4_J3] = Res4_J3[k,1]
            
            N_Res4_J3[i,kRes4_J3] = Res4_J3[k,2]
            
            kRes4_J3 = kRes4_J3 + 1
            
        if Res4_J3[k,0]>flag:
            
            break
    
    kRel4_J3=0
    for k in range(len(Rel4_J3[:,0])):
        
        if ((Rel4_J3[k,0]==flag) and (kRel4_J3<len_Rel4_J3)):
            
            t4Rel_J3[i,kRel4_J3] = Rel4_J3[k,1]
            
            N_Rel4_J3[i,kRel4_J3] = Rel4_J3[k,2]
            
            kRel4_J3 = kRel4_J3 + 1
            
        if Rel4_J3[k,0]>flag:
            
            break
        
    kRel4_J3_fixed=0
    for k in range(len(Rel4_J3_fixed[:,0])):
        
        if ((Rel4_J3_fixed[k,0]==flag) and (kRel4_J3_fixed<len_Rel4_J3_fixed)):
            
            t4Rel_J3_fixed[i,kRel4_J3_fixed] = Rel4_J3_fixed[k,1]
            
            N_Rel4_J3_fixed[i,kRel4_J3_fixed] = Rel4_J3_fixed[k,2]
            
            kRel4_J3_fixed = kRel4_J3_fixed + 1
            
        if  Rel4_J3_fixed[k,0]>flag:
            
            break      

# N=8 
        
for i in range(nrep):
    
    flag = float(i+1)
    
    kRes8_J3=0
    for k in range(len(Res8_J3[:,0])):
        
        if ((Res8_J3[k,0]==flag) and (kRes8_J3<len_Res8_J3)):
            
            t8Res_J3[i,kRes8_J3] = Res8_J3[k,1]
            
            N_Res8_J3[i,kRes8_J3] = Res8_J3[k,2]
            
            kRes8_J3 = kRes8_J3 + 1
            
        if Res8_J3[k,0]>flag:
            
            break
    
    kRel8_J3=0
    for k in range(len(Rel8_J3[:,0])):
        
        if ((Rel8_J3[k,0]==flag) and (kRel8_J3<len_Rel8_J3)):
            
            t8Rel_J3[i,kRel8_J3] = Rel8_J3[k,1]
            
            N_Rel8_J3[i,kRel8_J3] = Rel8_J3[k,2]
            
            kRel8_J3 = kRel8_J3 + 1
            
        if Rel8_J3[k,0]>flag:
            
            break
        
    kRel8_J3_fixed=0
    for k in range(len(Rel8_J3_fixed[:,0])):
        
        if ((Rel8_J3_fixed[k,0]==flag) and (kRel8_J3_fixed<len_Rel8_J3_fixed)):
            
            t8Rel_J3_fixed[i,kRel8_J3_fixed] = Rel8_J3_fixed[k,1]
            
            N_Rel8_J3_fixed[i,kRel8_J3_fixed] = Rel8_J3_fixed[k,2]
            
            kRel8_J3_fixed = kRel8_J3_fixed + 1
            
        if  Rel8_J3_fixed[k,0]>flag:
            
            break
        
        
       
        
#%% FIT RESURRECTION

n0=0.0

for i in range(nrep):
    poResHF_J3, pcResHF_J3 = curve_fit(N,tHF_J3[i,:],N_ResHF_J3[i,:])
    
    tResHF_J3[i] = poResHF_J3[0]
    errResHF_J3 = np.sqrt(np.diag(pcResHF_J3))
    errtResHF_J3[i] = errResHF_J3[0]
    
    
    poRes4_J3, pcRes4_J3 = curve_fit(N,t4Res_J3[i,:],N_Res4_J3[i,:])
    
    tRes4_J3[i] = poRes4_J3[0]
    errRes4_J3 = np.sqrt(np.diag(pcRes4_J3))
    errtRes4_J3[i] = errRes4_J3[0]
    
    
    poRes8_J3, pcRes8_J3 = curve_fit(N,t8Res_J3[i,:],N_Res8_J3[i,:])
    
    tRes8_J3[i] = poRes8_J3[0]
    errRes8_J3 = np.sqrt(np.diag(pcRes8_J3))
    errtRes8_J3[i] = errRes8_J3[0]




#%%

n0=9.0

for i in range(nrep):
    
    poRelHF_J3, pcRelHF_J3 = curve_fit(N,tHFRel_J3[i,:],N_RelHF_J3[i,:])
    poRelHF_J3_fixed, pcRelHF_J3_fixed = curve_fit(N,tHFRel_J3_fixed[i,:],N_RelHF_J3_fixed[i,:])
    
    tRelHF_J3[i] = poRelHF_J3[0]
    errRelHF_J3 = np.sqrt(np.diag(pcRelHF_J3))
    errtRelHF_J3[i] = errRelHF_J3[0]
    
    tRelHF_J3_fixed[i] = poRelHF_J3_fixed[0]
    errRelHF_J3_fixed = np.sqrt(np.diag(pcRelHF_J3_fixed))
    errtRelHF_J3_fixed[i] = errRelHF_J3_fixed[0]
    
    
    poRel8_J3, pcRel8_J3 = curve_fit(N,t8Rel_J3[i,:],N_Rel8_J3[i,:])
    poRel8_J3_fixed, pcRel8_J3_fixed = curve_fit(N,t8Rel_J3_fixed[i,:],N_Rel8_J3_fixed[i,:])
    
    tRel8_J3[i] = poRel8_J3[0]
    errRel8_J3 = np.sqrt(np.diag(pcRel8_J3))
    errtRel8_J3[i] = errRel8_J3[0]
    
    tRel8_J3_fixed[i] = poRel8_J3_fixed[0]
    errRel8_J3_fixed = np.sqrt(np.diag(pcRel8_J3_fixed))
    errtRel8_J3_fixed[i] = errRel8_J3_fixed[0]

#%%

n0=6.0

for i in range(nrep):
    
    poRel4_J3, pcRel4_J3 = curve_fit(N,t4Rel_J3[i,:],N_Rel4_J3[i,:])
    poRel4_J3_fixed, pcRel4_J3_fixed = curve_fit(N,t4Rel_J3_fixed[i,:],N_Rel4_J3_fixed[i,:])
    
    tRel4_J3[i] = poRel4_J3[0]
    errRel4_J3 = np.sqrt(np.diag(pcRel4_J3))
    errtRel4_J3[i] = errRel4_J3[0]
    
    tRel4_J3_fixed[i] = poRel4_J3_fixed[0]
    errRel4_J3_fixed = np.sqrt(np.diag(pcRel4_J3_fixed))
    errtRel4_J3_fixed[i] = errRel4_J3_fixed[0]

#%%

for i in range(nrep):
    
    DtHF_J3[i]=(tRelHF_J3[i]-tResHF_J3[i])/(tRelHF_J3[i]+tResHF_J3[i])
    errDtHF_J3[i]=2./(tRelHF_J3[i]+tResHF_J3[i])**2*(tResHF_J3[i]*errtRelHF_J3[i] + tRelHF_J3[i]*errtResHF_J3[i])
    
    DtHF_J3_fixed[i]=(tRelHF_J3_fixed[i]-tResHF_J3[i])/(tRelHF_J3_fixed[i]+tResHF_J3[i])
    errDtHF_J3_fixed[i]=2./(tRelHF_J3_fixed[i]+tResHF_J3[i])**2*(tResHF_J3[i]*errtRelHF_J3_fixed[i] + tRelHF_J3_fixed[i]*errtResHF_J3[i])
    
    
    Dt4_J3[i]=(tRel4_J3[i]-tRes4_J3[i])/(tRel4_J3[i]+tRes4_J3[i])
    errDt4_J3[i]=2./(tRel4_J3[i]+tRes4_J3[i])**2*(tRes4_J3[i]*errtRel4_J3[i] + tRel4_J3[i]*errtRes4_J3[i])
    
    Dt4_J3_fixed[i]=(tRel4_J3_fixed[i]-tRes4_J3[i])/(tRel4_J3_fixed[i]+tRes4_J3[i])
    errDt4_J3_fixed[i]=2./(tRel4_J3_fixed[i]+tRes4_J3[i])**2*(tRes4_J3[i]*errtRel4_J3_fixed[i] + tRel4_J3_fixed[i]*errtRes4_J3[i])


    Dt8_J3[i]=(tRel8_J3[i]-tRes8_J3[i])/(tRel8_J3[i]+tRes8_J3[i])
    errDt8_J3[i]=2./(tRel8_J3[i]+tRes8_J3[i])**2*(tRes8_J3[i]*errtRel8_J3[i] + tRel8_J3[i]*errtRes8_J3[i])
    
    Dt8_J3_fixed[i]=(tRel8_J3_fixed[i]-tRes8_J3[i])/(tRel8_J3_fixed[i]+tRes8_J3[i])
    errDt8_J3_fixed[i]=2./(tRel8_J3_fixed[i]+tRes8_J3[i])**2*(tRes8_J3[i]*errtRel8_J3_fixed[i] + tRel8_J3_fixed[i]*errtRes8_J3[i])



#%%

J=np.linspace(0,5,6)

t100 = np.linspace(0,100,100)


plt.title(r'J=3, <$\phi$>=1/2, <$n_0$>=$n_0$=9',fontsize=16)
plt.xlabel('# KMC simulation',fontsize=14)
plt.ylabel(r'$\Delta$',fontsize=14)
plt.grid()
plt.errorbar(t100,DtHF_J3,errDtHF_J3,ls='',marker='o',capsize=3,label='IC Distribution')
plt.errorbar(t100,DtHF_J3_fixed,errDtHF_J3_fixed,ls='',marker='^',capsize=3, label='IC Fixed')

plt.legend(prop={'size': 12})

plt.show()


plt.title(r'J=3, <$\phi$>=4/13, <$n_0$>=$n_0$=6',fontsize=16)
plt.xlabel('# KMC simulation',fontsize=14)
plt.ylabel(r'$\Delta$',fontsize=14)
plt.grid()
plt.errorbar(t100,Dt4_J3,errDt4_J3,ls='',marker='o',capsize=3,label='IC Distribution')
plt.errorbar(t100,Dt4_J3_fixed,errDt4_J3_fixed,ls='',marker='^',capsize=3,label='IC Fixed')

plt.legend(prop={'size': 12})

plt.show()


plt.title(r'J=3, <$\phi$>=8/13, <$n_0$>=$n_0$=9',fontsize=16)
plt.xlabel('# KMC Simulation',fontsize=14)
plt.ylabel(r'$\Delta$',fontsize=14)
plt.grid()
plt.errorbar(t100[1:100],Dt8_J3[1:100],errDt8_J3[1:100],ls='',marker='o',capsize=3,label='IC Distribution')
plt.errorbar(t100[1:100],Dt8_J3_fixed[1:100],errDt8_J3_fixed[1:100],ls='',marker='^',capsize=3,label='IC Fixed')

plt.legend(prop={'size': 12})

plt.show()

#%%
                            ########### SAME FOR J=4 ################

#%%

FileRHF_J4=['/home/procyon/mfranco/KMC_simulations/Traces_100_Res_J=4.00_N=6.50.dat']
DataRHF_J4=[]

FileSHF_J4=['/home/procyon/mfranco/KMC_simulations/Traces_100_Stall_J=4.00_N=6.50.dat']
DataSHF_J4=[]

FileSHF_J4_fixed=['/home/procyon/mfranco/KMC_simulations/Traces_100_Stall_J=4.00_N=6.50_fixed.dat']
DataSHF_J4_fixed=[]


for data_file in FileRHF_J4:
    DataRHF_J4.append(np.loadtxt(data_file))
    
for data_file in FileSHF_J4:
    DataSHF_J4.append(np.loadtxt(data_file))

for data_file in FileSHF_J4_fixed:
    DataSHF_J4_fixed.append(np.loadtxt(data_file))
    
ResHF_J4=DataRHF_J4[0]

RelHF_J4=DataSHF_J4[0]

RelHF_J4_fixed=DataSHF_J4_fixed[0]


#%%

nmax=13

def N(x,tau,nss):
    return nss + (n0-nss)*np.exp(-x/tau)


nrep=100

tResHF_J4=np.zeros(nrep)
errtResHF_J4=np.zeros(nrep)
tRelHF_J4=np.zeros(nrep)
errtRelHF_J4=np.zeros(nrep)
tRelHF_J4_fixed=np.zeros(nrep)
errtRelHF_J4_fixed=np.zeros(nrep)
DtHF_J4=np.zeros(nrep)
errDtHF_J4=np.zeros(nrep)
DtHF_J4_fixed=np.zeros(nrep)
errDtHF_J4_fixed=np.zeros(nrep)

#%%

#ARRAYS FOR HALF-FILLING
        
j=0
len_ResHF_J4=0
while(ResHF_J4[j,0]==1.):
        
    len_ResHF_J4 = len_ResHF_J4 + 1
        
    j = j + 1 

tHF_J4 = np.zeros((nrep,len_ResHF_J4))
N_ResHF_J4 = np.zeros((nrep,len_ResHF_J4))



j=0
len_RelHF_J4=0
while(RelHF_J4[j,0]==1.):
        
    len_RelHF_J4 = len_RelHF_J4 + 1
        
    j = j + 1 

tHFRel_J4 = np.zeros((nrep,len_RelHF_J4))
N_RelHF_J4 = np.zeros((nrep,len_RelHF_J4))

j=0
len_RelHF_J4_fixed=0
while(RelHF_J4_fixed[j,0]==1.):
        
    len_RelHF_J4_fixed = len_RelHF_J4_fixed + 1
        
    j = j + 1 

tHFRel_J4_fixed = np.zeros((nrep,len_RelHF_J4_fixed))
N_RelHF_J4_fixed = np.zeros((nrep,len_RelHF_J4_fixed))



#%%

# HALF-FILLING 

for i in range(nrep):
    
    flag = float(i+1)
    
    kResHF_J4=0
    for k in range(len(ResHF_J4[:,0])):
        
        if ((ResHF_J4[k,0]==flag) and (kResHF_J4<len_ResHF_J4)):
            
            tHF_J4[i,kResHF_J4] = ResHF_J4[k,1]
            
            N_ResHF_J4[i,kResHF_J4] = ResHF_J4[k,2]
            
            kResHF_J4 = kResHF_J4 + 1
            
        if ResHF_J4[k,0]>flag:
            
            break
    
    kRelHF_J4=0
    for k in range(len(RelHF_J4[:,0])):
        
        if ((RelHF_J4[k,0]==flag) and (kRelHF_J4<len_RelHF_J4)):
            
            tHFRel_J4[i,kRelHF_J4] = RelHF_J4[k,1]
            
            N_RelHF_J4[i,kRelHF_J4] = RelHF_J4[k,2]
            
            kRelHF_J4 = kRelHF_J4 + 1
            
        if RelHF_J4[k,0]>flag:
            
            break
        
    kRelHF_J4_fixed=0
    for k in range(len(RelHF_J4_fixed[:,0])):
        
        if ((RelHF_J4_fixed[k,0]==flag) and (kRelHF_J4_fixed<len_RelHF_J4_fixed)):
            
            tHFRel_J4_fixed[i,kRelHF_J4_fixed] = RelHF_J4_fixed[k,1]
            
            N_RelHF_J4_fixed[i,kRelHF_J4_fixed] = RelHF_J4_fixed[k,2]
            
            kRelHF_J4_fixed = kRelHF_J4_fixed + 1
            
        if  RelHF_J4_fixed[k,0]>flag:
            
            break
        
#%% FIT RESURRECTION

n0=0.0

for i in range(nrep):
    poResHF_J4, pcResHF_J4 = curve_fit(N,tHF_J4[i,:],N_ResHF_J4[i,:])
    
    tResHF_J4[i] = poResHF_J4[0]
    errResHF_J4 = np.sqrt(np.diag(pcResHF_J4))
    errtResHF_J4[i] = errResHF_J4[0]


#%%

n0=9.0

for i in range(nrep):
    
    poRelHF_J4, pcRelHF_J4 = curve_fit(N,tHFRel_J4[i,:],N_RelHF_J4[i,:])
    poRelHF_J4_fixed, pcRelHF_J4_fixed = curve_fit(N,tHFRel_J4_fixed[i,:],N_RelHF_J4_fixed[i,:])
    
    tRelHF_J4[i] = poRelHF_J4[0]
    errRelHF_J4 = np.sqrt(np.diag(pcRelHF_J4))
    errtRelHF_J4[i] = errRelHF_J4[0]
    
    tRelHF_J4_fixed[i] = poRelHF_J4_fixed[0]
    errRelHF_J4_fixed = np.sqrt(np.diag(pcRelHF_J4_fixed))
    errtRelHF_J4_fixed[i] = errRelHF_J4_fixed[0]


#%%

for i in range(nrep):
    
    DtHF_J4[i]=(tRelHF_J4[i]-tResHF_J4[i])/(tRelHF_J4[i]+tResHF_J4[i])
    errDtHF_J4[i]=2./(tRelHF_J4[i]+tResHF_J4[i])**2*(tResHF_J4[i]*errtRelHF_J4[i] + tRelHF_J4[i]*errtResHF_J4[i])
    
    DtHF_J4_fixed[i]=(tRelHF_J4_fixed[i]-tResHF_J4[i])/(tRelHF_J4_fixed[i]+tResHF_J4[i])
    errDtHF_J4_fixed[i]=2./(tRelHF_J4_fixed[i]+tResHF_J4[i])**2*(tResHF_J4[i]*errtRelHF_J4_fixed[i] + tRelHF_J4_fixed[i]*errtResHF_J4[i])


#%%

J=np.linspace(0,5,6)

t100 = np.linspace(0,100,100)


plt.title(r'J=4, <$\phi$>=1/2, <$n_0$>=$n_0$=9',fontsize=16)
plt.xlabel('# KMC simulation',fontsize=14)
plt.ylabel(r'$\Delta$',fontsize=14)
plt.grid()
plt.errorbar(t100,DtHF_J4,errDtHF_J4,ls='',marker='o',capsize=3,label='IC Distribution')
plt.errorbar(t100,DtHF_J4_fixed,errDtHF_J4_fixed,ls='',marker='^',capsize=3, label='IC Fixed')

plt.legend(prop={'size': 12})

plt.show()



#%%
                     ############### SAME FOR J=5 ###################

#%%

FileRHF_J5=['/home/procyon/mfranco/KMC_simulations/Traces_100_Res_J=5.00_N=6.50.dat']
DataRHF_J5=[]

FileSHF_J5=['/home/procyon/mfranco/KMC_simulations/Traces_100_Stall_J=5.00_N=6.50.dat']
DataSHF_J5=[]

FileSHF_J5_fixed=['/home/procyon/mfranco/KMC_simulations/Traces_100_Stall_J=5.00_N=6.50_fixed.dat']
DataSHF_J5_fixed=[]


for data_file in FileRHF_J5:
    DataRHF_J5.append(np.loadtxt(data_file))
    
for data_file in FileSHF_J5:
    DataSHF_J5.append(np.loadtxt(data_file))

for data_file in FileSHF_J5_fixed:
    DataSHF_J5_fixed.append(np.loadtxt(data_file))
    
ResHF_J5=DataRHF_J5[0]

RelHF_J5=DataSHF_J5[0]

RelHF_J5_fixed=DataSHF_J5_fixed[0]


#%%

nmax=13

def N(x,tau,nss):
    return nss + (n0-nss)*np.exp(-x/tau)


nrep=100

tResHF_J5=np.zeros(nrep)
errtResHF_J5=np.zeros(nrep)
tRelHF_J5=np.zeros(nrep)
errtRelHF_J5=np.zeros(nrep)
tRelHF_J5_fixed=np.zeros(nrep)
errtRelHF_J5_fixed=np.zeros(nrep)
DtHF_J5=np.zeros(nrep)
errDtHF_J5=np.zeros(nrep)
DtHF_J5_fixed=np.zeros(nrep)
errDtHF_J5_fixed=np.zeros(nrep)

#%%

#ARRAYS FOR HALF-FILLING
        
j=0
len_ResHF_J5=0
while(ResHF_J5[j,0]==1.):
        
    len_ResHF_J5 = len_ResHF_J5 + 1
        
    j = j + 1 

tHF_J5 = np.zeros((nrep,len_ResHF_J5))
N_ResHF_J5 = np.zeros((nrep,len_ResHF_J5))



j=0
len_RelHF_J5=0
while(RelHF_J5[j,0]==1.):
        
    len_RelHF_J5 = len_RelHF_J5 + 1
        
    j = j + 1 

tHFRel_J5 = np.zeros((nrep,len_RelHF_J5))
N_RelHF_J5 = np.zeros((nrep,len_RelHF_J5))

j=0
len_RelHF_J5_fixed=0
while(RelHF_J5_fixed[j,0]==1.):
        
    len_RelHF_J5_fixed = len_RelHF_J5_fixed + 1
        
    j = j + 1 

tHFRel_J5_fixed = np.zeros((nrep,len_RelHF_J5_fixed))
N_RelHF_J5_fixed = np.zeros((nrep,len_RelHF_J5_fixed))



#%%

# HALF-FILLING 

for i in range(nrep):
    
    flag = float(i+1)
    
    kResHF_J5=0
    for k in range(len(ResHF_J5[:,0])):
        
        if ((ResHF_J5[k,0]==flag) and (kResHF_J5<len_ResHF_J5)):
            
            tHF_J5[i,kResHF_J5] = ResHF_J5[k,1]
            
            N_ResHF_J5[i,kResHF_J5] = ResHF_J5[k,2]
            
            kResHF_J5 = kResHF_J5 + 1
            
        if ResHF_J5[k,0]>flag:
            
            break
    
    kRelHF_J5=0
    for k in range(len(RelHF_J5[:,0])):
        
        if ((RelHF_J5[k,0]==flag) and (kRelHF_J5<len_RelHF_J5)):
            
            tHFRel_J5[i,kRelHF_J5] = RelHF_J5[k,1]
            
            N_RelHF_J5[i,kRelHF_J5] = RelHF_J5[k,2]
            
            kRelHF_J5 = kRelHF_J5 + 1
            
        if RelHF_J5[k,0]>flag:
            
            break
        
    kRelHF_J5_fixed=0
    for k in range(len(RelHF_J5_fixed[:,0])):
        
        if ((RelHF_J5_fixed[k,0]==flag) and (kRelHF_J5_fixed<len_RelHF_J5_fixed)):
            
            tHFRel_J5_fixed[i,kRelHF_J5_fixed] = RelHF_J5_fixed[k,1]
            
            N_RelHF_J5_fixed[i,kRelHF_J5_fixed] = RelHF_J5_fixed[k,2]
            
            kRelHF_J5_fixed = kRelHF_J5_fixed + 1
            
        if  RelHF_J5_fixed[k,0]>flag:
            
            break
        
#%% FIT RESURRECTION

n0=0.0

for i in range(nrep):
    poResHF_J5, pcResHF_J5 = curve_fit(N,tHF_J5[i,:],N_ResHF_J5[i,:])
    
    tResHF_J5[i] = poResHF_J5[0]
    errResHF_J5 = np.sqrt(np.diag(pcResHF_J5))
    errtResHF_J5[i] = errResHF_J5[0]


#%%

n0=9.0

for i in range(nrep):
    
    poRelHF_J5, pcRelHF_J5 = curve_fit(N,tHFRel_J5[i,:],N_RelHF_J5[i,:])
    poRelHF_J5_fixed, pcRelHF_J5_fixed = curve_fit(N,tHFRel_J5_fixed[i,:],N_RelHF_J5_fixed[i,:])
    
    tRelHF_J5[i] = poRelHF_J5[0]
    errRelHF_J5 = np.sqrt(np.diag(pcRelHF_J5))
    errtRelHF_J5[i] = errRelHF_J5[0]
    
    tRelHF_J5_fixed[i] = poRelHF_J5_fixed[0]
    errRelHF_J5_fixed = np.sqrt(np.diag(pcRelHF_J5_fixed))
    errtRelHF_J5_fixed[i] = errRelHF_J5_fixed[0]


#%%

for i in range(nrep):
    
    DtHF_J5[i]=(tRelHF_J5[i]-tResHF_J5[i])/(tRelHF_J5[i]+tResHF_J5[i])
    errDtHF_J5[i]=2./(tRelHF_J5[i]+tResHF_J5[i])**2*(tResHF_J5[i]*errtRelHF_J5[i] + tRelHF_J5[i]*errtResHF_J5[i])
    
    DtHF_J5_fixed[i]=(tRelHF_J5_fixed[i]-tResHF_J5[i])/(tRelHF_J5_fixed[i]+tResHF_J5[i])
    errDtHF_J5_fixed[i]=2./(tRelHF_J5_fixed[i]+tResHF_J5[i])**2*(tResHF_J5[i]*errtRelHF_J5_fixed[i] + tRelHF_J5_fixed[i]*errtResHF_J5[i])


#%%

J=np.linspace(0,5,6)

t100 = np.linspace(0,100,100)


plt.title(r'J=5, <$\phi$>=1/2, <$n_0$>=$n_0$=9',fontsize=16)
plt.xlabel('# KMC simulation',fontsize=14)
plt.ylabel(r'$\Delta$',fontsize=14)
plt.grid()
plt.errorbar(t100,DtHF_J5,errDtHF_J5,ls='',marker='o',capsize=3,label='IC Distribution')
plt.errorbar(t100,DtHF_J5_fixed,errDtHF_J5_fixed,ls='',marker='^',capsize=3, label='IC Fixed')

plt.legend(prop={'size': 12})

plt.show()







#%% AVERAGE OF RELAXATION TIMES

avtResHF_J0 = sum(tResHF_J0)/nrep
erravtResHF_J0=np.sqrt(sum((tResHF_J0-avtResHF_J0)**2)/nrep)

avtRelHF_J0 = sum(tRelHF_J0)/nrep
erravtRelHF_J0=np.sqrt(sum((tRelHF_J0-avtRelHF_J0)**2)/nrep)

avtRelHF_J0_fixed = sum(tRelHF_J0_fixed)/nrep
erravtRelHF_J0_fixed=np.sqrt(sum((tRelHF_J0_fixed-avtRelHF_J0_fixed)**2)/nrep)


avtResHF_J1 = sum(tResHF_J1)/nrep
erravtResHF_J1=np.sqrt(sum((tResHF_J1-avtResHF_J1)**2)/nrep)

avtRelHF_J1 = sum(tRelHF_J1)/nrep
erravtRelHF_J1=np.sqrt(sum((tRelHF_J1-avtRelHF_J1)**2)/nrep)

avtRelHF_J1_fixed = sum(tRelHF_J1_fixed)/nrep
erravtRelHF_J1_fixed=np.sqrt(sum((tRelHF_J1_fixed-avtRelHF_J1_fixed)**2)/nrep)


avtResHF_J2 = sum(tResHF_J2)/nrep
erravtResHF_J2=np.sqrt(sum((tResHF_J2-avtResHF_J2)**2)/nrep)

avtRelHF_J2 = sum(tRelHF_J2)/nrep
erravtRelHF_J2=np.sqrt(sum((tRelHF_J2-avtRelHF_J2)**2)/nrep)

avtRelHF_J2_fixed = sum(tRelHF_J2_fixed)/nrep
erravtRelHF_J2_fixed=np.sqrt(sum((tRelHF_J2_fixed-avtRelHF_J2_fixed)**2)/nrep)


avtResHF_J3 = sum(tResHF_J3)/nrep
erravtResHF_J3=np.sqrt(sum((tResHF_J3-avtResHF_J3)**2)/nrep)

avtRelHF_J3 = sum(tRelHF_J3)/nrep
erravtRelHF_J3=np.sqrt(sum((tRelHF_J3-avtRelHF_J3)**2)/nrep)

avtRelHF_J3_fixed = sum(tRelHF_J3_fixed)/nrep
erravtRelHF_J3_fixed=np.sqrt(sum((tRelHF_J3_fixed-avtRelHF_J3_fixed)**2)/nrep)


avtResHF_J4 = sum(tResHF_J4)/nrep
erravtResHF_J4=np.sqrt(sum((tResHF_J4-avtResHF_J4)**2)/nrep)

avtRelHF_J4 = sum(tRelHF_J4)/nrep
erravtRelHF_J4=np.sqrt(sum((tRelHF_J4-avtRelHF_J4)**2)/nrep)

avtRelHF_J4_fixed = sum(tRelHF_J4_fixed)/nrep
erravtRelHF_J4_fixed=np.sqrt(sum((tRelHF_J4_fixed-avtRelHF_J4_fixed)**2)/nrep)


avtResHF_J5 = sum(tResHF_J5)/nrep
erravtResHF_J5=np.sqrt(sum((tResHF_J5-avtResHF_J5)**2)/nrep)

avtRelHF_J5 = sum(tRelHF_J5)/nrep
erravtRelHF_J5=np.sqrt(sum((tRelHF_J5-avtRelHF_J5)**2)/nrep)

avtRelHF_J5_fixed = sum(tRelHF_J5_fixed)/nrep
erravtRelHF_J5_fixed=np.sqrt(sum((tRelHF_J5_fixed-avtRelHF_J5_fixed)**2)/nrep)

#%%

avtResHF = np.zeros(6)
avtRelHF = np.zeros(6)
avtRelHF_fixed = np.zeros(6)

erravtResHF = np.zeros(6)
erravtRelHF = np.zeros(6)
erravtRelHF_fixed = np.zeros(6)

#%%

avtResHF[0] = avtResHF_J0
avtResHF[1] = avtResHF_J1
avtResHF[2] = avtResHF_J2
avtResHF[3] = avtResHF_J3
avtResHF[4] = avtResHF_J4
avtResHF[5] = avtResHF_J5

avtRelHF[0] = avtResHF_J0
avtRelHF[1] = avtResHF_J1
avtRelHF[2] = avtResHF_J2
avtRelHF[3] = avtResHF_J3
avtRelHF[4] = avtResHF_J4
avtRelHF[5] = avtResHF_J5

avtRelHF_fixed[0] = avtResHF_J0
avtRelHF_fixed[1] = avtResHF_J1
avtRelHF_fixed[2] = avtResHF_J2
avtRelHF_fixed[3] = avtResHF_J3
avtRelHF_fixed[4] = avtResHF_J4
avtRelHF_fixed[5] = avtResHF_J5


erravtResHF[0] = avtResHF_J0
erravtResHF[1] = avtResHF_J1
erravtResHF[2] = avtResHF_J2
erravtResHF[3] = avtResHF_J3
erravtResHF[4] = avtResHF_J4
erravtResHF[5] = avtResHF_J5

erravtRelHF[0] = avtResHF_J0
erravtRelHF[1] = avtResHF_J1
erravtRelHF[2] = avtResHF_J2
erravtRelHF[3] = avtResHF_J3
erravtRelHF[4] = avtResHF_J4
erravtRelHF[5] = avtResHF_J5

erravtRelHF_fixed[0] = avtResHF_J0
erravtRelHF_fixed[1] = avtResHF_J1
erravtRelHF_fixed[2] = avtResHF_J2
erravtRelHF_fixed[3] = avtResHF_J3
erravtRelHF_fixed[4] = avtResHF_J4
erravtRelHF_fixed[5] = avtResHF_J5

#%%

def t_th(x):
    return np.exp(x)/2

J = np.arange(0,6,1)
cJ = np.arange(0,5,0.001)

plt.title(r'Relaxation times for $<\phi>=1/2$ (100 simulations)',fontsize=16)
plt.xlabel('J',fontsize=14)
plt.ylabel(r'$log(\tau)$',fontsize=14)
plt.grid()


plt.errorbar(J,np.log(avtResHF),erravtResHF/avtResHF,ls='',marker='s',label='Resurrection',capsize=3)
plt.errorbar(J,np.log(avtRelHF),erravtRelHF/avtRelHF,ls='',marker='o',label='Release ($<n_0>=9$)',capsize=3)
plt.errorbar(J,np.log(avtRelHF_fixed),erravtRelHF_fixed/avtRelHF_fixed,ls='',marker='^',label='Release ($n_0=9$)',capsize=3)
plt.plot(cJ,np.log(t_th(cJ)),label=r'$log(\tau) \propto J$ ')

plt.legend(prop={'size': 12})

#%% AVERAGE OF RELATIVE DIFFERENCE OF RELAXATION TIMES

avDtHF_J0 = sum(DtHF_J0)/nrep
avDtHF_J0_fixed = sum(DtHF_J0_fixed)/nrep

avDt4_J0 = sum(Dt4_J0)/nrep
avDt4_J0_fixed = sum(Dt4_J0_fixed)/nrep

avDt8_J0 = sum(Dt8_J0)/nrep
avDt8_J0_fixed = sum(Dt8_J0_fixed)/nrep

erravDtHF_J0=np.sqrt(sum((DtHF_J0-avDtHF_J0)**2)/nrep)
erravDtHF_J0_fixed=np.sqrt(sum((DtHF_J0_fixed-avDtHF_J0_fixed)**2)/nrep)

erravDt4_J0=np.sqrt(sum((Dt4_J0-avDt4_J0)**2)/nrep)
erravDt4_J0_fixed=np.sqrt(sum((Dt4_J0_fixed-avDt4_J0_fixed)**2)/nrep)

erravDt8_J0=np.sqrt(sum((Dt8_J0-avDt8_J0)**2)/nrep)
erravDt8_J0_fixed=np.sqrt(sum((Dt8_J0_fixed-avDt8_J0_fixed)**2)/nrep)

avDtHF_J3 = sum(DtHF_J3)/nrep
avDtHF_J3_fixed = sum(DtHF_J3_fixed)/nrep

avDt4_J3 = sum(Dt4_J3)/nrep
avDt4_J3_fixed = sum(Dt4_J3_fixed)/nrep

avDt8_J3 = sum(Dt8_J3)/nrep
avDt8_J3_fixed = sum(Dt8_J3_fixed)/nrep

erravDtHF_J3=np.sqrt(sum((DtHF_J3-avDtHF_J3)**2)/nrep)
erravDtHF_J3_fixed=np.sqrt(sum((DtHF_J3_fixed-avDtHF_J3_fixed)**2)/nrep)

erravDt4_J3=np.sqrt(sum((Dt4_J3-avDt4_J3)**2)/nrep)
erravDt4_J3_fixed=np.sqrt(sum((Dt4_J3_fixed-avDt4_J3_fixed)**2)/nrep)

erravDt8_J3=np.sqrt(sum((Dt8_J3-avDt8_J3)**2)/nrep)
erravDt8_J3_fixed=np.sqrt(sum((Dt8_J3_fixed-avDt8_J3_fixed)**2)/nrep)

#%%
avDtHF = np.zeros(2)
erravDtHF = np.zeros(2)
avDtHF_fixed = np.zeros(2)
erravDtHF_fixed = np.zeros(2)

avDt4=np.zeros(2)
erravDt4=np.zeros(2)

avDt8 = np.zeros(2)
erravDt8 = np.zeros(2)
#%%

avDtHF[0] = avDtHF_J0
avDtHF[1] = avDtHF_J3


#%%

J0=np.linspace(0,3,2)

plt.errorbar(J0,avDtHF_J0,erravDtHF_J0,ls='',marker='o',capsize=3,color='royalblue')
plt.errorbar(J0,avDtHF_J0_fixed,erravDtHF_J0_fixed,ls='',marker='^',capsize=3,color='royalblue')

plt.show()

plt.errorbar(J0,avDt4_J0,erravDt4_J0,ls='',marker='o',capsize=3,color='darkorange')
plt.errorbar(J0,avDt4_J0_fixed,erravDt4_J0_fixed,ls='',marker='^',capsize=3,color='darkorange')
plt.errorbar(J0,avDt8_J0,erravDt8_J0,ls='',marker='o',capsize=3,color='purple')
plt.errorbar(J0,avDt8_J0_fixed,erravDt8_J0_fixed,ls='',marker='^',capsize=3,color='purple')