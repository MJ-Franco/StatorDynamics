#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 13:06:14 2022

@author: mfranco
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
plt.rcParams["figure.figsize"] = [8.0,6.0]

nmax=13

h=10000
l=20000

def N_G(t,J,alpha,n0):
    gamma=np.tanh(J/2)
    
    return 1/2 + (n0 - 1/2)*np.exp(-alpha*(1-gamma)*t) 

t0=np.arange(0,15,1e-2)
t3=np.arange(0,100,1e-3)
t5=np.arange(0,1000,5e-3)

N_J0_d9 = np.zeros((h,len(t0)))
mN_J0_d9 = np.zeros(len(t0))

N_J3_d9 = np.zeros((h,len(t3)))
mN_J3_d9 = np.zeros(len(t3))

N_J5_d9 = np.zeros((h,len(t5)))
mN_J5_d9 = np.zeros(len(t5))


tRes=np.zeros(3)
tRel_f9=np.zeros(3)
tRel_d9=np.zeros(3)

tRel_f13=np.zeros(3)


#%%

files9_J0=['/media/mfranco/Elements/Simulations/Without _depletion/Glauber/Initial_conditions_J=0.00_N=9.dat']
data9_J0=[]
files9_J3=['/media/mfranco/Elements/Simulations/Without _depletion/Glauber/Initial_conditions_J=3.00_N=9.dat']
data9_J3=[]
files9_J5=['/media/mfranco/Elements/Simulations/Without _depletion/Glauber/Initial_conditions_J=5.00_N=9.dat']
data9_J5=[]

for data_file in files9_J0:
    data9_J0.append(np.loadtxt(data_file))
for data_file in files9_J3:
    data9_J3.append(np.loadtxt(data_file))
for data_file in files9_J5:
    data9_J5.append(np.loadtxt(data_file))
    
IC9_J0=data9_J0[0]
IC9_J3=data9_J3[0]
IC9_J5=data9_J5[0]



#%%

N_J0_0 = N_G(t0,0,1,0)
N_J0_f9 = N_G(t0,0,1,9/nmax)

N_J3_0 = N_G(t3,3,1,0)
N_J3_f9 = N_G(t3,3,1,9/nmax)

N_J5_0 = N_G(t5,5,1,0)
N_J5_f9 = N_G(t5,5,1,9/nmax)

#%%

for i in range(h):
    
    N_J0_d9[i] = N_G(t0,0,1,IC9_J0[i])
    

for j in range(len(t0)):
    sumN_J0=0
    for i in range(h):
        
        sumN_J0 = sumN_J0 + N_J0_d9[i,j] 
        
    mN_J0_d9[j] = sumN_J0/h
    
#%%
 
for i in range(h):
    
    N_J3_d9[i] = N_G(t3,3,1,IC9_J3[i])
    

for j in range(len(t3)):
    sumN_J3=0
    for i in range(h):
        
        sumN_J3 = sumN_J3 + N_J3_d9[i,j] 
            
    mN_J3_d9[j] = sumN_J3/h
    
#%%
 
for i in range(h):
    
    N_J5_d9[i] = N_G(t5,5,1,IC9_J5[i])
    

for j in range(len(t5)):
    sumN_J5=0
    for i in range(h):
        
        sumN_J5 = sumN_J5 + N_J5_d9[i,j] 
        
    mN_J5_d9[j] = sumN_J5/h    
    
#%%

plt.title(r'J=0, $\alpha=1$',fontsize=16)    
plt.ylabel(r'$\phi$',fontsize=14)
plt.xlabel('t',fontsize=14)
plt.grid()    
plt.plot(t0,mN_J0_d9,lw=2,label='Rel: $<n_0>=9$')
plt.plot(t0,N_J0_f9,ls='--',lw=2,label='Rel: $n_0=9$')
plt.plot(t0,N_J0_0,label='Res')

plt.legend()

plt.show()

plt.title(r'J=3, $\alpha=1$',fontsize=16)    
plt.ylabel(r'$\phi$',fontsize=14)
plt.xlabel('t',fontsize=14)
plt.grid()  
plt.plot(t3,mN_J3_d9,lw=2,label='Rel: $<n_0>=9$')
plt.plot(t3,N_J3_f9,ls='--',lw=2,label='Rel: $n_0=9$')
plt.plot(t3,N_J3_0,label='Res')

plt.legend()

plt.show()

plt.title(r'J=5, $\alpha=1$',fontsize=16)    
plt.ylabel(r'$\phi$',fontsize=14)
plt.xlabel('t',fontsize=14)
plt.grid()  
plt.plot(t5,mN_J5_d9,lw=2,label='Rel: $<n_0>=9$')
plt.plot(t5,N_J5_f9,ls='--',lw=2,label='Rel: $n_0=9$')
plt.plot(t5,N_J5_0,label='Res')

plt.legend()

plt.show()

#%%

def Nl(t,n0,tau,nss):
    
    return nss + (n0 - nss)*np.exp(-t/tau)

poRes_J0, pcRes_J0 = curve_fit(Nl,t0,N_J0_0)
poRel_J0_f9, pcRel_J0_f9 = curve_fit(Nl,t0,N_J0_f9)
poRel_J0_d9, pcRel_J0_d9 = curve_fit(Nl,t0,mN_J0_d9)

poRes_J3, pcRes_J3 = curve_fit(Nl,t3,N_J3_0)
poRel_J3_f9, pcRel_J3_f9 = curve_fit(Nl,t3,N_J3_f9)
poRel_J3_d9, pcRel_J3_d9 = curve_fit(Nl,t3,mN_J3_d9)

poRes_J5, pcRes_J5 = curve_fit(Nl,t5,N_J5_0)
poRel_J5_f9, pcRel_J5_f9 = curve_fit(Nl,t5,N_J5_f9)
poRel_J5_d9, pcRel_J5_d9 = curve_fit(Nl,t5,mN_J5_d9)


tRes[0] = poRes_J0[1]
tRes[1] = poRes_J3[1]
tRes[2] = poRes_J5[1]

tRel_f9[0]=poRel_J0_f9[1]
tRel_f9[1]=poRel_J3_f9[1]
tRel_f9[2]=poRel_J5_f9[1]

tRel_d9[0]=poRel_J0_d9[1]
tRel_d9[1]=poRel_J3_d9[1]
tRel_d9[2]=poRel_J5_d9[1]

Dt_f9 = (tRel_f9-tRes)/(tRel_f9+tRes)
Dt_d9 = (tRel_d9-tRes)/(tRel_d9+tRes)

#%%

J=([0,3,5])

plt.title(r'Langmuir fit',fontsize=16) 
plt.xlabel('J',fontsize=14)
plt.ylabel(r'$\Delta$',fontsize=14)
plt.grid()

plt.plot(J,Dt_f9,ls='',marker='o',label='$n_0=9$')
plt.plot(J,Dt_d9, ls='',marker='^',label='$<n_0>=9$')

plt.legend()
plt.show()

#%%

N_J0_f13 = N_G(t0,0,1,1)

N_J3_f13 = N_G(t3,3,1,1)

N_J5_f13 = N_G(t5,5,1,1)

poRel_J0_f13, pcRel_J0_f13 = curve_fit(Nl,t0,N_J0_f13)
poRel_J3_f13, pcRel_J3_f13 = curve_fit(Nl,t3,N_J3_f13)
poRel_J5_f13, pcRel_J5_f13 = curve_fit(Nl,t5,N_J5_f13)

tRel_f13[0]=poRel_J0_f13[1]
tRel_f13[1]=poRel_J3_f13[1]
tRel_f13[2]=poRel_J5_f13[1]

Dt_f13 = (tRel_f13-tRes)/(tRel_f13+tRes)

#%%

plt.title(r'Langmuir fit',fontsize=16) 
plt.xlabel('J',fontsize=14)
plt.ylabel(r'$\Delta$',fontsize=14)
plt.grid()

plt.plot(J,Dt_f9,ls='',marker='o',label=r'$n_0=9$')
plt.plot(J,Dt_f13,ls='',marker='o',alpha=0.5,label=r'$n_0=13$')

plt.legend()