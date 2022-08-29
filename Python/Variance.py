#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 17:50:56 2020

@author: mariajose
"""

import numpy as np
import matplotlib.pyplot as plt

#%% Data

files1=['/home/mariajose/Escritorio/Simulations/Without depletion/Variance vs Nmax/Var_Dependence_J=  0mu=  0.dat'] 
data1=[]
files2=['/home/mariajose/Escritorio/Simulations/Without depletion/Variance vs Nmax/Var_Dependence_J=  1mu= -1.dat'] 
data2=[]
files3=['/home/mariajose/Escritorio/Simulations/Without depletion/Variance vs Nmax/Var_Dependence_J=  2mu= -2.dat']
data3=[] 


for data_file in files1:
    data1.append(np.loadtxt(data_file))
    
for data_file in files2:
    data2.append(np.loadtxt(data_file))
    
for data_file in files3:
    data3.append(np.loadtxt(data_file))


J_0=data1[0]
J_1=data2[0]
J_2=data3[0]

#%%Plot variance

plt.title('Variance (with depletion, $n_{tot}=11$)')
plt.xlabel('$N_{max}$')
plt.ylabel('$<\phi^2>-<\phi>^2$')
#plt.plot(J_0[:,0]*10,J_0[:,2],label='J=0, $\mu$=0',marker='o')
plt.errorbar(J_0[:,0]*10,J_0[:,2],J_0[:,3],label='J=0, $\mu_0$=0',marker='x')
plt.errorbar(J_1[:,0]*10,J_1[:,2],J_1[:,3],label='J=1, $\mu_0$=-1',marker='x')
plt.errorbar(J_2[:,0]*10,J_2[:,2],J_2[:,3],label='J=2, $\mu_0$=-2',marker='x')
#plt.plot(J_1[:,0]*10,J_1[:,2],label='J=1, $\mu$=-1',marker='o')
#plt.plot(J_2[:,0]*10,J_2[:,2],label='J=2, $\mu$=-2',marker='o')
#plt.hlines(10.93,0, 200,alpha=0.2,linestyles='dashed')
#plt.hlines(6.03,0,200,alpha=0.2,linestyles='dashed', label='No depletion')
#plt.hlines(3.23, 0, 200,alpha=0.2,linestyles='dashed')


plt.legend(loc='upper right')


plt.show()

#%%Plot Nss

plt.title('Steady state (with depletion, $n_{tot}=11$)')
plt.xlabel('$N_{max}$')
plt.ylabel('$<\phi>$')
#plt.plot(J_0[:,0]*10,J_0[:,2],label='J=0, $\mu$=0',marker='o')
plt.errorbar(J_0[:,0]*10,J_0[:,1],J_0[:,2],label='J=0, $\mu_0$=0',marker='x')
plt.errorbar(J_1[:,0]*10,J_1[:,1],J_1[:,2],label='J=1, $\mu_0$=-1',marker='x')
plt.errorbar(J_2[:,0]*10,J_2[:,1],J_2[:,2],label='J=2, $\mu_0$=-2',marker='x')
#plt.plot(J_1[:,0]*10,J_1[:,2],label='J=1, $\mu$=-1',marker='o')
#plt.plot(J_2[:,0]*10,J_2[:,2],label='J=2, $\mu$=-2',marker='o')
#plt.hlines(10.93,0, 200,alpha=0.2,linestyles='dashed')
#plt.hlines(6.03,0,200,alpha=0.2,linestyles='dashed', label='No depletion')
#plt.hlines(3.23, 0, 200,alpha=0.2,linestyles='dashed')


plt.legend(loc='upper right')


plt.show()


