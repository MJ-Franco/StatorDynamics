#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 13:19:58 2022

@author: luke
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
plt.rcParams["figure.figsize"] = [10.0,8.0]

nmax=13

l=20000

#%%

files6_J0=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Initial_conditions_J=0.00_N=65.dat']
data6_J0=[]
files6_J3=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Initial_conditions_J=3.00_N=65.dat']
data6_J3=[]
files6_J5=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Initial_conditions_J=5.00_N=65.dat']
data6_J5=[]

for data_file in files6_J0:
    data6_J0.append(np.loadtxt(data_file))
for data_file in files6_J3:
    data6_J3.append(np.loadtxt(data_file))
for data_file in files6_J5:
    data6_J5.append(np.loadtxt(data_file))
    
IC6_J0=data6_J0[0]
IC6_J3=data6_J3[0]
IC6_J5=data6_J5[0]

files9_J0=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Initial_conditions_J=0.00_N=9.dat']
data9_J0=[]
files9_J3=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Initial_conditions_J=3.00_N=9.dat']
data9_J3=[]
files9_J5=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Initial_conditions_J=5.00_N=9.dat']
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

files13_J0=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Initial_conditions_J=0.00_N=13.dat']
data13_J0=[]
files13_J3=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Initial_conditions_J=3.00_N=13.dat']
data13_J3=[]
files13_J5=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Initial_conditions_J=5.00_N=13.dat']
data13_J5=[]

for data_file in files13_J0:
    data13_J0.append(np.loadtxt(data_file))
for data_file in files13_J3:
    data13_J3.append(np.loadtxt(data_file))
for data_file in files13_J5:
    data13_J5.append(np.loadtxt(data_file))
    
IC13_J0=data13_J0[0]
IC13_J3=data13_J3[0]
IC13_J5=data13_J5[0]

#%%

plt.xlabel('N',fontsize=16)
plt.ylabel('Frecuency',fontsize=16)
plt.title('$<N> = 6, N_{max}=13$ ', fontsize = 16)
plt.grid(alpha=0.7)
plt.hist(IC6_J0*nmax,ls='solid',lw=2,bins=11, edgecolor='orange',hatch='.',fc='None',alpha=0.4,label='J=0')
plt.hist(IC6_J3*nmax,ls='solid',lw=2,bins=14, edgecolor='green',hatch='-',fc='None',alpha=0.4,label='J=3')
plt.hist(IC6_J5*nmax,ls='solid',lw=2,bins=14, edgecolor='royalblue',hatch='/',fc='None',alpha=0.5,label='J=5')
plt.vlines(6.5,0,3500,color='red',alpha=0.7,ls='dashed')
plt.legend(prop={'size': 17}, loc='upper left')

plt.show()

#%%

plt.xlabel('N',fontsize=16)
plt.ylabel('f',fontsize=16)
plt.title('$<n_0> = 9,N_{max}=13$ ')
plt.grid(alpha=0.3)
plt.hist(IC9_J0*nmax, density=True,ls='solid',lw=2,bins=12, edgecolor='orange',fc='None',alpha=0.3,label='J=0')
plt.hist(IC9_J3*nmax, density=True,ls='solid',lw=2,bins=14, edgecolor='green',fc='None',alpha=0.3,label='J=3')
plt.hist(IC9_J5*nmax, density=True,ls='solid',lw=2,bins=14, edgecolor='royalblue',fc='None',alpha=0.3,label='J=5')

plt.vlines(9,0,0.55,color='red',alpha=0.3,ls='dashed')
plt.legend(prop={'size': 12})

plt.show()

plt.xlabel('N',fontsize=16)
plt.ylabel('f',fontsize=16)
plt.title('$<n_0> = 13,N_{max}=13$ ')
plt.grid(alpha=0.3)
plt.hist(IC13_J0*nmax, density=True,ls='solid',lw=2, bins=14,edgecolor='orange',fc='None',alpha=0.3,label='J=0')
plt.hist(IC13_J3*nmax, density=True,ls='solid',lw=2, bins=14,edgecolor='green',fc='None',alpha=0.3,label='J=3')
plt.hist(IC13_J5*nmax, density=True,ls='solid',lw=2,bins=14,edgecolor='royalblue',fc='None',alpha=0.3,label='J=5')

plt.vlines(13,0,11,color='red',alpha=0.3,ls='dashed')
plt.legend(prop={'size': 12})

#%%

files4_J0_nm20=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Initial_conditions_J=0.00_N=4_nm=20.dat']
data4_J0_nm20=[]
files4_J3_nm20=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Initial_conditions_J=3.00_N=4_nm=20.dat']
data4_J3_nm20=[]
files4_J5_nm20=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Initial_conditions_J=5.00_N=4_nm=20.dat']
data4_J5_nm20=[]

for data_file in files4_J0_nm20:
    data4_J0_nm20.append(np.loadtxt(data_file))
for data_file in files4_J3_nm20:
    data4_J3_nm20.append(np.loadtxt(data_file))
for data_file in files4_J5_nm20:
    data4_J5_nm20.append(np.loadtxt(data_file))
    
IC4_J0_nm20=data4_J0_nm20[0]
IC4_J3_nm20=data4_J3_nm20[0]
IC4_J5_nm20=data4_J5_nm20[0]

plt.xlabel('N',fontsize=16)
plt.ylabel('f',fontsize=16)
plt.title('$N = 4, N_{max}=20$ ')
plt.grid(alpha=0.3)
plt.hist(IC4_J0_nm20*20, density=True,ls='solid',lw=2,bins=12, edgecolor='orange',fc='None',alpha=0.3,label='J=0')
plt.hist(IC4_J3_nm20*20, density=True,ls='solid',lw=2,bins=21, edgecolor='green',fc='None',alpha=0.3,label='J=3')
plt.hist(IC4_J5_nm20*20, density=True,ls='solid',lw=2,bins=21, edgecolor='royalblue',fc='None',alpha=0.3,label='J=5')
plt.vlines(4,0,0.35,color='red',alpha=0.3,ls='dashed')
plt.legend(prop={'size': 12})

plt.show()

#%%

files10_J0_nm20=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Initial_conditions_J=0.00_N=10_nm=20.dat']
data10_J0_nm20=[]
files10_J3_nm20=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Initial_conditions_J=3.00_N=10_nm=20.dat']
data10_J3_nm20=[]
files10_J5_nm20=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Initial_conditions_J=5.00_N=10_nm=20.dat']
data10_J5_nm20=[]

for data_file in files10_J0_nm20:
    data10_J0_nm20.append(np.loadtxt(data_file))
for data_file in files10_J3_nm20:
    data10_J3_nm20.append(np.loadtxt(data_file))
for data_file in files10_J5_nm20:
    data10_J5_nm20.append(np.loadtxt(data_file))
    
IC10_J0_nm20=data10_J0_nm20[0]
IC10_J3_nm20=data10_J3_nm20[0]
IC10_J5_nm20=data10_J5_nm20[0]

plt.xlabel('N',fontsize=16)
plt.ylabel('f',fontsize=16)
plt.title('$N = 10, N_{max}=20$ ')
plt.grid(alpha=0.3)
plt.hist(IC10_J0_nm20*20, density=True,ls='solid',lw=2,bins=15, edgecolor='orange',fc='None',alpha=0.3,label='J=0')
plt.hist(IC10_J3_nm20*20, density=True,ls='solid',lw=2,bins=20, edgecolor='green',fc='None',alpha=0.3,label='J=3')
plt.hist(IC10_J5_nm20*20, density=True,ls='solid',lw=2,bins=20, edgecolor='royalblue',fc='None',alpha=0.3,label='J=5')
plt.vlines(10,0,0.25,color='red',alpha=0.3,ls='dashed')
plt.legend(prop={'size': 12})

plt.show()

#%%

files16_J0_nm20=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Initial_conditions_J=0.00_N=14_nm=20.dat']
data16_J0_nm20=[]
files16_J3_nm20=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Initial_conditions_J=3.00_N=14_nm=20.dat']
data16_J3_nm20=[]
files16_J5_nm20=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Initial_conditions_J=5.00_N=14_nm=20.dat']
data16_J5_nm20=[]

for data_file in files16_J0_nm20:
    data16_J0_nm20.append(np.loadtxt(data_file))
for data_file in files16_J3_nm20:
    data16_J3_nm20.append(np.loadtxt(data_file))
for data_file in files16_J5_nm20:
    data16_J5_nm20.append(np.loadtxt(data_file))
    
IC16_J0_nm20=data16_J0_nm20[0]
IC16_J3_nm20=data16_J3_nm20[0]
IC16_J5_nm20=data16_J5_nm20[0]

plt.xlabel('N',fontsize=16)
plt.ylabel('f',fontsize=16)
plt.title('$N = 16, N_{max}=20$ ')
plt.grid(alpha=0.3)
plt.hist(IC16_J0_nm20*20, density=True,ls='solid',lw=2,bins=12, edgecolor='orange',fc='None',alpha=0.3,label='J=0')
plt.hist(IC16_J3_nm20*20, density=True,ls='solid',lw=2,bins=20, edgecolor='green',fc='None',alpha=0.3,label='J=3')
plt.hist(IC16_J5_nm20*20, density=True,ls='solid',lw=2,bins=20, edgecolor='royalblue',fc='None',alpha=0.3,label='J=5')
plt.vlines(16,0,0.35,color='red',alpha=0.3,ls='dashed')
plt.legend(prop={'size': 12})

plt.show()

#%%


files20_J0_nm40=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Initial_conditions_J=0.00_N=20_nm=40.dat']
data20_J0_nm40=[]
files20_J3_nm40=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Initial_conditions_J=3.00_N=20_nm=40.dat']
data20_J3_nm40=[]
files20_J5_nm40=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Initial_conditions_J=5.00_N=20_nm=40.dat']
data20_J5_nm40=[]

for data_file in files20_J0_nm40:
    data20_J0_nm40.append(np.loadtxt(data_file))
for data_file in files20_J3_nm40:
    data20_J3_nm40.append(np.loadtxt(data_file))
for data_file in files20_J5_nm40:
    data20_J5_nm40.append(np.loadtxt(data_file))
    
IC20_J0_nm40=data20_J0_nm40[0]
IC20_J3_nm40=data20_J3_nm40[0]
IC20_J5_nm40=data20_J5_nm40[0]

plt.xlabel('N',fontsize=16)
plt.ylabel('f',fontsize=16)
plt.title('$N = 20, N_{max}=40$ ')
plt.grid(alpha=0.3)
plt.hist(IC20_J0_nm40*40, density=True,ls='solid',lw=2,bins=40, edgecolor='orange',fc='None',alpha=0.3,label='J=0')
plt.hist(IC20_J3_nm40*40, density=True,ls='solid',lw=2,bins=40, edgecolor='green',fc='None',alpha=0.3,label='J=3')
plt.hist(IC20_J5_nm40*40, density=True,ls='solid',lw=2,bins=40, edgecolor='royalblue',fc='None',alpha=0.3,label='J=5')
plt.vlines(20,0,0.2,color='red',alpha=0.3,ls='dashed')
plt.legend(prop={'size': 12})

plt.show()

