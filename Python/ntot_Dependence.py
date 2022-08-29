#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 16:30:21 2020

@author: mariajose
"""

import numpy as np
import matplotlib.pyplot as plt

l=200
nmax=13

#%% File for representing ntot vs Nss

files1=['/media/mfranco/Elements/Simulations/With_depletion/Nss_vs_ntot/ntot_Dependence_J=0.dat'] 
data1=[]
files2=['/media/mfranco/Elements/Simulations/With_depletion/Nss_vs_ntot/ntot_Dependence_J=1.dat'] 
data2=[]
files3=['/media/mfranco/Elements/Simulations/With_depletion/Nss_vs_ntot/ntot_Dependence_J=2.dat'] 
data3=[]
files4=['/media/mfranco/Elements/Simulations/With_depletion/Nss_vs_ntot/ntot_Dependence_J=3.dat'] 
data4=[]

for data_file in files1:
    data1.append(np.loadtxt(data_file))
    
for data_file in files2:
    data2.append(np.loadtxt(data_file))
    
for data_file in files3:
    data3.append(np.loadtxt(data_file))

for data_file in files4:
    data4.append(np.loadtxt(data_file))

J_0=data1[0]
J_1=data2[0]
J_2=data3[0]
J_3=data4[0]

#%% Plot

plt.title('$N=8$',fontsize=17)
plt.xlabel('$n_{tot}$',fontsize=14)
plt.ylabel('Steady state',fontsize=14)
plt.rcParams["figure.figsize"] = [8.0,6.0]

plt.errorbar(J_0[0:l,0],J_0[0:l,1]*nmax,J_0[0:l,2],label='J=0',marker='x',capsize=4)
plt.errorbar(J_1[0:l,0],J_1[0:l,1]*nmax,J_1[0:l,2],label='J=1',marker='x',capsize=4)
plt.errorbar(J_2[0:l,0],J_2[0:l,1]*nmax,J_2[0:l,2],label='J=2',marker='x',capsize=4)
plt.errorbar(J_3[0:l,0],J_3[0:l,1]*nmax,J_3[0:l,2],label='J=3',marker='x',capsize=4)
#plt.hlines(2.09/nmax,0, 200,alpha=0.5,linestyles='dashed')
#plt.hlines(11.51/nmax,0, 200,alpha=0.5,linestyles='dashed')
plt.hlines(8,0,400,alpha=0.5,linestyles='dashed', label='No depletion')
#plt.hlines(1.43/nmax, 0, 200,alpha=0.5,linestyles='dashed')


plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.,prop={'size': 15})


plt.show()

#%% Log plot
plt.xscale('log')
plt.yscale('log')

plt.title('$\mu=-1$')
plt.xlabel('log($n_{tot}$)')
plt.ylabel('log($<\phi>$)')
plt.errorbar(J_0[:,0],J_0[:,1],J_0[:,2],label='J=0',marker='o')
plt.errorbar(J_1[:,0],J_1[:,1],J_0[:,2],label='J=1',marker='o')
plt.errorbar(J_2[:,0],J_2[:,1],J_2[:,2],label='J=2',marker='o')

plt.legend()

plt.show()

#%% Semilog plot

plt.title('$\mu=-1$')
plt.xlabel('log($n_{tot}$)')
plt.ylabel('$N_{ss}$')
plt.semilogy(J_0[:,0],J_0[:,1],label='J=0',marker='o')
plt.semilogy(J_1[:,0],J_1[:,1],label='J=1',marker='o')
plt.semilogy(J_2[:,0],J_2[:,1],label='J=2',marker='o')

plt.legend()

plt.show()

