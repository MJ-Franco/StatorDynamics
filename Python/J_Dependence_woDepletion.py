#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 16:53:41 2020

@author: mariajose
"""
import numpy as np
import matplotlib.pyplot as plt


#%% File for representing ntot vs Nss

files1=['/home/mariajose/Escritorio/Simulations/Without depletion/Dependence on J/J_Dependence_mu=  0.dat'] 
data1=[]
files2=['/home/mariajose/Escritorio/Simulations/Without depletion/Dependence on J/J_Dependence_mu= -1.dat'] 
data2=[]
files3=['/home/mariajose/Escritorio/Simulations/Without depletion/Dependence on J/J_Dependence_mu= -2.dat'] 
data3=[]
files4=['/home/mariajose/Escritorio/Simulations/Without depletion/Dependence on J/J_Dependence_mu= -3.dat'] 
data4=[]

for data_file in files1:
    data1.append(np.loadtxt(data_file))
    
for data_file in files2:
    data2.append(np.loadtxt(data_file))
    
for data_file in files3:
    data3.append(np.loadtxt(data_file))

for data_file in files4:
    data4.append(np.loadtxt(data_file))

mu_0=data1[0]
mu_1=data2[0]
mu_2=data3[0]
mu_3=data4[0]

#%%Plot

plt.title('No depletion')
plt.xlabel('$J$')
plt.ylabel('$<\phi>$')
plt.plot(mu_0[:,0],mu_0[:,1],label='$\mu$=0',marker='o')
plt.plot(mu_1[:,0],mu_1[:,1],label='$\mu$=-1',marker='o')
plt.plot(mu_2[:,0],mu_2[:,1],label='$\mu$=-2',marker='o')
plt.plot(mu_3[:,0],mu_3[:,1],label='$\mu$=-3',marker='o')
#plt.hlines(10.93,0, 200,alpha=0.2,linestyles='dashed')
#plt.hlines(6.03,0,200,alpha=0.2,linestyles='dashed', label='No depletion')
#plt.hlines(3.23, 0, 200,alpha=0.2,linestyles='dashed')


plt.legend(loc='upper right',bbox_to_anchor=(1, 0.4))


plt.show()

#%%Functions necessary for occupation

def X(J,mu):
    return 0.5*(J+mu)

def B(J):
    return np.exp(-J)

def lp(J,mu):
    return np.exp(X(J,mu))*(np.cosh(X(J,mu))+np.sqrt((np.sinh(X(J,mu)))**2+B(J)))

def lm(J,mu):
    return np.exp(X(J,mu))*(np.cosh(X(J,mu))-np.sqrt((np.sinh(X(J,mu)))**2+B(J)))

def vm(J,mu):
    return np.exp(-0.5*mu)*(lm(J,mu)-1)

#%% Mean occupation

M=12.0

#Exact
def N(J,mu):
    return ((lp(J,mu)**M+lm(J,mu)**M*vm(J,mu)**2)/((lp(J,mu)**M + lm(J,mu)**M)*(1+vm(J,mu)**2)))

#Thermo limit
def Nth(J,mu):
    return 1/(1+(vm(J,mu))**2)


#%% Exact theoretical plot

Js=np.linspace(1, 10,10)

plt.title('Exact solution')
plt.xlabel('J / $k_BT$')
plt.ylabel('$N_{ss}$')
plt.plot(Js,N(Js,0),marker='o',label='$\mu$=0')
plt.plot(Js,N(Js,-1),marker='o',label='$\mu$=-1')
plt.plot(Js,N(Js,-2),marker='o',label='$\mu$=-2')
plt.plot(Js,N(Js,-3),marker='o',label='$\mu$=-3')

plt.legend(loc='upper right',bbox_to_anchor=(1, 0.4))

plt.show()

#%% Theoretical plot thermo limit

Js=np.linspace(1, 10,10)

plt.title('Thermodynamic limit')
plt.xlabel('J / $k_BT$')
plt.ylabel('$N_{ss}$')
plt.plot(Js,Nth(Js,0),marker='o',label='$\mu$=0')
plt.plot(Js,Nth(Js,-1),marker='o',label='$\mu$=-1')
plt.plot(Js,Nth(Js,-2),marker='o',label='$\mu$=-2')
plt.plot(Js,Nth(Js,-3),marker='o',label='$\mu$=-3')

plt.legend(loc='upper right',bbox_to_anchor=(1, 0.4))

plt.show()

#%%

plt.title('Exact solution / Stall ($n_0=10$)')
plt.xlabel("J / $K_BT$")
plt.ylabel("$<\phi>$")
plt.plot(Js,N(Js,0),marker='o',label='$\mu$=0',color='blue')
plt.errorbar(mu_0[:,0],mu_0[:,1],mu_0[:,2],marker='x',ls='dashed',color='blue')
plt.plot(Js,N(Js,-1),marker='o',label='$\mu$=-1',color='orange')
plt.errorbar(mu_1[:,0],mu_1[:,1],mu_1[:,2],marker='x',ls='dashed',color='orange')
plt.plot(Js,N(Js,-2),marker='o',label='$\mu$=-2',color='green')
plt.errorbar(mu_2[:,0],mu_2[:,1],mu_2[:,2],marker='x',ls='dashed',color='green')
plt.plot(Js,N(Js,-3),marker='o',label='$\mu$=-3',color='red')
plt.errorbar(mu_3[:,0],mu_3[:,1],mu_3[:,2],marker='x',ls='dashed',color='red')

plt.legend(loc='upper right',bbox_to_anchor=(1, 0.4))

plt.show()

#%% Comparison between both thermolimit

plt.title('Thermo limit')
plt.xlabel("J / $K_BT$")
plt.ylabel("$N_{ss}$")
plt.plot(Js,Nth(Js,0),marker='o',label='$\mu$=0',color='blue')
plt.errorbar(mu_0[:,0],mu_0[:,1],mu_0[:,2],marker='x',ls='dashed',color='blue')
plt.plot(Js,Nth(Js,-1),marker='o',label='$\mu$=-1',color='orange')
plt.errorbar(mu_1[:,0],mu_1[:,1],mu_1[:,2],marker='x',ls='dashed',color='orange')
plt.plot(Js,Nth(Js,-2),marker='o',label='$\mu$=-2',color='green')
plt.errorbar(mu_2[:,0],mu_2[:,1],mu_2[:,2],marker='x',ls='dashed',color='green')
plt.plot(Js,Nth(Js,-3),marker='o',label='$\mu$=-3',color='red')
plt.errorbar(mu_3[:,0],mu_3[:,1],mu_3[:,2],marker='x',ls='dashed',color='red')

plt.legend(loc='upper right',bbox_to_anchor=(1, 0.4))

plt.show()