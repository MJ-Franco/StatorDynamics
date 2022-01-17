#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 13:18:36 2021

@author: mariajose
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

#%% Useful variables and functions

nmax=13
l=5000
#Function used for fits
def N(x,tau,nss):
    return nss + (n0-nss)*np.exp(-x/tau) 

#%%

files1=['/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/RelaxTimes_vs_Nmax/RelaxTime_vs_Nmax= 10_Res.dat']
data1=[]
files2=['/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/RelaxTimes_vs_Nmax/RelaxTime_vs_Nmax= 20_Res.dat']
data2=[]
files3=['/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/RelaxTimes_vs_Nmax/RelaxTime_vs_Nmax= 30_Res.dat']
data3=[]
files4=['/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/RelaxTimes_vs_Nmax/RelaxTime_vs_Nmax= 40_Res.dat']
data4=[]
files5=['/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/RelaxTimes_vs_Nmax/RelaxTime_vs_Nmax= 50_Res.dat']
data5=[]
files6=['/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/RelaxTimes_vs_Nmax/RelaxTime_vs_Nmax= 60_Res.dat']
data6=[]
files7=['/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/RelaxTimes_vs_Nmax/RelaxTime_vs_Nmax= 70_Res.dat']
data7=[]
files8=['/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/RelaxTimes_vs_Nmax/RelaxTime_vs_Nmax= 80_Res.dat']
data8=[]
files9=['/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/RelaxTimes_vs_Nmax/RelaxTime_vs_Nmax= 80_Res.dat']
data9=[]
files10=['/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/RelaxTimes_vs_Nmax/RelaxTime_vs_Nmax=100_Res.dat']
data10=[]


for data_file in files1:
    data1.append(np.loadtxt(data_file))

for data_file in files2:
    data2.append(np.loadtxt(data_file))
    
for data_file in files3:
    data3.append(np.loadtxt(data_file)) 
    
for data_file in files4:
    data4.append(np.loadtxt(data_file))
    
for data_file in files5:
    data5.append(np.loadtxt(data_file))

for data_file in files6:
    data6.append(np.loadtxt(data_file))

for data_file in files7:
    data7.append(np.loadtxt(data_file))
    
for data_file in files8:
    data8.append(np.loadtxt(data_file))

for data_file in files9:
    data9.append(np.loadtxt(data_file))
    
for data_file in files10:
    data10.append(np.loadtxt(data_file))
    

Nm10=data1[0]
Nm20=data2[0]
Nm30=data3[0]
Nm40=data4[0]
Nm50=data5[0]
Nm60=data6[0]
Nm70=data7[0]
Nm80=data8[0]
Nm90=data9[0]
Nm100=data10[0]

#%%

n0=0

poRes10, pcRes10 = curve_fit(N,Nm10[0:l,0],Nm10[0:l,1]*nmax)
tau_Res10=poRes10[0]

poRes20, pcRes20 = curve_fit(N,Nm20[0:l,0],Nm20[0:l,1]*nmax)
tau_Res20=poRes20[0]

poRes30, pcRes30 = curve_fit(N,Nm30[0:l,0],Nm30[0:l,1]*nmax)
tau_Res30=poRes30[0]

poRes40, pcRes40 = curve_fit(N,Nm40[0:l,0],Nm40[0:l,1]*nmax)
tau_Res40=poRes40[0]

poRes50, pcRes50 = curve_fit(N,Nm50[0:l,0],Nm50[0:l,1]*nmax)
tau_Res50=poRes50[0]

poRes60, pcRes60 = curve_fit(N,Nm60[0:l,0],Nm60[0:l,1]*nmax)
tau_Res60=poRes60[0]

poRes70, pcRes70 = curve_fit(N,Nm70[0:l,0],Nm70[0:l,1]*nmax)
tau_Res70=poRes70[0]

poRes80, pcRes80 = curve_fit(N,Nm80[0:l,0],Nm80[0:l,1]*nmax)
tau_Res80=poRes80[0]

poRes90, pcRes90 = curve_fit(N,Nm90[0:l,0],Nm90[0:l,1]*nmax)
tau_Res90=poRes90[0]

poRes100, pcRes100 = curve_fit(N,Nm100[0:l,0],Nm100[0:l,1]*nmax)
tau_Res100=poRes100[0]


#%%

Nm = np.arange(10,110,10) 

tau_Res = np.array([tau_Res10,tau_Res20,tau_Res30,tau_Res40,tau_Res50,tau_Res60,tau_Res70,tau_Res80,tau_Res90,tau_Res100])

#%%

def r(p,a,b):
    return a*p + b

por, pcr = curve_fit(r,Nm,tau_Res)
tr = r(Nm,*por)

plt.xlabel('$N_{max}$')
plt.ylabel('Relaxation time')
plt.plot(Nm,tau_Res,marker='o',linestyle='')
plt.plot(Nm,tr)

slope = por[0]
constant= por[1]

print(slope,constant)

plt.grid()