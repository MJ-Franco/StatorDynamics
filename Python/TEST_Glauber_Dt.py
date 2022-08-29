#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 16:09:38 2022

@author: mfranco
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
plt.rcParams["figure.figsize"] = [8.0,6.0]

def N(x,tau,nss):
    return nss + (n0-nss)*np.exp(-x/tau)



nmax=13


tRes4=np.zeros(6)
tRes6=np.zeros(6)
tRes8=np.zeros(6)
errtRes4=np.zeros(6)
errtRes6=np.zeros(6)
errtRes8=np.zeros(6)

tStall4d=np.zeros(6)
tStall8d=np.zeros(6)
errtStall4d=np.zeros(6)
errtStall8d=np.zeros(6)

tStall4f_13=np.zeros(6)
tStall8f_13=np.zeros(6)
tStall6f_13=np.zeros(6)
errtStall4f_13=np.zeros(6)
errtStall8f_13=np.zeros(6)
errtStall6f_13=np.zeros(6)

tStall4f_6=np.zeros(6)
tStall8f_10=np.zeros(6)
tStall6f_8=np.zeros(6)
errtStall4f_6=np.zeros(6)
errtStall8f_10=np.zeros(6)
errtStall6f_8=np.zeros(6)

tStall8f_9=np.zeros(6)
tStall6d_9=np.zeros(6)
errtStall8f_9=np.zeros(6)
errtStall6d_9=np.zeros(6)

tStall4d_13=np.zeros(6)
tStall8d_13=np.zeros(6)
tStall6d_13=np.zeros(6)
errtStall4d_13=np.zeros(6)
errtStall8d_13=np.zeros(6)
errtStall6d_13=np.zeros(6)


#%%RELAXATION TIME DATA WITH DISTRIBUTION

FiletS4D=['D:/Simulations/Without _depletion/Glauber/RelaxTimes/StallTime_4_n06_dist.dat']
DatatS4D=[]
FiletS8D=['D:/Simulations/Without _depletion/Glauber/RelaxTimes/StallTime_8_n09_dist.dat']
DatatS8D=[]

FiletR4=['D:/Simulations/Without _depletion/Glauber/RelaxTimes/ResTime_4_dist.dat']
DatatR4=[]
FiletR8=['D:/Simulations/Without _depletion/Glauber/RelaxTimes/ResTime_8_dist.dat']
DatatR8=[]

for data_file in FiletS4D:
    DatatS4D.append(np.loadtxt(data_file))
for data_file in FiletS8D:
    DatatS8D.append(np.loadtxt(data_file))
    
for data_file in FiletR4:
    DatatR4.append(np.loadtxt(data_file))
for data_file in FiletR8:
    DatatR8.append(np.loadtxt(data_file))



tRes4_old=DatatR4[0]
tRes8_old=DatatR8[0]
tStall4d_old=DatatS4D[0]
tStall8d_old=DatatS8D[0]

Dt4d_old = (tStall4d_old[:,1] - tRes4_old[:,1])/(tStall4d_old[:,1] + tRes4_old[:,1])
Dt8d_old = (tStall8d_old[:,1] - tRes8_old[:,1])/(tStall8d_old[:,1] + tRes8_old[:,1])


plt.xlabel('J',fontsize=16)
plt.ylabel('$\Delta$',fontsize=16)
plt.title('$N_{max}=13$ (old data)')
plt.grid()
plt.plot(tRes4_old[:,0],Dt4d_old,ls='',marker='o',label='<N>=4, $<n_0>=6$')
plt.plot(tRes8_old[:,0],Dt8d_old,ls='',marker='^',label='<N>=8, $<n_0>=9$')

plt.legend()



#%%RELAXATION TIME DATA WITHOUT DISTRIBUTION

FiletS4F=['D:/Simulations/Without _depletion/Glauber/RelaxTimes/StallTime_4_n06.dat']
DatatS4F=[]
FiletS8F=['D:/Simulations/Without _depletion/Glauber/RelaxTimes/StallTime_8_n09.dat']
DatatS8F=[]


for data_file in FiletS4F:
    DatatS4F.append(np.loadtxt(data_file))
for data_file in FiletS8F:
    DatatS8F.append(np.loadtxt(data_file))
    

tStall4f=DatatS4F[0]
tStall8f=DatatS8F[0]



Dt4f = (tStall4f[:,1] - tRes4_old[0:17,1])/(tStall4f[:,1] + tRes4_old[0:17,1])
Dt8f = (tStall8f[:,1] - tRes8_old[0:17,1])/(tStall8f[:,1] + tRes8_old[0:17,1])


plt.xlabel('J',fontsize=16)
plt.ylabel('$\Delta$',fontsize=16)
plt.title('$N_{max}=13$ (old data)')
plt.grid()
plt.plot(tRes4_old[0:17,0],Dt4f,ls='',marker='o',label='<N>=4, $n_0=6$')
plt.plot(tRes8_old[0:17,0],Dt8f,ls='',marker='^',label='<N>=8, $n_0=9$')

plt.legend()

#%% RECALCULATION OF RELAX TIMES

l=10000


FileR4_J0=['D:/Simulations/Without _depletion/Glauber/Traces_Res_J=0.00_N=4.dat']
DataR4_J0=[]
FileR4_J1=['D:/Simulations/Without _depletion/Glauber/Traces_Res_J=1.00_N=4.dat']
DataR4_J1=[]
FileR4_J2=['D:/Simulations/Without _depletion/Glauber/Traces_Res_J=2.00_N=4.dat']
DataR4_J2=[]
FileR4_J3=['D:/Simulations/Without _depletion/Glauber/Traces_Res_J=3.00_N=4.dat']
DataR4_J3=[]
FileR4_J4=['D:/Simulations/Without _depletion/Glauber/Traces_Res_J=4.00_N=4.dat']
DataR4_J4=[]
FileR4_J5=['D:/Simulations/Without _depletion/Glauber/Traces_Res_J=5.00_N=4.dat']
DataR4_J5=[]

FileS4d_J0=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=0.00_N=4.dat']
DataS4d_J0=[]
FileS4d_J1=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=1.00_N=4.dat']
DataS4d_J1=[]
FileS4d_J2=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=2.00_N=4.dat']
DataS4d_J2=[]
FileS4d_J3=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=3.00_N=4.dat']
DataS4d_J3=[]
FileS4d_J4=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=4.00_N=4.dat']
DataS4d_J4=[]
FileS4d_J5=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=5.00_N=4.dat']
DataS4d_J5=[]

for data_file in FileR4_J0:
    DataR4_J0.append(np.loadtxt(data_file))
for data_file in FileR4_J1:
    DataR4_J1.append(np.loadtxt(data_file))
for data_file in FileR4_J2:
    DataR4_J2.append(np.loadtxt(data_file))
for data_file in FileR4_J3:
    DataR4_J3.append(np.loadtxt(data_file))
for data_file in FileR4_J4:
    DataR4_J4.append(np.loadtxt(data_file))
for data_file in FileR4_J5:
    DataR4_J5.append(np.loadtxt(data_file))
    
for data_file in FileS4d_J0:
    DataS4d_J0.append(np.loadtxt(data_file))
for data_file in FileS4d_J1:
    DataS4d_J1.append(np.loadtxt(data_file))
for data_file in FileS4d_J2:
    DataS4d_J2.append(np.loadtxt(data_file))
for data_file in FileS4d_J3:
    DataS4d_J3.append(np.loadtxt(data_file))
for data_file in FileS4d_J4:
    DataS4d_J4.append(np.loadtxt(data_file))
for data_file in FileS4d_J5:
    DataS4d_J5.append(np.loadtxt(data_file))
    
Res4_J0=DataR4_J0[0]
Res4_J1=DataR4_J1[0]
Res4_J2=DataR4_J2[0]
Res4_J3=DataR4_J3[0]
Res4_J4=DataR4_J4[0]
Res4_J5=DataR4_J5[0]

Stall4d_J0=DataS4d_J0[0]
Stall4d_J1=DataS4d_J1[0]
Stall4d_J2=DataS4d_J2[0]
Stall4d_J3=DataS4d_J3[0]
Stall4d_J4=DataS4d_J4[0]
Stall4d_J5=DataS4d_J5[0]    



FileR8_J0=['D:/Simulations/Without _depletion/Glauber/Traces_Res_J=0.00_N=8.dat']
DataR8_J0=[]
FileR8_J1=['D:/Simulations/Without _depletion/Glauber/Traces_Res_J=1.00_N=8.dat']
DataR8_J1=[]
FileR8_J2=['D:/Simulations/Without _depletion/Glauber/Traces_Res_J=2.00_N=8.dat']
DataR8_J2=[]
FileR8_J3=['D:/Simulations/Without _depletion/Glauber/Traces_Res_J=3.00_N=8.dat']
DataR8_J3=[]
FileR8_J4=['D:/Simulations/Without _depletion/Glauber/Traces_Res_J=4.00_N=8.dat']
DataR8_J4=[]
FileR8_J5=['D:/Simulations/Without _depletion/Glauber/Traces_Res_J=5.00_N=8.dat']
DataR8_J5=[]

FileS8d_J0=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=0.00_N=8.dat']
DataS8d_J0=[]
FileS8d_J1=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=1.00_N=8.dat']
DataS8d_J1=[]
FileS8d_J2=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=2.00_N=8.dat']
DataS8d_J2=[]
FileS8d_J3=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=3.00_N=8.dat']
DataS8d_J3=[]
FileS8d_J4=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=4.00_N=8.dat']
DataS8d_J4=[]
FileS8d_J5=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=5.00_N=8.dat']
DataS8d_J5=[]

for data_file in FileR8_J0:
    DataR8_J0.append(np.loadtxt(data_file))
for data_file in FileR8_J1:
    DataR8_J1.append(np.loadtxt(data_file))
for data_file in FileR8_J2:
    DataR8_J2.append(np.loadtxt(data_file))
for data_file in FileR8_J3:
    DataR8_J3.append(np.loadtxt(data_file))
for data_file in FileR8_J4:
    DataR8_J4.append(np.loadtxt(data_file))
for data_file in FileR8_J5:
    DataR8_J5.append(np.loadtxt(data_file))    

for data_file in FileS8d_J0:
    DataS8d_J0.append(np.loadtxt(data_file))
for data_file in FileS8d_J1:
    DataS8d_J1.append(np.loadtxt(data_file))
for data_file in FileS8d_J2:
    DataS8d_J2.append(np.loadtxt(data_file))
for data_file in FileS8d_J3:
    DataS8d_J3.append(np.loadtxt(data_file))
for data_file in FileS8d_J4:
    DataS8d_J4.append(np.loadtxt(data_file))
for data_file in FileS8d_J5:
    DataS8d_J5.append(np.loadtxt(data_file))
    
Res8_J0=DataR8_J0[0]
Res8_J1=DataR8_J1[0]
Res8_J2=DataR8_J2[0]
Res8_J3=DataR8_J3[0]
Res8_J4=DataR8_J4[0]
Res8_J5=DataR8_J5[0]

Stall8d_J0=DataS8d_J0[0]
Stall8d_J1=DataS8d_J1[0]
Stall8d_J2=DataS8d_J2[0]
Stall8d_J3=DataS8d_J3[0]
Stall8d_J4=DataS8d_J4[0]
Stall8d_J5=DataS8d_J5[0]

#%%

plt.rcParams["figure.figsize"] = [10.0,8.0]

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
#plt.yaxis.label.set_size(16)

plt.title(r'$J = 1$, $<N>_\infty = 8$  ($N_{max}=13$)',fontsize=18)
plt.grid()
plt.ylabel(r'<N>',fontsize=16)
plt.xlabel('Monte Carlo step',fontsize=16)
#plt.xlim(0,310)
plt.ylim(0,10)
plt.plot(Res8_J1[0:300,0],Res8_J1[0:300,1]*nmax)
plt.plot(Stall8d_J1[0:300,0],Stall8d_J1[0:300,1]*nmax)

#%%



#%% RESURRECTION 4 and 8

n0=0

t0=np.arange(0,300.,1e-3)
t3=np.arange(0,2000.,1e-4)
t5=np.arange(0,8000.,-1e4)

poRes4_J0, pcRes4_J0 = curve_fit(N,Res4_J0[0:l,0],Res4_J0[0:l,1]*nmax)
errRes4_J0=np.sqrt(np.diag(pcRes4_J0))
modelRes4_J0=N(t0,*poRes4_J0)
poRes4_J1, pcRes4_J1 = curve_fit(N,Res4_J1[0:l,0],Res4_J1[0:l,1]*nmax)
errRes4_J1=np.sqrt(np.diag(pcRes4_J1))
poRes4_J2, pcRes4_J2 = curve_fit(N,Res4_J2[0:l,0],Res4_J2[0:l,1]*nmax)
errRes4_J2=np.sqrt(np.diag(pcRes4_J2))
poRes4_J3, pcRes4_J3 = curve_fit(N,Res4_J3[0:l,0],Res4_J3[0:l,1]*nmax)
errRes4_J3=np.sqrt(np.diag(pcRes4_J3))
modelRes4_J3=N(t3,*poRes4_J3)
poRes4_J4, pcRes4_J4 = curve_fit(N,Res4_J4[0:l,0],Res4_J4[0:l,1]*nmax)
errRes4_J4=np.sqrt(np.diag(pcRes4_J4))
poRes4_J5, pcRes4_J5 = curve_fit(N,Res4_J5[0:l,0],Res4_J5[0:l,1]*nmax)
errRes4_J5=np.sqrt(np.diag(pcRes4_J5))
modelRes4_J5=N(t5,*poRes4_J5)


poRes8_J0, pcRes8_J0 = curve_fit(N,Res8_J0[0:l,0],Res8_J0[0:l,1]*nmax)
errRes8_J0=np.sqrt(np.diag(pcRes8_J0))
modelRes8_J0=N(t0,*poRes8_J0)
poRes8_J1, pcRes8_J1 = curve_fit(N,Res8_J1[0:l,0],Res8_J1[0:l,1]*nmax)
errRes8_J1=np.sqrt(np.diag(pcRes8_J1))
poRes8_J2, pcRes8_J2 = curve_fit(N,Res8_J2[0:l,0],Res8_J2[0:l,1]*nmax)
errRes8_J2=np.sqrt(np.diag(pcRes8_J2))
poRes8_J3, pcRes8_J3 = curve_fit(N,Res8_J3[0:l,0],Res8_J3[0:l,1]*nmax)
errRes8_J3=np.sqrt(np.diag(pcRes8_J3))
modelRes8_J3=N(t3,*poRes8_J3)
poRes8_J4, pcRes8_J4 = curve_fit(N,Res8_J4[0:l,0],Res8_J4[0:l,1]*nmax)
errRes8_J4=np.sqrt(np.diag(pcRes8_J4))
poRes8_J5, pcRes8_J5 = curve_fit(N,Res8_J5[0:l,0],Res8_J5[0:l,1]*nmax)
errRes8_J5=np.sqrt(np.diag(pcRes8_J5))
modelRes8_J5=N(t5,*poRes8_J5)

#%% STALL 4 N0=6

n0=6

poStall4d_J0, pcStall4d_J0 = curve_fit(N,Stall4d_J0[0:l,0],Stall4d_J0[0:l,1]*nmax)
errStall4d_J0=np.sqrt(np.diag(pcStall4d_J0))
modelStall4d_J0=N(t0,*poStall4d_J0)
poStall4d_J1, pcStall4d_J1 = curve_fit(N,Stall4d_J1[0:l,0],Stall4d_J1[0:l,1]*nmax)
errStall4d_J1=np.sqrt(np.diag(pcStall4d_J1))
poStall4d_J2, pcStall4d_J2 = curve_fit(N,Stall4d_J2[0:l,0],Stall4d_J2[0:l,1]*nmax)
errStall4d_J2=np.sqrt(np.diag(pcStall4d_J2))
poStall4d_J3, pcStall4d_J3 = curve_fit(N,Stall4d_J3[0:l,0],Stall4d_J3[0:l,1]*nmax)
errStall4d_J3=np.sqrt(np.diag(pcStall4d_J3))
modelStall4d_J3=N(t3,*poStall4d_J3)
poStall4d_J4, pcStall4d_J4 = curve_fit(N,Stall4d_J4[0:l,0],Stall4d_J4[0:l,1]*nmax)
errStall4d_J4=np.sqrt(np.diag(pcStall4d_J4))
poStall4d_J5, pcStall4d_J5 = curve_fit(N,Stall4d_J5[0:l,0],Stall4d_J5[0:l,1]*nmax)
errStall4d_J5=np.sqrt(np.diag(pcStall4d_J5))
modelStall4d_J5=N(t5,*poStall4d_J5)

#%% STALL 8 N0=9

n0=9

poStall8d_J0, pcStall8d_J0 = curve_fit(N,Stall8d_J0[0:l,0],Stall8d_J0[0:l,1]*nmax)
errStall8d_J0=np.sqrt(np.diag(pcStall8d_J0))
modelStall8d_J0=N(t0,*poStall8d_J0)
poStall8d_J1, pcStall8d_J1 = curve_fit(N,Stall8d_J1[0:l,0],Stall8d_J1[0:l,1]*nmax)
errStall8d_J1=np.sqrt(np.diag(pcStall8d_J1))
poStall8d_J2, pcStall8d_J2 = curve_fit(N,Stall8d_J2[0:l,0],Stall8d_J2[0:l,1]*nmax)
errStall8d_J2=np.sqrt(np.diag(pcStall8d_J2))
poStall8d_J3, pcStall8d_J3 = curve_fit(N,Stall8d_J3[0:l,0],Stall8d_J3[0:l,1]*nmax)
errStall8d_J3=np.sqrt(np.diag(pcStall8d_J3))
modelStall8d_J3=N(t3,*poStall8d_J3)
poStall8d_J4, pcStall8d_J4 = curve_fit(N,Stall8d_J4[0:l,0],Stall8d_J4[0:l,1]*nmax)
errStall8d_J4=np.sqrt(np.diag(pcStall8d_J4))
poStall8d_J5, pcStall8d_J5 = curve_fit(N,Stall8d_J5[0:l,0],Stall8d_J5[0:l,1]*nmax)
errStall8d_J5=np.sqrt(np.diag(pcStall8d_J5))
modelStall8d_J5=N(t5,*poStall8d_J5)

plt.plot(Res8_J5[:,0],Res8_J5[:,1]*nmax)
plt.plot(Stall8d_J5[:,0],Stall8d_J5[:,1]*nmax)

#%%

tRes4[0]=poRes4_J0[0]
errtRes4[0]=errRes4_J0[0]
tRes4[1]=poRes4_J1[0]
errtRes4[1]=errRes4_J1[0]
tRes4[2]=poRes4_J2[0]
errtRes4[2]=errRes4_J2[0]
tRes4[3]=poRes4_J3[0]
errtRes4[3]=errRes4_J3[0]
tRes4[4]=poRes4_J4[0]
errtRes4[4]=errRes4_J4[0]
tRes4[5]=poRes4_J5[0]
errtRes4[5]=errRes4_J5[0]

tStall4d[0]=poStall4d_J0[0]
errtStall4d[0]=errStall4d_J0[0]
tStall4d[1]=poStall4d_J1[0]
errtStall4d[1]=errStall4d_J1[0]
tStall4d[2]=poStall4d_J2[0]
errtStall4d[2]=errStall4d_J2[0]
tStall4d[3]=poStall4d_J3[0]
errtStall4d[3]=errStall4d_J3[0]
tStall4d[4]=poStall4d_J4[0]
errtStall4d[4]=errStall4d_J4[0]
tStall4d[5]=poStall4d_J5[0]
errtStall4d[5]=errStall4d_J5[0]

tRes8[0]=poRes8_J0[0]
errtRes8[0]=errRes8_J0[0]
tRes8[1]=poRes8_J1[0]
errtRes8[1]=errRes8_J1[0]
tRes8[2]=poRes8_J2[0]
errtRes8[2]=errRes8_J2[0]
tRes8[3]=poRes8_J3[0]
errtRes8[3]=errRes8_J3[0]
tRes8[4]=poRes8_J4[0]
errtRes8[4]=errRes8_J4[0]
tRes8[5]=poRes8_J5[0]
errtRes8[5]=errRes8_J5[0]

tStall8d[0]=poStall8d_J0[0]
errtStall8d[0]=errStall8d_J0[0]
tStall8d[1]=poStall8d_J1[0]
errtStall8d[1]=errStall8d_J1[0]
tStall8d[2]=poStall8d_J2[0]
errtStall8d[2]=errStall8d_J2[0]
tStall8d[3]=poStall8d_J3[0]
errtStall8d[3]=errStall8d_J3[0]
tStall8d[4]=poStall8d_J4[0]
errtStall8d[4]=errStall8d_J4[0]
tStall8d[5]=poStall8d_J5[0]
errtStall8d[5]=errStall8d_J5[0]


Dt4d_new=(tStall4d-tRes4)/(tStall4d+tRes4)
errDt4d_new=2/(tStall4d+tRes4)**2*(tRes4*errtStall4d + tStall4d*errtRes4)

Dt8d_new=(tStall8d-tRes8)/(tStall8d+tRes8)
errDt8d_new=2/(tStall8d+tRes8)**2*(tRes8*errtStall8d + tStall8d*errtRes8)

J=np.arange(0,6,1)

plt.xlabel('J',fontsize=16)
plt.ylabel('$\Delta$',fontsize=16)
plt.title('$N_{max}=13$ ')
plt.grid()
plt.errorbar(J,Dt4d_new,errDt4d_new,ls='',marker='o',label='<N>=4, $<n_0>=6$',capsize=3)
plt.errorbar(J,Dt8d_new,errDt8d_new,ls='',marker='^',label='<N>=8, $<n_0>=9$',capsize=3)

plt.legend()

#%%

plt.xlabel('J',fontsize=16)
plt.ylabel('$\Delta$',fontsize=16)
plt.title('$N_{max}=13$ ')
plt.grid()
plt.plot(J,Dt4d_new,ls='',marker='o',color='blue',label='<N>=4, $<n_0>=6$')
plt.plot(J,Dt8d_new,ls='',marker='^',color='orange',label='<N>=8, $<n_0>=9$')
plt.plot(tRes4_old[:,0],Dt4d_old,ls='',marker='o',color='pink',label='<N>=4 (old), $<n_0>=6$')
plt.plot(tRes8_old[:,0],Dt8d_old,ls='',marker='^',color='green',label='<N>=8 (old), $<n_0>=9$')

plt.legend()


#%% DIFFERENT INITIAL CONDITIONS FOR STALL: N0=13 FIXED

FileS4f_13_J0=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=0.00_N=4_n0=13.dat']
DataS4f_13_J0=[]
FileS4f_13_J1=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=1.00_N=4_n0=13.dat']
DataS4f_13_J1=[]
FileS4f_13_J2=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=2.00_N=4_n0=13.dat']
DataS4f_13_J2=[]
FileS4f_13_J3=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=3.00_N=4_n0=13.dat']
DataS4f_13_J3=[]
FileS4f_13_J4=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=4.00_N=4_n0=13.dat']
DataS4f_13_J4=[]
FileS4f_13_J5=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=5.00_N=4_n0=13.dat']
DataS4f_13_J5=[]

for data_file in FileS4f_13_J0:
    DataS4f_13_J0.append(np.loadtxt(data_file))
for data_file in FileS4f_13_J1:
    DataS4f_13_J1.append(np.loadtxt(data_file))
for data_file in FileS4f_13_J2:
    DataS4f_13_J2.append(np.loadtxt(data_file))
for data_file in FileS4f_13_J3:
    DataS4f_13_J3.append(np.loadtxt(data_file))
for data_file in FileS4f_13_J4:
    DataS4f_13_J4.append(np.loadtxt(data_file))
for data_file in FileS4f_13_J5:
    DataS4f_13_J5.append(np.loadtxt(data_file))

Stall4f_13_J0=DataS4f_13_J0[0]
Stall4f_13_J1=DataS4f_13_J1[0]
Stall4f_13_J2=DataS4f_13_J2[0]
Stall4f_13_J3=DataS4f_13_J3[0]
Stall4f_13_J4=DataS4f_13_J4[0]
Stall4f_13_J5=DataS4f_13_J5[0]    


FileS8f_13_J0=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=0.00_N=8_n0=13.dat']
DataS8f_13_J0=[]
FileS8f_13_J1=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=1.00_N=8_n0=13.dat']
DataS8f_13_J1=[]
FileS8f_13_J2=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=2.00_N=8_n0=13.dat']
DataS8f_13_J2=[]
FileS8f_13_J3=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=3.00_N=8_n0=13.dat']
DataS8f_13_J3=[]
FileS8f_13_J4=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=4.00_N=8_n0=13.dat']
DataS8f_13_J4=[]
FileS8f_13_J5=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=5.00_N=8_n0=13.dat']
DataS8f_13_J5=[]

for data_file in FileS8f_13_J0:
    DataS8f_13_J0.append(np.loadtxt(data_file))
for data_file in FileS8f_13_J1:
    DataS8f_13_J1.append(np.loadtxt(data_file))
for data_file in FileS8f_13_J2:
    DataS8f_13_J2.append(np.loadtxt(data_file))
for data_file in FileS8f_13_J3:
    DataS8f_13_J3.append(np.loadtxt(data_file))
for data_file in FileS8f_13_J4:
    DataS8f_13_J4.append(np.loadtxt(data_file))
for data_file in FileS8f_13_J5:
    DataS8f_13_J5.append(np.loadtxt(data_file))
    
Stall8f_13_J0=DataS8f_13_J0[0]
Stall8f_13_J1=DataS8f_13_J1[0]
Stall8f_13_J2=DataS8f_13_J2[0]
Stall8f_13_J3=DataS8f_13_J3[0]
Stall8f_13_J4=DataS8f_13_J4[0]
Stall8f_13_J5=DataS8f_13_J5[0]


n0=13

poStall4f_13_J0, pcStall4f_13_J0 = curve_fit(N,Stall4f_13_J0[0:l,0],Stall4f_13_J0[0:l,1]*nmax)
errStall4f_13_J0=np.sqrt(np.diag(pcStall4f_13_J0))
poStall4f_13_J1, pcStall4f_13_J1 = curve_fit(N,Stall4f_13_J1[0:l,0],Stall4f_13_J1[0:l,1]*nmax)
errStall4f_13_J1=np.sqrt(np.diag(pcStall4f_13_J1))
poStall4f_13_J2, pcStall4f_13_J2 = curve_fit(N,Stall4f_13_J2[0:l,0],Stall4f_13_J2[0:l,1]*nmax)
errStall4f_13_J2=np.sqrt(np.diag(pcStall4f_13_J2))
poStall4f_13_J3, pcStall4f_13_J3 = curve_fit(N,Stall4f_13_J3[0:l,0],Stall4f_13_J3[0:l,1]*nmax)
errStall4f_13_J3=np.sqrt(np.diag(pcStall4f_13_J3))
poStall4f_13_J4, pcStall4f_13_J4 = curve_fit(N,Stall4f_13_J4[0:l,0],Stall4f_13_J4[0:l,1]*nmax)
errStall4f_13_J4=np.sqrt(np.diag(pcStall4f_13_J4))
poStall4f_13_J5, pcStall4f_13_J5 = curve_fit(N,Stall4f_13_J5[0:l,0],Stall4f_13_J5[0:l,1]*nmax)
errStall4f_13_J5=np.sqrt(np.diag(pcStall4f_13_J5))

poStall8f_13_J0, pcStall8f_13_J0 = curve_fit(N,Stall8f_13_J0[0:l,0],Stall8f_13_J0[0:l,1]*nmax)
errStall8f_13_J0=np.sqrt(np.diag(pcStall8f_13_J0))
poStall8f_13_J1, pcStall8f_13_J1 = curve_fit(N,Stall8f_13_J1[0:l,0],Stall8f_13_J1[0:l,1]*nmax)
errStall8f_13_J1=np.sqrt(np.diag(pcStall8f_13_J1))
poStall8f_13_J2, pcStall8f_13_J2 = curve_fit(N,Stall8f_13_J2[0:l,0],Stall8f_13_J2[0:l,1]*nmax)
errStall8f_13_J2=np.sqrt(np.diag(pcStall8f_13_J2))
poStall8f_13_J3, pcStall8f_13_J3 = curve_fit(N,Stall8f_13_J3[0:l,0],Stall8f_13_J3[0:l,1]*nmax)
errStall8f_13_J3=np.sqrt(np.diag(pcStall8f_13_J3))
poStall8f_13_J4, pcStall8f_13_J4 = curve_fit(N,Stall8f_13_J4[0:l,0],Stall8f_13_J4[0:l,1]*nmax)
errStall8f_13_J4=np.sqrt(np.diag(pcStall8f_13_J4))
poStall8f_13_J5, pcStall8f_13_J5 = curve_fit(N,Stall8f_13_J5[0:l,0],Stall8f_13_J5[0:l,1]*nmax)
errStall8f_13_J5=np.sqrt(np.diag(pcStall8f_13_J5))

tStall4f_13[0]=poStall4f_13_J0[0]
errtStall4f_13[0]=errStall4f_13_J0[0]
tStall4f_13[1]=poStall4f_13_J1[0]
errtStall4f_13[1]=errStall4f_13_J1[0]
tStall4f_13[2]=poStall4f_13_J2[0]
errtStall4f_13[2]=errStall4f_13_J2[0]
tStall4f_13[3]=poStall4f_13_J3[0]
errtStall4f_13[3]=errStall4f_13_J3[0]
tStall4f_13[4]=poStall4f_13_J4[0]
errtStall4f_13[4]=errStall4f_13_J4[0]
tStall4f_13[5]=poStall4f_13_J5[0]
errtStall4f_13[5]=errStall4f_13_J5[0]

tStall8f_13[0]=poStall8f_13_J0[0]
errtStall8f_13[0]=errStall8f_13_J0[0]
tStall8f_13[1]=poStall8f_13_J1[0]
errtStall8f_13[1]=errStall8f_13_J1[0]
tStall8f_13[2]=poStall8f_13_J2[0]
errtStall8f_13[2]=errStall8f_13_J2[0]
tStall8f_13[3]=poStall8f_13_J3[0]
errtStall8f_13[3]=errStall8f_13_J3[0]
tStall8f_13[4]=poStall8f_13_J4[0]
errtStall8f_13[4]=errStall8f_13_J4[0]
tStall8f_13[5]=poStall8f_13_J5[0]
errtStall8f_13[5]=errStall8f_13_J5[0]

Dt4f_13=(tStall4f_13-tRes4)/(tStall4f_13+tRes4)
errDt4f_13=2/(tStall4f_13+tRes4)**2*(tRes4*errtStall4f_13 + tStall4f_13*errtRes4)

Dt8f_13=(tStall8f_13-tRes8)/(tStall8f_13+tRes8)
errDt8f_13=2/(tStall8f_13+tRes8)**2*(tRes4*errtStall8f_13 + tStall8f_13*errtRes8)


plt.xlabel('J',fontsize=16)
plt.ylabel('$\Delta$',fontsize=16)
plt.title('$N_{max}=13$ ')
plt.grid()
plt.errorbar(J,Dt4f_13,errDt4f_13,ls='',marker='o',label='<N>=4, $n_0=13$',capsize=3)
plt.errorbar(J,Dt8f_13,errDt8f_13,ls='',marker='^',label='<N>=8, $n_0=13$',capsize=3)

plt.legend()

#%% HALF FILLING


FileR6_J0=['D:/Simulations/Without _depletion/Glauber/Traces_Res_J=0.00_N=6.dat']
DataR6_J0=[]
FileR6_J1=['D:/Simulations/Without _depletion/Glauber/Traces_Res_J=1.00_N=6.dat']
DataR6_J1=[]
FileR6_J2=['D:/Simulations/Without _depletion/Glauber/Traces_Res_J=2.00_N=6.dat']
DataR6_J2=[]
FileR6_J3=['D:/Simulations/Without _depletion/Glauber/Traces_Res_J=3.00_N=6.dat']
DataR6_J3=[]
FileR6_J4=['D:/Simulations/Without _depletion/Glauber/Traces_Res_J=4.00_N=6.dat']
DataR6_J4=[]
FileR6_J5=['D:/Simulations/Without _depletion/Glauber/Traces_Res_J=5.00_N=6.dat']
DataR6_J5=[]

for data_file in FileR6_J0:
    DataR6_J0.append(np.loadtxt(data_file))
for data_file in FileR6_J1:
    DataR6_J1.append(np.loadtxt(data_file))
for data_file in FileR6_J2:
    DataR6_J2.append(np.loadtxt(data_file))
for data_file in FileR6_J3:
    DataR6_J3.append(np.loadtxt(data_file))
for data_file in FileR6_J4:
    DataR6_J4.append(np.loadtxt(data_file))
for data_file in FileR6_J5:
    DataR6_J5.append(np.loadtxt(data_file))
    
Res6_J0=DataR6_J0[0]
Res6_J1=DataR6_J1[0]
Res6_J2=DataR6_J2[0]
Res6_J3=DataR6_J3[0]
Res6_J4=DataR6_J4[0]
Res6_J5=DataR6_J5[0]

FileS6f_13_J0=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=0.00_N=6_n0=13.dat']
DataS6f_13_J0=[]
FileS6f_13_J1=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=1.00_N=6_n0=13.dat']
DataS6f_13_J1=[]
FileS6f_13_J2=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=2.00_N=6_n0=13.dat']
DataS6f_13_J2=[]
FileS6f_13_J3=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=3.00_N=6_n0=13.dat']
DataS6f_13_J3=[]
FileS6f_13_J4=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=4.00_N=6_n0=13.dat']
DataS6f_13_J4=[]
FileS6f_13_J5=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=5.00_N=6_n0=13.dat']
DataS6f_13_J5=[]

for data_file in FileS6f_13_J0:
    DataS6f_13_J0.append(np.loadtxt(data_file))
for data_file in FileS6f_13_J1:
    DataS6f_13_J1.append(np.loadtxt(data_file))
for data_file in FileS6f_13_J2:
    DataS6f_13_J2.append(np.loadtxt(data_file))
for data_file in FileS6f_13_J3:
    DataS6f_13_J3.append(np.loadtxt(data_file))
for data_file in FileS6f_13_J4:
    DataS6f_13_J4.append(np.loadtxt(data_file))
for data_file in FileS6f_13_J5:
    DataS6f_13_J5.append(np.loadtxt(data_file))

Stall6f_13_J0=DataS6f_13_J0[0]
Stall6f_13_J1=DataS6f_13_J1[0]
Stall6f_13_J2=DataS6f_13_J2[0]
Stall6f_13_J3=DataS6f_13_J3[0]
Stall6f_13_J4=DataS6f_13_J4[0]
Stall6f_13_J5=DataS6f_13_J5[0] 

#%%
n0=0

poRes6_J0, pcRes6_J0 = curve_fit(N,Res6_J0[0:l,0],Res6_J0[0:l,1]*nmax)
errRes6_J0=np.sqrt(np.diag(pcRes6_J0))
modelRes6_J0=N(t0,*poRes6_J0)
poRes6_J1, pcRes6_J1 = curve_fit(N,Res6_J1[0:l,0],Res6_J1[0:l,1]*nmax)
errRes6_J1=np.sqrt(np.diag(pcRes6_J1))
poRes6_J2, pcRes6_J2 = curve_fit(N,Res6_J2[0:l,0],Res6_J2[0:l,1]*nmax)
errRes6_J2=np.sqrt(np.diag(pcRes6_J2))
poRes6_J3, pcRes6_J3 = curve_fit(N,Res6_J3[0:l,0],Res6_J3[0:l,1]*nmax)
errRes6_J3=np.sqrt(np.diag(pcRes6_J3))
modelRes6_J3=N(t3,*poRes6_J3)
poRes6_J4, pcRes6_J4 = curve_fit(N,Res6_J4[0:l,0],Res6_J4[0:l,1]*nmax)
errRes6_J4=np.sqrt(np.diag(pcRes6_J4))
poRes6_J5, pcRes6_J5 = curve_fit(N,Res6_J5[0:l,0],Res6_J5[0:l,1]*nmax)
errRes6_J5=np.sqrt(np.diag(pcRes6_J5))
modelRes6_J5=N(t5,*poRes6_J5)

#%%
n0 = 13

poStall6f_13_J0, pcStall6f_13_J0 = curve_fit(N,Stall6f_13_J0[0:l,0],Stall6f_13_J0[0:l,1]*nmax)
errStall6f_13_J0=np.sqrt(np.diag(pcStall6f_13_J0))
poStall6f_13_J1, pcStall6f_13_J1 = curve_fit(N,Stall6f_13_J1[0:l,0],Stall6f_13_J1[0:l,1]*nmax)
errStall6f_13_J1=np.sqrt(np.diag(pcStall6f_13_J1))
poStall6f_13_J2, pcStall6f_13_J2 = curve_fit(N,Stall6f_13_J2[0:l,0],Stall6f_13_J2[0:l,1]*nmax)
errStall6f_13_J2=np.sqrt(np.diag(pcStall6f_13_J2))
poStall6f_13_J3, pcStall6f_13_J3 = curve_fit(N,Stall6f_13_J3[0:l,0],Stall6f_13_J3[0:l,1]*nmax)
errStall6f_13_J3=np.sqrt(np.diag(pcStall6f_13_J3))
poStall6f_13_J4, pcStall6f_13_J4 = curve_fit(N,Stall6f_13_J4[0:l,0],Stall6f_13_J4[0:l,1]*nmax)
errStall6f_13_J4=np.sqrt(np.diag(pcStall6f_13_J4))
poStall6f_13_J5, pcStall6f_13_J5 = curve_fit(N,Stall6f_13_J5[0:l,0],Stall6f_13_J5[0:l,1]*nmax)
errStall6f_13_J5=np.sqrt(np.diag(pcStall6f_13_J5))


tRes6[0]=poRes6_J0[0]
errtRes6[0]=errRes6_J0[0]
tRes6[1]=poRes6_J1[0]
errtRes6[1]=errRes6_J1[0]
tRes6[2]=poRes6_J2[0]
errtRes6[2]=errRes6_J2[0]
tRes6[3]=poRes6_J3[0]
errtRes6[3]=errRes6_J3[0]
tRes6[4]=poRes6_J4[0]
errtRes6[4]=errRes6_J4[0]
tRes6[5]=poRes6_J5[0]
errtRes6[5]=errRes6_J5[0]

tStall6f_13[0]=poStall6f_13_J0[0]
errtStall6f_13[0]=errStall6f_13_J0[0]
tStall6f_13[1]=poStall6f_13_J1[0]
errtStall6f_13[1]=errStall6f_13_J1[0]
tStall6f_13[2]=poStall6f_13_J2[0]
errtStall6f_13[2]=errStall6f_13_J2[0]
tStall6f_13[3]=poStall6f_13_J3[0]
errtStall6f_13[3]=errStall6f_13_J3[0]
tStall6f_13[4]=poStall6f_13_J4[0]
errtStall6f_13[4]=errStall6f_13_J4[0]
tStall6f_13[5]=poStall6f_13_J5[0]
errtStall6f_13[5]=errStall6f_13_J5[0]


Dt6f_13=(tStall6f_13-tRes6)/(tStall6f_13+tRes6)
errDt6f_13=2/(tStall6f_13+tRes6)**2*(tRes4*errtStall6f_13 + tStall6f_13*errtRes6)

plt.xlabel('J',fontsize=16)
plt.ylabel('$\Delta$',fontsize=16)
plt.title('$N_{max}=13$ ')
plt.grid()
plt.errorbar(J,Dt4f_13,errDt4f_13,ls='',marker='o',label='<N>=4, $n_0=13$',capsize=3)
plt.errorbar(J,Dt8f_13,errDt8f_13,ls='',marker='^',label='<N>=8, $n_0=13$',capsize=3)
plt.errorbar(J,Dt6f_13,errDt6f_13,ls='',marker='s',label='<N>=6.5, $n_0=13$',capsize=3)

plt.legend()


#%% N=4 WITH N0=6, N=8 WITH N0=10 AND N=6 WITH N0=8

FileS4f_6_J0=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=0.00_N=4_n0=6.dat']
DataS4f_6_J0=[]
FileS4f_6_J1=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=1.00_N=4_n0=6.dat']
DataS4f_6_J1=[]
FileS4f_6_J2=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=2.00_N=4_n0=6.dat']
DataS4f_6_J2=[]
FileS4f_6_J3=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=3.00_N=4_n0=6.dat']
DataS4f_6_J3=[]
FileS4f_6_J4=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=4.00_N=4_n0=6.dat']
DataS4f_6_J4=[]
FileS4f_6_J5=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=5.00_N=4_n0=6.dat']
DataS4f_6_J5=[]

for data_file in FileS4f_6_J0:
    DataS4f_6_J0.append(np.loadtxt(data_file))
for data_file in FileS4f_6_J1:
    DataS4f_6_J1.append(np.loadtxt(data_file))
for data_file in FileS4f_6_J2:
    DataS4f_6_J2.append(np.loadtxt(data_file))
for data_file in FileS4f_6_J3:
    DataS4f_6_J3.append(np.loadtxt(data_file))
for data_file in FileS4f_6_J4:
    DataS4f_6_J4.append(np.loadtxt(data_file))
for data_file in FileS4f_6_J5:
    DataS4f_6_J5.append(np.loadtxt(data_file))

Stall4f_6_J0=DataS4f_6_J0[0]
Stall4f_6_J1=DataS4f_6_J1[0]
Stall4f_6_J2=DataS4f_6_J2[0]
Stall4f_6_J3=DataS4f_6_J3[0]
Stall4f_6_J4=DataS4f_6_J4[0]
Stall4f_6_J5=DataS4f_6_J5[0] 

FileS6f_8_J0=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=0.00_N=6_n0=8.dat']
DataS6f_8_J0=[]
FileS6f_8_J1=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=1.00_N=6_n0=8.dat']
DataS6f_8_J1=[]
FileS6f_8_J2=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=2.00_N=6_n0=8.dat']
DataS6f_8_J2=[]
FileS6f_8_J3=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=3.00_N=6_n0=8.dat']
DataS6f_8_J3=[]
FileS6f_8_J4=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=4.00_N=6_n0=8.dat']
DataS6f_8_J4=[]
FileS6f_8_J5=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=5.00_N=6_n0=8.dat']
DataS6f_8_J5=[]

for data_file in FileS6f_8_J0:
    DataS6f_8_J0.append(np.loadtxt(data_file))
for data_file in FileS6f_8_J1:
    DataS6f_8_J1.append(np.loadtxt(data_file))
for data_file in FileS6f_8_J2:
    DataS6f_8_J2.append(np.loadtxt(data_file))
for data_file in FileS6f_8_J3:
    DataS6f_8_J3.append(np.loadtxt(data_file))
for data_file in FileS6f_8_J4:
    DataS6f_8_J4.append(np.loadtxt(data_file))
for data_file in FileS6f_8_J5:
    DataS6f_8_J5.append(np.loadtxt(data_file))

Stall6f_8_J0=DataS6f_8_J0[0]
Stall6f_8_J1=DataS6f_8_J1[0]
Stall6f_8_J2=DataS6f_8_J2[0]
Stall6f_8_J3=DataS6f_8_J3[0]
Stall6f_8_J4=DataS6f_8_J4[0]
Stall6f_8_J5=DataS6f_8_J5[0]


FileS8f_10_J0=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=0.00_N=8_n0=10.dat']
DataS8f_10_J0=[]
FileS8f_10_J1=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=1.00_N=8_n0=10.dat']
DataS8f_10_J1=[]
FileS8f_10_J2=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=2.00_N=8_n0=10.dat']
DataS8f_10_J2=[]
FileS8f_10_J3=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=3.00_N=8_n0=10.dat']
DataS8f_10_J3=[]
FileS8f_10_J4=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=4.00_N=8_n0=10.dat']
DataS8f_10_J4=[]
FileS8f_10_J5=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=5.00_N=8_n0=10.dat']
DataS8f_10_J5=[]

for data_file in FileS8f_10_J0:
    DataS8f_10_J0.append(np.loadtxt(data_file))
for data_file in FileS8f_10_J1:
    DataS8f_10_J1.append(np.loadtxt(data_file))
for data_file in FileS8f_10_J2:
    DataS8f_10_J2.append(np.loadtxt(data_file))
for data_file in FileS8f_10_J3:
    DataS8f_10_J3.append(np.loadtxt(data_file))
for data_file in FileS8f_10_J4:
    DataS8f_10_J4.append(np.loadtxt(data_file))
for data_file in FileS8f_10_J5:
    DataS8f_10_J5.append(np.loadtxt(data_file))

Stall8f_10_J0=DataS8f_10_J0[0]
Stall8f_10_J1=DataS8f_10_J1[0]
Stall8f_10_J2=DataS8f_10_J2[0]
Stall8f_10_J3=DataS8f_10_J3[0]
Stall8f_10_J4=DataS8f_10_J4[0]
Stall8f_10_J5=DataS8f_10_J5[0]

#%%
n0=6

poStall4f_6_J0, pcStall4f_6_J0 = curve_fit(N,Stall4f_6_J0[0:l,0],Stall4f_6_J0[0:l,1]*nmax)
errStall4f_6_J0=np.sqrt(np.diag(pcStall4f_6_J0))
modelStall4f_6_J0=N(t0,*poStall4f_6_J0)
poStall4f_6_J1, pcStall4f_6_J1 = curve_fit(N,Stall4f_6_J1[0:l,0],Stall4f_6_J1[0:l,1]*nmax)
errStall4f_6_J1=np.sqrt(np.diag(pcStall4f_6_J1))
poStall4f_6_J2, pcStall4f_6_J2 = curve_fit(N,Stall4f_6_J2[0:l,0],Stall4f_6_J2[0:l,1]*nmax)
errStall4f_6_J2=np.sqrt(np.diag(pcStall4f_6_J2))
poStall4f_6_J3, pcStall4f_6_J3 = curve_fit(N,Stall4f_6_J3[0:l,0],Stall4f_6_J3[0:l,1]*nmax)
errStall4f_6_J3=np.sqrt(np.diag(pcStall4f_6_J3))
modelStall4f_6_J3=N(t3,*poStall4f_6_J3)
poStall4f_6_J4, pcStall4f_6_J4 = curve_fit(N,Stall4f_6_J4[0:l,0],Stall4f_6_J4[0:l,1]*nmax)
errStall4f_6_J4=np.sqrt(np.diag(pcStall4f_6_J4))
poStall4f_6_J5, pcStall4f_6_J5 = curve_fit(N,Stall4f_6_J5[0:l,0],Stall4f_6_J5[0:l,1]*nmax)
errStall4f_6_J5=np.sqrt(np.diag(pcStall4f_6_J5))
modelStall4f_6_J5=N(t5,*poStall4f_6_J5)

#%%
n0=9

poStall6f_8_J0, pcStall6f_8_J0 = curve_fit(N,Stall6f_8_J0[0:l,0],Stall6f_8_J0[0:l,1]*nmax)
errStall6f_8_J0=np.sqrt(np.diag(pcStall6f_8_J0))
modelStall6f_J0 = N(t0,*poStall6f_8_J0)
poStall6f_8_J1, pcStall6f_8_J1 = curve_fit(N,Stall6f_8_J1[0:l,0],Stall6f_8_J1[0:l,1]*nmax)
errStall6f_8_J1=np.sqrt(np.diag(pcStall6f_8_J1))
poStall6f_8_J2, pcStall6f_8_J2 = curve_fit(N,Stall6f_8_J2[0:l,0],Stall6f_8_J2[0:l,1]*nmax)
errStall6f_8_J2=np.sqrt(np.diag(pcStall6f_8_J2))
poStall6f_8_J3, pcStall6f_8_J3 = curve_fit(N,Stall6f_8_J3[0:l,0],Stall6f_8_J3[0:l,1]*nmax)
errStall6f_8_J3=np.sqrt(np.diag(pcStall6f_8_J3))
modelStall6f_J3 = N(t3,*poStall6f_8_J3)
poStall6f_8_J4, pcStall6f_8_J4 = curve_fit(N,Stall6f_8_J4[0:l,0],Stall6f_8_J4[0:l,1]*nmax)
errStall6f_8_J4=np.sqrt(np.diag(pcStall6f_8_J4))
poStall6f_8_J5, pcStall6f_8_J5 = curve_fit(N,Stall6f_8_J5[0:l,0],Stall6f_8_J5[0:l,1]*nmax)
errStall6f_8_J5=np.sqrt(np.diag(pcStall6f_8_J5))
modelStall6f_J5 = N(t5,*poStall6f_8_J5)

#%%
n0=10

poStall8f_10_J0, pcStall8f_10_J0 = curve_fit(N,Stall8f_10_J0[0:l,0],Stall8f_10_J0[0:l,1]*nmax)
poStall8f_10_J1, pcStall8f_10_J1 = curve_fit(N,Stall8f_10_J1[0:l,0],Stall8f_10_J1[0:l,1]*nmax)
poStall8f_10_J2, pcStall8f_10_J2 = curve_fit(N,Stall8f_10_J2[0:l,0],Stall8f_10_J2[0:l,1]*nmax)
poStall8f_10_J3, pcStall8f_10_J3 = curve_fit(N,Stall8f_10_J3[0:l,0],Stall8f_10_J3[0:l,1]*nmax)
poStall8f_10_J4, pcStall8f_10_J4 = curve_fit(N,Stall8f_10_J4[0:l,0],Stall8f_10_J4[0:l,1]*nmax)
poStall8f_10_J5, pcStall8f_10_J5 = curve_fit(N,Stall8f_10_J5[0:l,0],Stall8f_10_J5[0:l,1]*nmax)


#%%

#

tStall4f_6[0]=poStall4f_6_J0[0]
errtStall4f_6[0]=errStall4f_6_J0[0]
tStall4f_6[1]=poStall4f_6_J1[0]
errtStall4f_6[1]=errStall4f_6_J1[0]
tStall4f_6[2]=poStall4f_6_J2[0]
errtStall4f_6[2]=errStall4f_6_J2[0]
tStall4f_6[3]=poStall4f_6_J3[0]
errtStall4f_6[3]=errStall4f_6_J3[0]
tStall4f_6[4]=poStall4f_6_J4[0]
errtStall4f_6[4]=errStall4f_6_J4[0]
tStall4f_6[5]=poStall4f_6_J5[0]
errtStall4f_6[5]=errStall4f_6_J5[0]

tStall6f_8[0]=poStall6f_8_J0[0]
errtStall6f_8[0]=errStall6f_8_J0[0]
tStall6f_8[1]=poStall6f_8_J1[0]
errtStall6f_8[1]=errStall6f_8_J1[0]
tStall6f_8[2]=poStall6f_8_J2[0]
errtStall6f_8[2]=errStall6f_8_J2[0]
tStall6f_8[3]=poStall6f_8_J3[0]
errtStall6f_8[3]=errStall6f_8_J3[0]
tStall6f_8[4]=poStall6f_8_J4[0]
errtStall6f_8[4]=errStall6f_8_J4[0]
tStall6f_8[5]=poStall6f_8_J5[0]
errtStall6f_8[5]=errStall6f_8_J5[0]

tStall8f_10[0]=poStall8f_10_J0[0]
tStall8f_10[1]=poStall8f_10_J1[0]
tStall8f_10[2]=poStall8f_10_J2[0]
tStall8f_10[3]=poStall8f_10_J3[0]
tStall8f_10[4]=poStall8f_10_J4[0]
tStall8f_10[5]=poStall8f_10_J5[0]


Dt4f_6=(tStall4f_6-tRes4)/(tStall4f_6+tRes4)
errDt4f_6=2/(tStall4f_6+tRes4)**2*(tRes4*errtStall4f_6 + tStall4f_6*errtRes4)

Dt6f_8=(tStall6f_8-tRes6)/(tStall6f_8+tRes6)
errDt6f_8=2/(tStall6f_8+tRes6)**2*(tRes6*errtStall6f_8 + tStall6f_8*errtRes6)

Dt8f_10=(tStall8f_10-tRes8)/(tStall8f_10+tRes8) 

plt.xlabel('J',fontsize=16)
plt.ylabel('$\Delta$',fontsize=16)
plt.title('$N_{max}=13$ ')
plt.grid()
plt.errorbar(J,Dt4f_6,errDt4f_6,ls='',marker='o',label='<N>=4, $n_0=6$',capsize=3)
plt.plot(J,Dt8f_10,ls='',marker='^',label='<N>=8, $n_0=10$')
plt.plot(J,Dt6f_8,ls='',marker='s',label='<N>=6.5, $n_0=9$')

plt.legend()

#%%

FileS8f_9_J0=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=0.00_N=8_n0=9.dat']
DataS8f_9_J0=[]
FileS8f_9_J1=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=1.00_N=8_n0=9.dat']
DataS8f_9_J1=[]
FileS8f_9_J2=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=2.00_N=8_n0=9.dat']
DataS8f_9_J2=[]
FileS8f_9_J3=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=3.00_N=8_n0=9.dat']
DataS8f_9_J3=[]
FileS8f_9_J4=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=4.00_N=8_n0=9.dat']
DataS8f_9_J4=[]
FileS8f_9_J5=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=5.00_N=8_n0=9.dat']
DataS8f_9_J5=[]

for data_file in FileS8f_9_J0:
    DataS8f_9_J0.append(np.loadtxt(data_file))
for data_file in FileS8f_9_J1:
    DataS8f_9_J1.append(np.loadtxt(data_file))
for data_file in FileS8f_9_J2:
    DataS8f_9_J2.append(np.loadtxt(data_file))
for data_file in FileS8f_9_J3:
    DataS8f_9_J3.append(np.loadtxt(data_file))
for data_file in FileS8f_9_J4:
    DataS8f_9_J4.append(np.loadtxt(data_file))
for data_file in FileS8f_9_J5:
    DataS8f_9_J5.append(np.loadtxt(data_file))

Stall8f_9_J0=DataS8f_9_J0[0]
Stall8f_9_J1=DataS8f_9_J1[0]
Stall8f_9_J2=DataS8f_9_J2[0]
Stall8f_9_J3=DataS8f_9_J3[0]
Stall8f_9_J4=DataS8f_9_J4[0]
Stall8f_9_J5=DataS8f_9_J5[0]

n0=9

poStall8f_9_J0, pcStall8f_9_J0 = curve_fit(N,Stall8f_9_J0[0:l,0],Stall8f_9_J0[0:l,1]*nmax)
errStall8f_9_J0=np.sqrt(np.diag(pcStall8f_9_J0))
modelStall8f_9_J0=N(t0,*poStall8f_9_J0)
poStall8f_9_J1, pcStall8f_9_J1 = curve_fit(N,Stall8f_9_J1[0:l,0],Stall8f_9_J1[0:l,1]*nmax)
errStall8f_9_J1=np.sqrt(np.diag(pcStall8f_9_J1))
poStall8f_9_J2, pcStall8f_9_J2 = curve_fit(N,Stall8f_9_J2[0:l,0],Stall8f_9_J2[0:l,1]*nmax)
errStall8f_9_J2=np.sqrt(np.diag(pcStall8f_9_J2))
poStall8f_9_J3, pcStall8f_9_J3 = curve_fit(N,Stall8f_9_J3[0:l,0],Stall8f_9_J3[0:l,1]*nmax)
errStall8f_9_J3=np.sqrt(np.diag(pcStall8f_9_J3))
modelStall8f_9_J3=N(t3,*poStall8f_9_J3)
poStall8f_9_J4, pcStall8f_9_J4 = curve_fit(N,Stall8f_9_J4[0:l,0],Stall8f_9_J4[0:l,1]*nmax)
errStall8f_9_J4=np.sqrt(np.diag(pcStall8f_9_J4))
poStall8f_9_J5, pcStall8f_9_J5 = curve_fit(N,Stall8f_9_J5[0:l,0],Stall8f_9_J5[0:l,1]*nmax)
errStall8f_9_J5=np.sqrt(np.diag(pcStall8f_9_J5))
modelStall8f_9_J5=N(t5,*poStall8f_9_J5)

tStall8f_9[0]=poStall8f_9_J0[0]
errtStall8f_9[0]=errStall8f_9_J0[0]
tStall8f_9[1]=poStall8f_9_J1[0]
errtStall8f_9[1]=errStall8f_9_J1[0]
tStall8f_9[2]=poStall8f_9_J2[0]
errtStall8f_9[2]=errStall8f_9_J2[0]
tStall8f_9[3]=poStall8f_9_J3[0]
errtStall8f_9[3]=errStall8f_9_J3[0]
tStall8f_9[4]=poStall8f_9_J4[0]
errtStall8f_9[4]=errStall8f_9_J4[0]
tStall8f_9[5]=poStall8f_10_J5[0]
errtStall8f_9[5]=errStall8f_9_J5[0]

Dt8f_9=(tStall8f_9-tRes8)/(tStall8f_9+tRes8) 
errDt8f_9=2/(tStall8f_9+tRes8)**2*(tRes8*errtStall8f_9 + tStall8f_9*errtRes8)

plt.plot(Res4_J5[0:l,0],Res4_J5[0:l,1]*nmax)
plt.plot(Stall4d_J5[0:l,0],Stall4d_J5[0:l,1]*nmax)
plt.plot(Stall4f_6_J5[0:l,0],Stall4f_6_J5[0:l,1]*nmax)

#%%

FileS6d_9_J0=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=0.00_N=6.dat']
DataS6d_9_J0=[]
FileS6d_9_J1=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=1.00_N=6.dat']
DataS6d_9_J1=[]
FileS6d_9_J2=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=2.00_N=6.dat']
DataS6d_9_J2=[]
FileS6d_9_J3=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=3.00_N=6.dat']
DataS6d_9_J3=[]
FileS6d_9_J4=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=4.00_N=6.dat']
DataS6d_9_J4=[]
FileS6d_9_J5=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=5.00_N=6.dat']
DataS6d_9_J5=[]

for data_file in FileS6d_9_J0:
    DataS6d_9_J0.append(np.loadtxt(data_file))
for data_file in FileS6d_9_J1:
    DataS6d_9_J1.append(np.loadtxt(data_file))
for data_file in FileS6d_9_J2:
    DataS6d_9_J2.append(np.loadtxt(data_file))
for data_file in FileS6d_9_J3:
    DataS6d_9_J3.append(np.loadtxt(data_file))
for data_file in FileS6d_9_J4:
    DataS6d_9_J4.append(np.loadtxt(data_file))
for data_file in FileS6d_9_J5:
    DataS6d_9_J5.append(np.loadtxt(data_file))

Stall6d_9_J0=DataS6d_9_J0[0]
Stall6d_9_J1=DataS6d_9_J1[0]
Stall6d_9_J2=DataS6d_9_J2[0]
Stall6d_9_J3=DataS6d_9_J3[0]
Stall6d_9_J4=DataS6d_9_J4[0]
Stall6d_9_J5=DataS6d_9_J5[0]

n0=9

poStall6d_9_J0, pcStall6d_9_J0 = curve_fit(N,Stall6d_9_J0[0:l,0],Stall6d_9_J0[0:l,1]*nmax)
errStall6d_9_J0=np.sqrt(np.diag(pcStall6d_9_J0))
modelStall6d_9_J0=N(t0,*poStall6d_9_J0)
poStall6d_9_J1, pcStall6d_9_J1 = curve_fit(N,Stall6d_9_J1[0:l,0],Stall6d_9_J1[0:l,1]*nmax)
errStall6d_9_J1=np.sqrt(np.diag(pcStall6d_9_J1))
poStall6d_9_J2, pcStall6d_9_J2 = curve_fit(N,Stall6d_9_J2[0:l,0],Stall6d_9_J2[0:l,1]*nmax)
errStall6d_9_J2=np.sqrt(np.diag(pcStall6d_9_J2))
poStall6d_9_J3, pcStall6d_9_J3 = curve_fit(N,Stall6d_9_J3[0:l,0],Stall6d_9_J3[0:l,1]*nmax)
errStall6d_9_J3=np.sqrt(np.diag(pcStall6d_9_J3))
modelStall6d_9_J3=N(t3,*poStall6d_9_J3)
poStall6d_9_J4, pcStall6d_9_J4 = curve_fit(N,Stall6d_9_J4[0:l,0],Stall6d_9_J4[0:l,1]*nmax)
errStall6d_9_J4=np.sqrt(np.diag(pcStall6d_9_J4))
poStall6d_9_J5, pcStall6d_9_J5 = curve_fit(N,Stall6d_9_J5[0:l,0],Stall6d_9_J5[0:l,1]*nmax)
errStall6d_9_J5=np.sqrt(np.diag(pcStall6d_9_J5))
modelStall6d_9_J5=N(t5,*poStall6d_9_J5)

tStall6d_9[0]=poStall6d_9_J0[0]
errtStall6d_9[0]=errStall6d_9_J0[0]
tStall6d_9[1]=poStall6d_9_J1[0]
errtStall6d_9[1]=errStall6d_9_J1[0]
tStall6d_9[2]=poStall6d_9_J2[0]
errtStall6d_9[2]=errStall6d_9_J2[0]
tStall6d_9[3]=poStall6d_9_J3[0]
errtStall6d_9[3]=errStall6d_9_J3[0]
tStall6d_9[4]=poStall6d_9_J4[0]
errtStall6d_9[4]=errStall6d_9_J4[0]
tStall6d_9[5]=poStall6d_9_J5[0]
errtStall6d_9[5]=errStall6d_9_J5[0]

Dt6d_9=(tStall6d_9-tRes6)/(tStall6d_9+tRes6)
errDt6d_9=2/(tStall6d_9+tRes6)**2*(tRes6*errtStall6d_9 + tStall6d_9*errtRes6)


#%%

def t_th(x):
    return np.exp(x)/2

cJ = np.linspace(0,5,100000)

plt.xlabel('J',fontsize=16)
plt.ylabel(r'$\tau$',fontsize=16)
plt.title('$<N> = 8,N_{max}=13$ ')
plt.grid()
plt.errorbar(J,np.log(tRes8),errtRes8/tRes8,marker='o',ls='',label='Resurrection',capsize=3)
plt.errorbar(J,np.log(tStall8d),errtStall8d/tStall8d,marker='^',ls='',label='Stall, $<n_0>=9$',capsize=3)
plt.errorbar(J,np.log(tStall8f_9),errtStall8f_9/tStall8f_9,marker='s',ls='',label='Stall, $n_0=9$',capsize=3)
plt.plot(cJ,np.log(15*t_th(cJ)),label=r'$log(\tau)\propto$ J')

plt.legend(prop={'size': 12})

plt.show()

plt.xlabel('J',fontsize=16)
plt.ylabel('$\Delta$',fontsize=16)
plt.title('$<N> = 8,N_{max}=13$ ')
plt.grid()
plt.errorbar(J,Dt8d_new,errDt8d_new,marker='^',ls='',label='$<n_0>=9$',color='darkorange',capsize=3)
plt.errorbar(J,Dt8f_9,errDt8f_9,marker='s',ls='',label='$n_0=9$',color='green',capsize=3)

plt.legend(prop={'size': 12})

plt.show()


#%%

plt.xlabel('J',fontsize=16)
plt.ylabel(r'$\tau$',fontsize=16)
plt.title('$<N> = 4,N_{max}=13$ ')
plt.grid()
plt.errorbar(J,np.log(tRes4),errtRes4/tRes4,marker='o',ls='',label='Resurrection',capsize=3)
plt.errorbar(J,np.log(tStall4d),errtStall4d/tStall4d,marker='^',ls='',label='Stall, $<n_0>=6$',capsize=3)
plt.errorbar(J,np.log(tStall4f_6),errtStall4f_6/tStall4f_6,marker='s',ls='',label='Stall, $n_0=6$',capsize=3)
plt.plot(cJ,np.log(15*t_th(cJ)),label=r'$log(\tau)\propto$ J')

plt.legend(prop={'size': 12})

plt.show()

plt.xlabel('J',fontsize=16)
plt.ylabel('$\Delta$',fontsize=16)
plt.title('$<N> = 4,N_{max}=13$ ')
plt.grid()
plt.errorbar(J,Dt4d_new,errDt4d_new,marker='^',ls='',label='$<n_0>=6$',color='darkorange',capsize=3)
plt.errorbar(J,Dt4f_6,errDt4f_6,marker='s',ls='',label='$n_0=6$',color='green',capsize=3)

plt.legend(prop={'size': 12})

plt.show()
#%%

plt.xlabel('J',fontsize=16)
plt.ylabel(r'log($\tau$)',fontsize=16)
plt.title('$<N> = 6.5,N_{max}=13$ ')
plt.grid()
plt.errorbar(J,np.log(tRes6),errtRes6/tRes6,marker='o',ls='',label='Resurrection',capsize=3)
plt.errorbar(J,np.log(tStall6d_9),errtStall6d_9/tStall6d_9,marker='^',ls='',label='Release, $<n_0>=9$',capsize=3)
plt.errorbar(J,np.log(tStall6f_8),errtStall6f_8/tStall6f_8,marker='s',ls='',label='Release, $n_0=9$',capsize=3)
plt.plot(cJ,np.log(15*t_th(cJ)),label=r'$log(\tau)\propto$ J')

plt.legend(prop={'size': 12})

plt.show()

plt.xlabel('J',fontsize=16)
plt.ylabel('$\Delta$',fontsize=16)
plt.title('$<N> = 6.5,N_{max}=13$ ')
plt.grid()
plt.errorbar(J,Dt6d_9,errDt6d_9,marker='^',ls='',label='$<n_0>=9$',color='darkorange',capsize=3)
plt.errorbar(J,Dt6f_8,errDt6f_8,marker='s',ls='',label='$n_0=9$',color='green',capsize=3)

plt.legend(prop={'size': 12})

plt.show()

#%% N0=13 distribution

FileS4d_13_J0=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=0.00_N=4_n0=13_d.dat']
DataS4d_13_J0=[]
FileS4d_13_J1=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=1.00_N=4_n0=13_d.dat']
DataS4d_13_J1=[]
FileS4d_13_J2=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=2.00_N=4_n0=13_d.dat']
DataS4d_13_J2=[]
FileS4d_13_J3=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=3.00_N=4_n0=13_d.dat']
DataS4d_13_J3=[]
FileS4d_13_J4=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=4.00_N=4_n0=13_d.dat']
DataS4d_13_J4=[]
FileS4d_13_J5=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=5.00_N=4_n0=13_d.dat']
DataS4d_13_J5=[]

for data_file in FileS4d_13_J0:
    DataS4d_13_J0.append(np.loadtxt(data_file))
for data_file in FileS4d_13_J1:
    DataS4d_13_J1.append(np.loadtxt(data_file))
for data_file in FileS4d_13_J2:
    DataS4d_13_J2.append(np.loadtxt(data_file))
for data_file in FileS4d_13_J3:
    DataS4d_13_J3.append(np.loadtxt(data_file))
for data_file in FileS4d_13_J4:
    DataS4d_13_J4.append(np.loadtxt(data_file))
for data_file in FileS4d_13_J5:
    DataS4d_13_J5.append(np.loadtxt(data_file))

Stall4d_13_J0=DataS4d_13_J0[0]
Stall4d_13_J1=DataS4d_13_J1[0]
Stall4d_13_J2=DataS4d_13_J2[0]
Stall4d_13_J3=DataS4d_13_J3[0]
Stall4d_13_J4=DataS4d_13_J4[0]
Stall4d_13_J5=DataS4d_13_J5[0]    


FileS8d_13_J0=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=0.00_N=8_n0=13_d.dat']
DataS8d_13_J0=[]
FileS8d_13_J1=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=1.00_N=8_n0=13_d.dat']
DataS8d_13_J1=[]
FileS8d_13_J2=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=2.00_N=8_n0=13_d.dat']
DataS8d_13_J2=[]
FileS8d_13_J3=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=3.00_N=8_n0=13_d.dat']
DataS8d_13_J3=[]
FileS8d_13_J4=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=4.00_N=8_n0=13_d.dat']
DataS8d_13_J4=[]
FileS8d_13_J5=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=5.00_N=8_n0=13_d.dat']
DataS8d_13_J5=[]

for data_file in FileS8d_13_J0:
    DataS8d_13_J0.append(np.loadtxt(data_file))
for data_file in FileS8d_13_J1:
    DataS8d_13_J1.append(np.loadtxt(data_file))
for data_file in FileS8d_13_J2:
    DataS8d_13_J2.append(np.loadtxt(data_file))
for data_file in FileS8d_13_J3:
    DataS8d_13_J3.append(np.loadtxt(data_file))
for data_file in FileS8d_13_J4:
    DataS8d_13_J4.append(np.loadtxt(data_file))
for data_file in FileS8d_13_J5:
    DataS8d_13_J5.append(np.loadtxt(data_file))
    
Stall8d_13_J0=DataS8d_13_J0[0]
Stall8d_13_J1=DataS8d_13_J1[0]
Stall8d_13_J2=DataS8d_13_J2[0]
Stall8d_13_J3=DataS8d_13_J3[0]
Stall8d_13_J4=DataS8d_13_J4[0]
Stall8d_13_J5=DataS8d_13_J5[0]

FileS6d_13_J0=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=0.00_N=6_n0=13_d.dat']
DataS6d_13_J0=[]
FileS6d_13_J1=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=1.00_N=6_n0=13_d.dat']
DataS6d_13_J1=[]
FileS6d_13_J2=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=2.00_N=6_n0=13_d.dat']
DataS6d_13_J2=[]
FileS6d_13_J3=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=3.00_N=6_n0=13_d.dat']
DataS6d_13_J3=[]
FileS6d_13_J4=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=4.00_N=6_n0=13_d.dat']
DataS6d_13_J4=[]
FileS6d_13_J5=['D:/Simulations/Without _depletion/Glauber/Traces_Stall_J=5.00_N=6_n0=13_d.dat']
DataS6d_13_J5=[]

for data_file in FileS6d_13_J0:
    DataS6d_13_J0.append(np.loadtxt(data_file))
for data_file in FileS6d_13_J1:
    DataS6d_13_J1.append(np.loadtxt(data_file))
for data_file in FileS6d_13_J2:
    DataS6d_13_J2.append(np.loadtxt(data_file))
for data_file in FileS6d_13_J3:
    DataS6d_13_J3.append(np.loadtxt(data_file))
for data_file in FileS6d_13_J4:
    DataS6d_13_J4.append(np.loadtxt(data_file))
for data_file in FileS6d_13_J5:
    DataS6d_13_J5.append(np.loadtxt(data_file))

Stall6d_13_J0=DataS6d_13_J0[0]
Stall6d_13_J1=DataS6d_13_J1[0]
Stall6d_13_J2=DataS6d_13_J2[0]
Stall6d_13_J3=DataS6d_13_J3[0]
Stall6d_13_J4=DataS6d_13_J4[0]
Stall6d_13_J5=DataS6d_13_J5[0] 


n0=13

poStall4d_13_J0, pcStall4d_13_J0 = curve_fit(N,Stall4d_13_J0[0:l,0],Stall4d_13_J0[0:l,1]*nmax)
errStall4d_13_J0=np.sqrt(np.diag(pcStall4d_13_J0))
poStall4d_13_J1, pcStall4d_13_J1 = curve_fit(N,Stall4d_13_J1[0:l,0],Stall4d_13_J1[0:l,1]*nmax)
errStall4d_13_J1=np.sqrt(np.diag(pcStall4d_13_J1))
poStall4d_13_J2, pcStall4d_13_J2 = curve_fit(N,Stall4d_13_J2[0:l,0],Stall4d_13_J2[0:l,1]*nmax)
errStall4d_13_J2=np.sqrt(np.diag(pcStall4d_13_J2))
poStall4d_13_J3, pcStall4d_13_J3 = curve_fit(N,Stall4d_13_J3[0:l,0],Stall4d_13_J3[0:l,1]*nmax)
errStall4d_13_J3=np.sqrt(np.diag(pcStall4d_13_J3))
poStall4d_13_J4, pcStall4d_13_J4 = curve_fit(N,Stall4d_13_J4[0:l,0],Stall4d_13_J4[0:l,1]*nmax)
errStall4d_13_J4=np.sqrt(np.diag(pcStall4d_13_J4))
poStall4d_13_J5, pcStall4d_13_J5 = curve_fit(N,Stall4d_13_J5[0:l,0],Stall4d_13_J5[0:l,1]*nmax)
errStall4d_13_J5=np.sqrt(np.diag(pcStall4d_13_J5))

poStall8d_13_J0, pcStall8d_13_J0 = curve_fit(N,Stall8d_13_J0[0:l,0],Stall8d_13_J0[0:l,1]*nmax)
errStall8d_13_J0=np.sqrt(np.diag(pcStall8d_13_J0))
poStall8d_13_J1, pcStall8d_13_J1 = curve_fit(N,Stall8d_13_J1[0:l,0],Stall8d_13_J1[0:l,1]*nmax)
errStall8d_13_J1=np.sqrt(np.diag(pcStall8d_13_J1))
poStall8d_13_J2, pcStall8d_13_J2 = curve_fit(N,Stall8d_13_J2[0:l,0],Stall8d_13_J2[0:l,1]*nmax)
errStall8d_13_J2=np.sqrt(np.diag(pcStall8d_13_J2))
poStall8d_13_J3, pcStall8d_13_J3 = curve_fit(N,Stall8d_13_J3[0:l,0],Stall8d_13_J3[0:l,1]*nmax)
errStall8d_13_J3=np.sqrt(np.diag(pcStall8d_13_J3))
poStall8d_13_J4, pcStall8d_13_J4 = curve_fit(N,Stall8d_13_J4[0:l,0],Stall8d_13_J4[0:l,1]*nmax)
errStall8d_13_J4=np.sqrt(np.diag(pcStall8d_13_J4))
poStall8d_13_J5, pcStall8d_13_J5 = curve_fit(N,Stall8d_13_J5[0:l,0],Stall8d_13_J5[0:l,1]*nmax)
errStall8d_13_J5=np.sqrt(np.diag(pcStall8d_13_J5))

poStall6d_13_J0, pcStall6d_13_J0 = curve_fit(N,Stall6d_13_J0[0:l,0],Stall6d_13_J0[0:l,1]*nmax)
errStall6d_13_J0=np.sqrt(np.diag(pcStall6d_13_J0))
poStall6d_13_J1, pcStall6d_13_J1 = curve_fit(N,Stall6d_13_J1[0:l,0],Stall6d_13_J1[0:l,1]*nmax)
errStall6d_13_J1=np.sqrt(np.diag(pcStall6d_13_J1))
poStall6d_13_J2, pcStall6d_13_J2 = curve_fit(N,Stall6d_13_J2[0:l,0],Stall6d_13_J2[0:l,1]*nmax)
errStall6d_13_J2=np.sqrt(np.diag(pcStall6d_13_J2))
poStall6d_13_J3, pcStall6d_13_J3 = curve_fit(N,Stall6d_13_J3[0:l,0],Stall6d_13_J3[0:l,1]*nmax)
errStall6d_13_J3=np.sqrt(np.diag(pcStall6d_13_J3))
poStall6d_13_J4, pcStall6d_13_J4 = curve_fit(N,Stall6d_13_J4[0:l,0],Stall6d_13_J4[0:l,1]*nmax)
errStall6d_13_J4=np.sqrt(np.diag(pcStall6d_13_J4))
poStall6d_13_J5, pcStall6d_13_J5 = curve_fit(N,Stall6d_13_J5[0:l,0],Stall6d_13_J5[0:l,1]*nmax)
errStall6d_13_J5=np.sqrt(np.diag(pcStall6d_13_J5))

tStall4d_13[0]=poStall4d_13_J0[0]
errtStall4d_13[0]=errStall4d_13_J0[0]
tStall4d_13[1]=poStall4d_13_J1[0]
errtStall4d_13[1]=errStall4d_13_J1[0]
tStall4d_13[2]=poStall4d_13_J2[0]
errtStall4d_13[2]=errStall4d_13_J2[0]
tStall4d_13[3]=poStall4d_13_J3[0]
errtStall4d_13[3]=errStall4d_13_J3[0]
tStall4d_13[4]=poStall4d_13_J4[0]
errtStall4d_13[4]=errStall4d_13_J4[0]
tStall4d_13[5]=poStall4d_13_J5[0]
errtStall4d_13[5]=errStall4d_13_J5[0]

tStall8d_13[0]=poStall8d_13_J0[0]
errtStall8d_13[0]=errStall8d_13_J0[0]
tStall8d_13[1]=poStall8d_13_J1[0]
errtStall8d_13[1]=errStall8d_13_J1[0]
tStall8d_13[2]=poStall8d_13_J2[0]
errtStall8d_13[2]=errStall8d_13_J2[0]
tStall8d_13[3]=poStall8d_13_J3[0]
errtStall8d_13[3]=errStall8d_13_J3[0]
tStall8d_13[4]=poStall8d_13_J4[0]
errtStall8d_13[4]=errStall8d_13_J4[0]
tStall8d_13[5]=poStall8d_13_J5[0]
errtStall8d_13[5]=errStall8d_13_J5[0]

tStall6d_13[0]=poStall6d_13_J0[0]
errtStall6d_13[0]=errStall6d_13_J0[0]
tStall6d_13[1]=poStall6d_13_J1[0]
errtStall6d_13[1]=errStall6d_13_J1[0]
tStall6d_13[2]=poStall6d_13_J2[0]
errtStall6d_13[2]=errStall6d_13_J2[0]
tStall6d_13[3]=poStall6d_13_J3[0]
errtStall6d_13[3]=errStall6d_13_J3[0]
tStall6d_13[4]=poStall6d_13_J4[0]
errtStall6d_13[4]=errStall6d_13_J4[0]
tStall6d_13[5]=poStall6d_13_J5[0]
errtStall6d_13[5]=errStall6d_13_J5[0]


Dt4d_13=(tStall4d_13-tRes4)/(tStall4d_13+tRes4)
errDt4d_13=2/(tStall4d_13+tRes4)**2*(tRes4*errtStall4d_13 + tStall4d_13*errtRes4)

Dt8d_13=(tStall8d_13-tRes8)/(tStall8d_13+tRes8)
errDt8d_13=2/(tStall8d_13+tRes8)**2*(tRes8*errtStall8d_13 + tStall8d_13*errtRes8)

Dt6d_13=(tStall6d_13-tRes6)/(tStall6d_13+tRes6)
errDt6d_13=2/(tStall6d_13+tRes6)**2*(tRes6*errtStall6d_13 + tStall6d_13*errtRes6)

#%%

plt.xlabel('J',fontsize=16)
plt.ylabel(r'$\tau$',fontsize=16)
plt.title('$<N> = 4,N_{max}=13$ ')
plt.grid()
plt.errorbar(J,np.log(tRes4),errtRes4/tRes4,marker='o',ls='',label='Resurrection',capsize=3)
plt.errorbar(J,np.log(tStall4d_13),errtStall4d_13/tStall4d_13,marker='^',ls='',label='Stall, $<n_0>=13$',capsize=3)
plt.errorbar(J,np.log(tStall4f_13),errtStall4f_13/tStall4f_13,marker='s',ls='',label='Stall, $n_0=13$',alpha=0.5,capsize=3)
plt.plot(cJ,np.log(15*t_th(cJ)),label=r'$log(\tau)\propto$ J')

plt.legend(prop={'size': 12})

plt.show()

plt.xlabel('J',fontsize=16)
plt.ylabel('$\Delta$',fontsize=16)
plt.title('$<N> = 4,N_{max}=13$ ')
plt.grid()
plt.errorbar(J,Dt4d_13,errDt4d_13,marker='^',ls='',label='$<n_0>=13$',color='darkorange',capsize=3)
plt.errorbar(J,Dt4f_13,errDt4f_13,marker='s',ls='',label='$n_0=13$',color='green',alpha=0.5,capsize=3)

plt.legend(prop={'size': 12})

plt.show()

#%%

plt.xlabel('J',fontsize=16)
plt.ylabel(r'$\tau$',fontsize=16)
plt.title('$<N> = 6.5,N_{max}=13$ ')
plt.grid()
plt.errorbar(J,np.log(tRes6),errtRes6/tRes6,marker='o',ls='',label='Resurrection',capsize=3)
plt.errorbar(J,np.log(tStall6d_13),errtStall6d_13/tStall6d_13,marker='^',ls='',label='Stall, $<n_0>=13$',capsize=3)
plt.errorbar(J,np.log(tStall6f_13),errtStall6f_13/tStall6f_13,marker='s',ls='',label='Stall, $n_0=13$',alpha=0.5,capsize=3)
plt.plot(cJ,np.log(15*t_th(cJ)),label=r'$log(\tau)\propto$ J')

plt.legend(prop={'size': 12})

plt.show()

plt.xlabel('J',fontsize=16)
plt.ylabel('$\Delta$',fontsize=16)
plt.title('$<N> = 6.5,N_{max}=13$ ')
plt.grid()
plt.errorbar(J,Dt6d_13,errDt6d_13,marker='^',ls='',label='$<n_0>=13$',color='darkorange',capsize=3)
plt.errorbar(J,Dt6f_13,errDt6d_13,marker='s',ls='',label='$n_0=13$',color='green',alpha=0.5,capsize=3)

plt.legend(prop={'size': 12})

plt.show()


#%%

plt.xlabel('J',fontsize=16)
plt.ylabel(r'$\tau$',fontsize=16)
plt.title('$<N> = 8,N_{max}=13$ ')
plt.grid()
plt.errorbar(J,np.log(tRes8),errtRes8/tRes8,marker='o',ls='',label='Resurrection',capsize=3)
plt.errorbar(J,np.log(tStall8d_13),errtStall8d_13/tStall8d_13,marker='^',ls='',label='Stall, $<n_0>=13$',capsize=3)
plt.errorbar(J,np.log(tStall8f_13),errtStall8f_13/tStall8f_13,marker='s',ls='',label='Stall, $n_0=13$',alpha=0.5,capsize=3)
plt.plot(cJ,np.log(15*t_th(cJ)),label=r'$log(\tau)\propto$ J')

plt.legend(prop={'size': 12})

plt.show()

plt.xlabel('J',fontsize=16)
plt.ylabel('$\Delta$',fontsize=16)
plt.title('$<N> = 8,N_{max}=13$ ')
plt.grid()
plt.errorbar(J,Dt8d_13,errDt8d_13,marker='^',ls='',label='$<n_0>=13$',color='darkorange',capsize=3)
plt.errorbar(J,Dt8f_13,errDt8d_13,marker='s',ls='',label='$n_0=13$',color='green',alpha=0.5,capsize=3)

plt.legend(prop={'size': 12})

plt.show()


#%% PLOTS OF FITS

plt.title(r'$<N_\infty>=4$, J=0', fontsize=16)
plt.xlabel('t (steps)',fontsize=14)
plt.ylabel('<N>',fontsize=14)
plt.grid()

plt.plot(Res4_J0[0:l-9700,0],Res4_J0[0:l-9700,1]*nmax,ls='dashed',lw=3,color='black',alpha=0.5,label='GMC res')
plt.plot(Stall4d_J0[0:l-9700,0],Stall4d_J0[0:l-9700,1]*nmax,ls='dashed',lw=3,color='darkorange',label='GMC rel ($<n_0>=6$)')
plt.plot(Stall4f_6_J0[0:l-9700,0],Stall4f_6_J0[0:l-9700,1]*nmax,ls='dashed',lw=3,color='green',label='GMC rel ($n_0=6$)')


plt.plot(t0,modelRes4_J0,color='red',label='Fit res')
plt.plot(t0,modelStall4d_J0,color='royalblue', label='Fit rel ($<n_0>=6$)')
plt.plot(t0,modelStall4f_6_J0,color='mediumvioletred',label='Fit rel ($n_0=6$)')

plt.legend(loc='lower right')

plt.show()

plt.title(r'$<N_\infty>=6.5$, J=0', fontsize=16)
plt.xlabel('t (steps)',fontsize=14)
plt.ylabel('<N>',fontsize=14)
plt.grid()

plt.plot(Res6_J0[0:l-9700,0],Res6_J0[0:l-9700,1]*nmax,ls='dashed',lw=3,color='black',alpha=0.5,label='GMC res')
plt.plot(Stall6d_9_J0[0:l-9700,0],Stall6d_9_J0[0:l-9700,1]*nmax,ls='dashed',lw=3,color='darkorange',label='GMC rel ($<n_0>=9$)')
plt.plot(Stall6f_8_J0[0:l-9700,0],Stall6f_8_J0[0:l-9700,1]*nmax,ls='dashed',lw=3,color='green',label='GMC rel ($n_0=9$)')


plt.plot(t0,modelRes6_J0,color='red',label='Fit res')
plt.plot(t0,modelStall6d_9_J0,color='royalblue', label='Fit rel ($<n_0>=9$)')
plt.plot(t0,modelStall6f_J0,color='mediumvioletred',label='Fit rel ($n_0=9$)')

plt.legend(loc='lower right')

plt.show()


plt.title(r'$<N_\infty>=8$, J=0', fontsize=16)
plt.xlabel('t (steps)',fontsize=14)
plt.ylabel('<N>',fontsize=14)
plt.grid()

plt.plot(Res8_J0[0:l-9700,0],Res8_J0[0:l-9700,1]*nmax,ls='dashed',lw=3,color='black',alpha=0.5,label='GMC res')
plt.plot(Stall8d_J0[0:l-9700,0],Stall8d_J0[0:l-9700,1]*nmax,ls='dashed',lw=3,color='darkorange',label='GMC rel ($<n_0>=9$)')
plt.plot(Stall8f_9_J0[0:l-9700,0],Stall8f_9_J0[0:l-9700,1]*nmax,ls='dashed',lw=3,color='green',label='GMC rel ($n_0=9$)')


plt.plot(t0,modelRes8_J0,color='red',label='Fit res')
plt.plot(t0,modelStall8d_J0,color='royalblue', label='Fit rel ($<n_0>=9$)')
plt.plot(t0,modelStall8f_9_J0,color='mediumvioletred',label='Fit rel ($n_0=9$)')

plt.legend(loc='lower right')

plt.show()

#%%

plt.title(r'$<N_\infty>=4$, J=3', fontsize=16)
plt.xlabel('t (steps)',fontsize=14)
plt.ylabel('<N>',fontsize=14)
plt.grid()

plt.plot(Res4_J3[0:l-8000,0],Res4_J3[0:l-8000,1]*nmax,ls='dashed',lw=3,color='black',alpha=0.5,label='GMC res')
plt.plot(Stall4d_J3[0:l-8000,0],Stall4d_J3[0:l-8000,1]*nmax,ls='dashed',lw=3,color='darkorange',label='GMC rel ($<n_0>=6$)')
plt.plot(Stall4f_6_J3[0:l-8000,0],Stall4f_6_J3[0:l-8000,1]*nmax,ls='dashed',lw=3,color='green',label='GMC rel ($n_0=6$)')


plt.plot(t3,modelRes4_J3,color='red',label='Fit res')
plt.plot(t3,modelStall4d_J3,color='royalblue', label='Fit rel ($<n_0>=6$)')
plt.plot(t3,modelStall4f_6_J3,color='mediumvioletred',label='Fit rel ($n_0=6$)')

plt.legend(loc='lower right')

plt.show()

plt.title(r'$<N_\infty>=6.5$, J=3', fontsize=16)
plt.xlabel('t (steps)',fontsize=14)
plt.ylabel('<N>',fontsize=14)
plt.grid()

plt.plot(Res6_J3[0:l-8000,0],Res6_J3[0:l-8000,1]*nmax,ls='dashed',lw=3,color='black',alpha=0.5,label='GMC res')
plt.plot(Stall6d_9_J3[0:l-8000,0],Stall6d_9_J3[0:l-8000,1]*nmax,ls='dashed',lw=3,color='darkorange',label='GMC rel ($<n_0>=9$)')
plt.plot(Stall6f_8_J3[0:l-8000,0],Stall6f_8_J3[0:l-8000,1]*nmax,ls='dashed',lw=3,color='green',label='GMC rel ($n_0=9$)')


plt.plot(t3,modelRes6_J3,color='red',label='Fit res')
plt.plot(t3,modelStall6d_9_J3,color='royalblue', label='Fit rel ($<n_0>=9$)')
plt.plot(t3,modelStall6f_J3,color='mediumvioletred',label='Fit rel ($n_0=9$)')

plt.legend(loc='lower right')

plt.show()


plt.title(r'$<N_\infty>=8$, J=3', fontsize=16)
plt.xlabel('t (steps)',fontsize=14)
plt.ylabel('<N>',fontsize=14)
plt.grid()

plt.plot(Res8_J3[0:l-8000,0],Res8_J3[0:l-8000,1]*nmax,ls='dashed',lw=3,color='black',alpha=0.5,label='GMC res')
plt.plot(Stall8d_J3[0:l-8000,0],Stall8d_J3[0:l-8000,1]*nmax,ls='dashed',lw=3,color='darkorange',label='GMC rel ($<n_0>=9$)')
plt.plot(Stall8f_9_J3[0:l-8000,0],Stall8f_9_J3[0:l-8000,1]*nmax,ls='dashed',lw=3,color='green',label='GMC rel ($n_0=9$)')


plt.plot(t3,modelRes8_J3,color='red',label='Fit res')
plt.plot(t3,modelStall8d_J3,color='royalblue', label='Fit rel ($<n_0>=9$)')
plt.plot(t3,modelStall8f_9_J3,color='mediumvioletred',label='Fit rel ($n_0=9$)')

plt.legend(loc='lower right')

plt.show()

#%%

plt.title(r'$<N_\infty>=4$, J=5', fontsize=16)
plt.xlabel('t (steps)',fontsize=14)
plt.ylabel('<N>',fontsize=14)
plt.grid()

plt.plot(Res4_J5[0:l-2000,0],Res4_J5[0:l-2000,1]*nmax,ls='dashed',lw=3,color='black',alpha=0.5,label='GMC res')
plt.plot(Stall4d_J5[0:l-2000,0],Stall4d_J5[0:l-2000,1]*nmax,ls='dashed',lw=3,color='darkorange',label='GMC rel ($<n_0>=6$)')
plt.plot(Stall4f_6_J5[0:l-2000,0],Stall4f_6_J5[0:l-2000,1]*nmax,ls='dashed',lw=3,color='green',label='GMC rel ($n_0=6$)')


plt.plot(t5,modelRes4_J5,color='red',label='Fit res')
plt.plot(t5,modelStall4d_J5,color='royalblue', label='Fit rel ($<n_0>=6$)')
plt.plot(t5,modelStall4f_6_J5,color='mediumvioletred',label='Fit rel ($n_0=6$)')

plt.legend(loc='lower right')

plt.show()

plt.title(r'$<N_\infty>=6.5$, J=5', fontsize=16)
plt.xlabel('t (steps)',fontsize=14)
plt.ylabel('<N>',fontsize=14)
plt.grid()

plt.plot(Res6_J5[0:l-2000,0],Res6_J5[0:l-2000,1]*nmax,ls='dashed',lw=3,color='black',alpha=0.5,label='GMC res')
plt.plot(Stall6d_9_J5[0:l-2000,0],Stall6d_9_J5[0:l-2000,1]*nmax,ls='dashed',lw=3,color='darkorange',label='GMC rel ($<n_0>=9$)')
plt.plot(Stall6f_8_J5[0:l-2000,0],Stall6f_8_J5[0:l-2000,1]*nmax,ls='dashed',lw=3,color='green',label='GMC rel ($n_0=9$)')


plt.plot(t5,modelRes6_J5,color='red',label='Fit res')
plt.plot(t5,modelStall6d_9_J5,color='royalblue', label='Fit rel ($<n_0>=9$)')
plt.plot(t5,modelStall6f_J5,color='mediumvioletred',label='Fit rel ($n_0=9$)')

plt.legend(loc='lower right')

plt.show()


plt.title(r'$<N_\infty>=8$, J=5', fontsize=16)
plt.xlabel('t (steps)',fontsize=14)
plt.ylabel('<N>',fontsize=14)
plt.grid()

plt.plot(Res8_J5[0:l-2000,0],Res8_J5[0:l-2000,1]*nmax,ls='dashed',lw=3,color='black',alpha=0.5,label='GMC res')
plt.plot(Stall8d_J5[0:l-2000,0],Stall8d_J5[0:l-2000,1]*nmax,ls='dashed',lw=3,color='darkorange',label='GMC rel ($<n_0>=9$)')
plt.plot(Stall8f_9_J5[0:l-2000,0],Stall8f_9_J5[0:l-2000,1]*nmax,ls='dashed',lw=3,color='green',label='GMC rel ($n_0=9$)')


plt.plot(t5,modelRes8_J5,color='red',label='Fit res')
plt.plot(t5,modelStall8d_J5,color='royalblue', label='Fit rel ($<n_0>=9$)')
plt.plot(t5,modelStall8f_9_J5,color='mediumvioletred',label='Fit rel ($n_0=9$)')

plt.legend(loc='lower right')

plt.show()

#%% PLOT + ERROR N=4

sumres4_J0=0.0
errsumres4_J0=0.0
for i in range(150,300):
    
    sumres4_J0 = sumres4_J0+Res4_J0[i,1]
    errsumres4_J0 = errsumres4_J0 + Res4_J0[i,2]
    
av4res_J0 = (sumres4_J0/150)*nmax
errav4res_J0 = (errsumres4_J0/150)*nmax
    
sumstall4_J0=0.0
errsumstall4_J0=0.0
for i in range(150,300):
    
    sumstall4_J0 = sumstall4_J0+Stall4d_J0[i,1]
    errsumstall4_J0 = errsumstall4_J0+Stall4d_J0[i,2]
    
av4stall_J0 = (sumstall4_J0/150)*nmax
erravstall4res_J0=(errsumstall4_J0/150)*nmax

av4_J0 = (av4res_J0+av4stall_J0)/2
errav4_J0 = (errav4res_J0+erravstall4res_J0)/2



sumres4_J3=0.0
errsumres4_J3=0.0
for i in range(1000,2000):
    
    sumres4_J3 = sumres4_J3+Res4_J3[i,1]
    errsumres4_J3 = errsumres4_J3 + Res4_J3[i,2]
    
av4res_J3 = (sumres4_J3/1000)*nmax
errav4res_J3 = (errsumres4_J3/1000)*nmax
    
sumstall4_J3=0.0
errsumstall4_J3=0.0
for i in range(1000,2000):
    
    sumstall4_J3 = sumstall4_J3+Stall4d_J3[i,1]
    errsumstall4_J3 = errsumstall4_J3+Stall4d_J3[i,2]
    
av4stall_J3 = (sumstall4_J3/1000)*nmax
erravstall4res_J3=(errsumstall4_J3/1000)*nmax

av4_J3 = (av4res_J3+av4stall_J3)/2
errav4_J3 = (errav4res_J3+erravstall4res_J3)/2


sumres4_J5=0.0
errsumres4_J5=0.0
for i in range(5000,8000):
    
    sumres4_J5 = sumres4_J5+Res4_J5[i,1]
    errsumres4_J5 = errsumres4_J5 + Res4_J5[i,2]
    
av4res_J5 = (sumres4_J5/3000)*nmax
errav4res_J5 = (errsumres4_J5/3000)*nmax
    
sumstall4_J5=0.0
errsumstall4_J5=0.0
for i in range(5000,8000):
    
    sumstall4_J5 = sumstall4_J5+Stall4d_J5[i,1]
    errsumstall4_J5 = errsumstall4_J5+Stall4d_J5[i,2]
    
av4stall_J5 = (sumstall4_J5/3000)*nmax
erravstall4res_J5=(errsumstall4_J5/3000)*nmax

av4_J5 = (av4res_J5+av4stall_J5)/2
errav4_J5 = (errav4res_J5+erravstall4res_J5)/2


plt.title(r'J=0 ,$<N_\infty>={:1.2f} \pm {:1.2f}$'.format(av4_J0,errav4_J0), fontsize=16)
plt.xlabel('t (steps)',fontsize=14)
plt.ylabel('<N>',fontsize=14)
plt.grid()

plt.plot(Res4_J0[0:l-9700,0],Res4_J0[0:l-9700,1]*nmax,ls='dashed',lw=3,color='black',alpha=0.5,label='GMC res')
plt.fill_between(Res4_J0[0:l-9700,0], Res4_J0[0:l-9700,1]*nmax-Res4_J0[0:l-9700,2]*nmax, 
                 Res4_J0[0:l-9700,1]*nmax+Res4_J0[0:l-9700,2]*nmax,alpha=0.3,color='black')

plt.plot(Stall4d_J0[0:l-9700,0],Stall4d_J0[0:l-9700,1]*nmax,ls='dashed',lw=3,color='darkorange',label='GMC rel ($<n_0>=6$)')
plt.fill_between(Stall4d_J0[0:l-9700,0], Stall4d_J0[0:l-9700,1]*nmax-Stall4d_J0[0:l-9700,2]*nmax, 
                 Stall4d_J0[0:l-9700,1]*nmax+Stall4d_J0[0:l-9700,2]*nmax,alpha=0.3,color='orange')

plt.show()


plt.title(r'J=3,$<N_\infty>={:1.2f} \pm {:1.2f}$'.format(av4_J3,errav4_J3), fontsize=16)
plt.xlabel('t (steps)',fontsize=14)
plt.ylabel('<N>',fontsize=14)
plt.grid()

plt.plot(Res4_J3[0:l-8000,0],Res4_J3[0:l-8000,1]*nmax,ls='dashed',lw=3,color='black',alpha=0.5,label='GMC res')
plt.fill_between(Res4_J3[0:l-8000,0], Res4_J3[0:l-8000,1]*nmax-Res4_J3[0:l-8000,2]*nmax, 
                 Res4_J3[0:l-8000,1]*nmax+Res4_J3[0:l-8000,2]*nmax,alpha=0.3,color='black')


plt.plot(Stall4d_J3[0:l-8000,0],Stall4d_J3[0:l-8000,1]*nmax,ls='dashed',lw=3,color='darkorange',label='GMC rel ($<n_0>=6$)')
plt.fill_between(Stall4d_J3[0:l-8000,0], Stall4d_J3[0:l-8000,1]*nmax-Stall4d_J3[0:l-8000,2]*nmax, 
                 Stall4d_J3[0:l-8000,1]*nmax+Stall4d_J3[0:l-8000,2]*nmax,alpha=0.3,color='orange')


plt.show()

plt.title(r'J=5,$<N_\infty>={:1.2f} \pm {:1.2f}$'.format(av4_J5,errav4_J5), fontsize=16)
plt.xlabel('t (steps)',fontsize=14)
plt.ylabel('<N>',fontsize=14)
plt.grid()

plt.plot(Res4_J5[0:l-2000,0],Res4_J5[0:l-2000,1]*nmax,ls='dashed',lw=3,color='black',alpha=0.5,label='GMC res')
plt.fill_between(Res4_J5[0:l-2000,0], Res4_J5[0:l-2000,1]*nmax-Res4_J5[0:l-2000,2]*nmax, 
                 Res4_J5[0:l-2000,1]*nmax+Res4_J5[0:l-2000,2]*nmax,alpha=0.3,color='black')

plt.plot(Stall4d_J5[0:l-2000,0],Stall4d_J5[0:l-2000,1]*nmax,ls='dashed',lw=3,color='darkorange',label='GMC rel ($<n_0>=6$)')
plt.fill_between(Stall4d_J5[0:l-2000,0], Stall4d_J5[0:l-2000,1]*nmax-Stall4d_J5[0:l-2000,2]*nmax, 
                 Stall4d_J5[0:l-2000,1]*nmax+Stall4d_J5[0:l-2000,2]*nmax,alpha=0.3,color='orange')

#%%

sumres6_J0=0.0
errsumres6_J0=0.0
for i in range(150,300):
    
    sumres6_J0 = sumres6_J0+Res6_J0[i,1]
    errsumres6_J0 = errsumres6_J0 + Res6_J0[i,2]
    
av6res_J0 = (sumres6_J0/150)*nmax
errav6res_J0 = (errsumres6_J0/150)*nmax
    
sumstall6_J0=0.0
errsumstall6_J0=0.0
for i in range(150,300):
    
    sumstall6_J0 = sumstall6_J0+Stall6d_9_J0[i,1]
    errsumstall6_J0 = errsumstall6_J0+Stall6d_9_J0[i,2]
    
av6stall_J0 = (sumstall6_J0/150)*nmax
erravstall6res_J0=(errsumstall6_J0/150)*nmax

av6_J0 = (av6res_J0+av6stall_J0)/2
errav6_J0 = (errav6res_J0+erravstall6res_J0)/2



sumres6_J3=0.0
errsumres6_J3=0.0
for i in range(1000,2000):
    
    sumres6_J3 = sumres6_J3+Res6_J3[i,1]
    errsumres4_J3 = errsumres6_J3 + Res6_J3[i,2]
    
av6res_J3 = (sumres6_J3/1000)*nmax
errav6res_J3 = (errsumres6_J3/1000)*nmax
    
sumstall6_J3=0.0
errsumstall6_J3=0.0
for i in range(1000,2000):
    
    sumstall6_J3 = sumstall6_J3+Stall6d_9_J3[i,1]
    errsumstall6_J3 = errsumstall6_J3+Stall6d_9_J3[i,2]
    
av6stall_J3 = (sumstall6_J3/1000)*nmax
erravstall6res_J3=(errsumstall6_J3/1000)*nmax

av6_J3 = (av6res_J3+av6stall_J3)/2
errav6_J3 = (errav6res_J3+erravstall6res_J3)/2


sumres6_J5=0.0
errsumres6_J5=0.0
for i in range(5000,8000):
    
    sumres6_J5 = sumres6_J5+Res6_J5[i,1]
    errsumres6_J5 = errsumres6_J5 + Res6_J5[i,2]
    
av6res_J5 = (sumres6_J5/3000)*nmax
errav6res_J5 = (errsumres6_J5/3000)*nmax
    
sumstall6_J5=0.0
errsumstall6_J5=0.0
for i in range(5000,8000):
    
    sumstall6_J5 = sumstall6_J5+Stall6d_9_J5[i,1]
    errsumstall6_J5 = errsumstall6_J5+Stall6d_9_J5[i,2]
    
av6stall_J5 = (sumstall6_J5/3000)*nmax
erravstall6res_J5=(errsumstall6_J5/3000)*nmax

av6_J5 = (av6res_J5+av6stall_J5)/2
errav6_J5 = (errav6res_J5+erravstall6res_J5)/2


plt.title(r'J=0 ,$<N_\infty>={:1.2f} \pm {:1.2f}$'.format(av6_J0,errav6_J0), fontsize=16)
plt.xlabel('t (steps)',fontsize=14)
plt.ylabel('<N>',fontsize=14)
plt.grid()

plt.plot(Res6_J0[0:l-9700,0],Res6_J0[0:l-9700,1]*nmax,ls='dashed',lw=3,color='black',alpha=0.5,label='GMC res')
plt.fill_between(Res6_J0[0:l-9700,0], Res6_J0[0:l-9700,1]*nmax-Res6_J0[0:l-9700,2]*nmax, 
                 Res6_J0[0:l-9700,1]*nmax+Res6_J0[0:l-9700,2]*nmax,alpha=0.3,color='black')

plt.plot(Stall6d_9_J0[0:l-9700,0],Stall6d_9_J0[0:l-9700,1]*nmax,ls='dashed',lw=3,color='darkorange',label='GMC rel ($<n_0>=6$)')
plt.fill_between(Stall6d_9_J0[0:l-9700,0], Stall6d_9_J0[0:l-9700,1]*nmax-Stall6d_9_J0[0:l-9700,2]*nmax, 
                 Stall6d_9_J0[0:l-9700,1]*nmax+Stall6d_9_J0[0:l-9700,2]*nmax,alpha=0.3,color='orange')

plt.show()


plt.title(r'J=3,$<N_\infty>={:1.2f} \pm {:1.2f}$'.format(av6_J3,errav6_J3), fontsize=16)
plt.xlabel('t (steps)',fontsize=14)
plt.ylabel('<N>',fontsize=14)
plt.grid()

plt.plot(Res6_J3[0:l-8000,0],Res6_J3[0:l-8000,1]*nmax,ls='dashed',lw=3,color='black',alpha=0.5,label='GMC res')
plt.fill_between(Res6_J3[0:l-8000,0], Res6_J3[0:l-8000,1]*nmax-Res6_J3[0:l-8000,2]*nmax, 
                 Res6_J3[0:l-8000,1]*nmax+Res6_J3[0:l-8000,2]*nmax,alpha=0.3,color='black')


plt.plot(Stall6d_9_J3[0:l-8000,0],Stall6d_9_J3[0:l-8000,1]*nmax,ls='dashed',lw=3,color='darkorange',label='GMC rel ($<n_0>=6$)')
plt.fill_between(Stall6d_9_J3[0:l-8000,0], Stall6d_9_J3[0:l-8000,1]*nmax-Stall6d_9_J3[0:l-8000,2]*nmax, 
                 Stall6d_9_J3[0:l-8000,1]*nmax+Stall6d_9_J3[0:l-8000,2]*nmax,alpha=0.3,color='orange')


plt.show()

plt.title(r'J=5,$<N_\infty>={:1.2f} \pm {:1.2f}$'.format(av6_J5,errav6_J5), fontsize=16)
plt.xlabel('t (steps)',fontsize=14)
plt.ylabel('<N>',fontsize=14)
plt.grid()

plt.plot(Res6_J5[0:l-2000,0],Res6_J5[0:l-2000,1]*nmax,ls='dashed',lw=3,color='black',alpha=0.5,label='GMC res')
plt.fill_between(Res6_J5[0:l-2000,0], Res6_J5[0:l-2000,1]*nmax-Res6_J5[0:l-2000,2]*nmax, 
                 Res6_J5[0:l-2000,1]*nmax+Res6_J5[0:l-2000,2]*nmax,alpha=0.3,color='black')

plt.plot(Stall6d_9_J5[0:l-2000,0],Stall6d_9_J5[0:l-2000,1]*nmax,ls='dashed',lw=3,color='darkorange',label='GMC rel ($<n_0>=6$)')
plt.fill_between(Stall6d_9_J5[0:l-2000,0], Stall6d_9_J5[0:l-2000,1]*nmax-Stall6d_9_J5[0:l-2000,2]*nmax, 
                 Stall6d_9_J5[0:l-2000,1]*nmax+Stall6d_9_J5[0:l-2000,2]*nmax,alpha=0.3,color='orange')

#%%

sumres8_J0=0.0
errsumres8_J0=0.0
for i in range(150,300):
    
    sumres8_J0 = sumres8_J0+Res8_J0[i,1]
    errsumres8_J0 = errsumres8_J0 + Res8_J0[i,2]
    
av8res_J0 = (sumres8_J0/150)*nmax
errav8res_J0 = (errsumres8_J0/150)*nmax
    
sumstall8_J0=0.0
errsumstall8_J0=0.0
for i in range(150,300):
    
    sumstall8_J0 = sumstall8_J0+Stall8d_J0[i,1]
    errsumstall8_J0 = errsumstall8_J0+Stall8d_J0[i,2]
    
av8stall_J0 = (sumstall8_J0/150)*nmax
erravstall8res_J0=(errsumstall8_J0/150)*nmax

av8_J0 = (av8res_J0+av8stall_J0)/2
errav8_J0 = (errav8res_J0+erravstall8res_J0)/2



sumres8_J3=0.0
errsumres8_J3=0.0
for i in range(1000,2000):
    
    sumres8_J3 = sumres8_J3+Res8_J3[i,1]
    errsumres8_J3 = errsumres8_J3 + Res8_J3[i,2]
    
av8res_J3 = (sumres8_J3/1000)*nmax
errav8res_J3 = (errsumres8_J3/1000)*nmax
    
sumstall8_J3=0.0
errsumstall8_J3=0.0
for i in range(1000,2000):
    
    sumstall8_J3 = sumstall8_J3+Stall8d_J3[i,1]
    errsumstall8_J3 = errsumstall8_J3+Stall8d_J3[i,2]
    
av8stall_J3 = (sumstall8_J3/1000)*nmax
erravstall8res_J3=(errsumstall8_J3/1000)*nmax

av8_J3 = (av8res_J3+av8stall_J3)/2
errav8_J3 = (errav8res_J3+erravstall8res_J3)/2


sumres8_J5=0.0
errsumres8_J5=0.0
for i in range(5000,8000):
    
    sumres8_J5 = sumres8_J5+Res8_J5[i,1]
    errsumres8_J5 = errsumres8_J5 + Res8_J5[i,2]
    
av8res_J5 = (sumres8_J5/3000)*nmax
errav8res_J5 = (errsumres8_J5/3000)*nmax
    
sumstall8_J5=0.0
errsumstall8_J5=0.0
for i in range(5000,8000):
    
    sumstall8_J5 = sumstall8_J5+Stall8d_J5[i,1]
    errsumstall8_J5 = errsumstall8_J5+Stall8d_J5[i,2]
    
av8stall_J5 = (sumstall8_J5/3000)*nmax
erravstall8res_J5=(errsumstall8_J5/3000)*nmax

av8_J5 = (av8res_J5+av8stall_J5)/2
errav8_J5 = (errav8res_J5+erravstall8res_J5)/2


plt.title(r'J=0 ,$<N_\infty>={:1.2f} \pm {:1.2f}$'.format(av8_J0,errav8_J0), fontsize=16)
plt.xlabel('t (steps)',fontsize=14)
plt.ylabel('<N>',fontsize=14)
plt.grid()

plt.plot(Res8_J0[0:l-9700,0],Res8_J0[0:l-9700,1]*nmax,ls='dashed',lw=3,color='black',alpha=0.5,label='GMC res')
plt.fill_between(Res8_J0[0:l-9700,0], Res8_J0[0:l-9700,1]*nmax-Res8_J0[0:l-9700,2]*nmax, 
                 Res8_J0[0:l-9700,1]*nmax+Res8_J0[0:l-9700,2]*nmax,alpha=0.3,color='black')

plt.plot(Stall8d_J0[0:l-9700,0],Stall8d_J0[0:l-9700,1]*nmax,ls='dashed',lw=3,color='darkorange',label='GMC rel ($<n_0>=6$)')
plt.fill_between(Stall8d_J0[0:l-9700,0], Stall8d_J0[0:l-9700,1]*nmax-Stall8d_J0[0:l-9700,2]*nmax, 
                 Stall8d_J0[0:l-9700,1]*nmax+Stall8d_J0[0:l-9700,2]*nmax,alpha=0.3,color='orange')

plt.show()


plt.title(r'J=3,$<N_\infty>={:1.2f} \pm {:1.2f}$'.format(av8_J3,errav8_J3), fontsize=16)
plt.xlabel('t (steps)',fontsize=14)
plt.ylabel('<N>',fontsize=14)
plt.grid()

plt.plot(Res8_J3[0:l-8000,0],Res8_J3[0:l-8000,1]*nmax,ls='dashed',lw=3,color='black',alpha=0.5,label='GMC res')
plt.fill_between(Res8_J3[0:l-8000,0], Res8_J3[0:l-8000,1]*nmax-Res8_J3[0:l-8000,2]*nmax, 
                 Res8_J3[0:l-8000,1]*nmax+Res8_J3[0:l-8000,2]*nmax,alpha=0.3,color='black')


plt.plot(Stall8d_J3[0:l-8000,0],Stall8d_J3[0:l-8000,1]*nmax,ls='dashed',lw=3,color='darkorange',label='GMC rel ($<n_0>=6$)')
plt.fill_between(Stall8d_J3[0:l-8000,0], Stall8d_J3[0:l-8000,1]*nmax-Stall8d_J3[0:l-8000,2]*nmax, 
                 Stall8d_J3[0:l-8000,1]*nmax+Stall8d_J3[0:l-8000,2]*nmax,alpha=0.3,color='orange')


plt.show()

plt.title(r'J=5,$<N_\infty>={:1.2f} \pm {:1.2f}$'.format(av8_J5,errav8_J5), fontsize=16)
plt.xlabel('t (steps)',fontsize=14)
plt.ylabel('<N>',fontsize=14)
plt.grid()

plt.plot(Res8_J5[0:l-2000,0],Res8_J5[0:l-2000,1]*nmax,ls='dashed',lw=3,color='black',alpha=0.5,label='GMC res')
plt.fill_between(Res8_J5[0:l-2000,0], Res8_J5[0:l-2000,1]*nmax-Res8_J5[0:l-2000,2]*nmax, 
                 Res8_J5[0:l-2000,1]*nmax+Res8_J5[0:l-2000,2]*nmax,alpha=0.3,color='black')

plt.plot(Stall8d_J5[0:l-2000,0],Stall8d_J5[0:l-2000,1]*nmax,ls='dashed',lw=3,color='darkorange',label='GMC rel ($<n_0>=6$)')
plt.fill_between(Stall8d_J5[0:l-2000,0], Stall8d_J5[0:l-2000,1]*nmax-Stall8d_J5[0:l-2000,2]*nmax, 
                 Stall8d_J5[0:l-2000,1]*nmax+Stall8d_J5[0:l-2000,2]*nmax,alpha=0.3,color='orange')
