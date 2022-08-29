#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 10:28:07 2022

@author: luke
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.lines as mlines

nmax=13

#Function used for fits
def N(x,tau,nss):
    return nss + (n0-nss)*np.exp(-x/tau) 

l=100000

#%%

files0Res_KMC=['/media/luke/Elements/Simulations/Without _depletion/KMC/CooperativityKMC/Traces_Res_KMC_Coop_Berg_J=0.00_N=4.dat']
data0Res_KMC=[]
files0Stall_KMC=['/media/luke/Elements/Simulations/Without _depletion/KMC/CooperativityKMC/Traces_Stall_KMC_Coop_Berg_J=0.00_N=4.dat']
data0Stall_KMC=[]

for data_file in files0Res_KMC:
    data0Res_KMC.append(np.loadtxt(data_file))
for data_file in files0Stall_KMC:
    data0Stall_KMC.append(np.loadtxt(data_file))
    
Res0_KMC = data0Res_KMC[0]
Stall0_KMC = data0Stall_KMC[0]


files1Res_KMC=['/media/luke/Elements/Simulations/Without _depletion/KMC/CooperativityKMC/Traces_Res_KMC_Coop_Berg_J=1.00_N=4.dat']
data1Res_KMC=[]
files1Stall_KMC=['/media/luke/Elements/Simulations/Without _depletion/KMC/CooperativityKMC/Traces_Stall_KMC_Coop_Berg_J=1.00_N=4.dat']
data1Stall_KMC=[]

for data_file in files1Res_KMC:
    data1Res_KMC.append(np.loadtxt(data_file))
for data_file in files1Stall_KMC:
    data1Stall_KMC.append(np.loadtxt(data_file))
    
Res1_KMC = data1Res_KMC[0]
Stall1_KMC = data1Stall_KMC[0]


files2Res_KMC=['/media/luke/Elements/Simulations/Without _depletion/KMC/CooperativityKMC/Traces_Res_KMC_Coop_Berg_J=2.00_N=4.dat']
data2Res_KMC=[]
files2Stall_KMC=['/media/luke/Elements/Simulations/Without _depletion/KMC/CooperativityKMC/Traces_Stall_KMC_Coop_Berg_J=2.00_N=4.dat']
data2Stall_KMC=[]

for data_file in files2Res_KMC:
    data2Res_KMC.append(np.loadtxt(data_file))
for data_file in files2Stall_KMC:
    data2Stall_KMC.append(np.loadtxt(data_file))
    
Res2_KMC = data2Res_KMC[0]
Stall2_KMC = data2Stall_KMC[0]


files3Res_KMC=['/media/luke/Elements/Simulations/Without _depletion/KMC/CooperativityKMC/Traces_Res_KMC_Coop_Berg_J=3.00_N=4.dat']
data3Res_KMC=[]
files3Stall_KMC=['/media/luke/Elements/Simulations/Without _depletion/KMC/CooperativityKMC/Traces_Stall_KMC_Coop_Berg_J=3.00_N=4.dat']
data3Stall_KMC=[]

for data_file in files3Res_KMC:
    data3Res_KMC.append(np.loadtxt(data_file))
for data_file in files3Stall_KMC:
    data3Stall_KMC.append(np.loadtxt(data_file))
    
Res3_KMC = data3Res_KMC[0]
Stall3_KMC = data3Stall_KMC[0]


files4Res_KMC=['/media/luke/Elements/Simulations/Without _depletion/KMC/CooperativityKMC/Traces_Res_KMC_Coop_Berg_J=4.00_N=4.dat']
data4Res_KMC=[]
files4Stall_KMC=['/media/luke/Elements/Simulations/Without _depletion/KMC/CooperativityKMC/Traces_Stall_KMC_Coop_Berg_J=4.00_N=4.dat']
data4Stall_KMC=[]

for data_file in files4Res_KMC:
    data4Res_KMC.append(np.loadtxt(data_file))
for data_file in files4Stall_KMC:
    data4Stall_KMC.append(np.loadtxt(data_file))
    
Res4_KMC = data4Res_KMC[0]
Stall4_KMC = data4Stall_KMC[0]


files5Res_KMC=['/media/luke/Elements/Simulations/Without _depletion/KMC/CooperativityKMC/Traces_Res_KMC_Coop_Berg_J=5.00_N=4.dat']
data5Res_KMC=[]
files5Stall_KMC=['/media/luke/Elements/Simulations/Without _depletion/KMC/CooperativityKMC/Traces_Stall_KMC_Coop_Berg_J=5.00_N=4.dat']
data5Stall_KMC=[]

for data_file in files5Res_KMC:
    data5Res_KMC.append(np.loadtxt(data_file))
for data_file in files5Stall_KMC:
    data5Stall_KMC.append(np.loadtxt(data_file))
    
Res5_KMC = data5Res_KMC[0]
Stall5_KMC = data5Stall_KMC[0]




files0Res_G=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Traces_Res_J=0.00_N=4.dat'] #All resurrection simulations
data0Res_G=[]
files0Stall_G=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Traces_Stall_J=0.00_N=4.dat'] #All resurrection simulations
data0Stall_G=[]

for data_file in files0Res_G:
    data0Res_G.append(np.loadtxt(data_file))
for data_file in files0Stall_G:
    data0Stall_G.append(np.loadtxt(data_file))
    
Res0_G=data0Res_G[0]
Stall0_G=data0Stall_G[0]


files1Res_G=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Traces_Res_J=1.00_N=4.dat'] #All resurrection simulations
data1Res_G=[]
files1Stall_G=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Traces_Stall_J=1.00_N=4.dat'] #All resurrection simulations
data1Stall_G=[]

for data_file in files1Res_G:
    data1Res_G.append(np.loadtxt(data_file))
for data_file in files1Stall_G:
    data1Stall_G.append(np.loadtxt(data_file))
    
Res1_G=data1Res_G[0]
Stall1_G=data1Stall_G[0]


files2Res_G=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Traces_Res_J=2.00_N=4.dat'] #All resurrection simulations
data2Res_G=[]
files2Stall_G=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Traces_Stall_J=2.00_N=4.dat'] #All resurrection simulations
data2Stall_G=[]

for data_file in files2Res_G:
    data2Res_G.append(np.loadtxt(data_file))
for data_file in files2Stall_G:
    data2Stall_G.append(np.loadtxt(data_file))
    
Res2_G=data2Res_G[0]
Stall2_G=data2Stall_G[0]


files3Res_G=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Traces_Res_J=3.00_N=4.dat'] #All resurrection simulations
data3Res_G=[]
files3Stall_G=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Traces_Stall_J=3.00_N=4.dat'] #All resurrection simulations
data3Stall_G=[]

for data_file in files3Res_G:
    data3Res_G.append(np.loadtxt(data_file))
for data_file in files3Stall_G:
    data3Stall_G.append(np.loadtxt(data_file))
    
Res3_G=data3Res_G[0]
Stall3_G=data3Stall_G[0]


files4Res_G=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Traces_Res_J=4.00_N=4.dat'] #All resurrection simulations
data4Res_G=[]
files4Stall_G=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Traces_Stall_J=4.00_N=4.dat'] #All resurrection simulations
data4Stall_G=[]

for data_file in files4Res_G:
    data4Res_G.append(np.loadtxt(data_file))
for data_file in files4Stall_G:
    data4Stall_G.append(np.loadtxt(data_file))
    
Res4_G=data4Res_G[0]
Stall4_G=data4Stall_G[0]


files5Res_G=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Traces_Res_J=5.00_N=4.dat'] #All resurrection simulations
data5Res_G=[]
files5Stall_G=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Traces_Stall_J=5.00_N=4.dat'] #All resurrection simulations
data5Stall_G=[]

for data_file in files5Res_G:
    data5Res_G.append(np.loadtxt(data_file))
for data_file in files5Stall_G:
    data5Stall_G.append(np.loadtxt(data_file))
    
Res5_G=data5Res_G[0]
Stall5_G=data5Stall_G[0]

#%%

nt=1000

r0= Res0_KMC[nt,0]/Res0_G[nt,0]
r1= Res1_KMC[nt,0]/Res1_G[nt,0]
r2= Res2_KMC[nt,0]/Res2_G[nt,0]
r3= Res3_KMC[nt,0]/Res3_G[nt,0]
r4= Res4_KMC[nt,0]/Res4_G[nt,0]

r0

plt.plot(r0*Res0_G[0:500,0],Res0_G[0:500,1]*nmax)
plt.plot(r0*Stall0_G[0:500,0],Stall0_G[0:500,1]*nmax)

plt.plot(Res0_KMC[0:500,0],Res0_KMC[0:500,1],ls='dashed')
plt.plot(Stall0_KMC[0:500,0],Stall0_KMC[0:500,1],ls='dashed')

plt.show()


plt.plot(r1*Res1_G[0:500,0],Res1_G[0:500,1]*nmax)
plt.plot(r1*Stall1_G[0:500,0],Stall1_G[0:500,1]*nmax)

plt.plot(Res1_KMC[0:500,0],Res1_KMC[0:500,1],ls='dashed')
plt.plot(Stall1_KMC[0:500,0],Stall1_KMC[0:500,1],ls='dashed')

plt.show()


plt.plot(r2*Res2_G[0:1000,0],Res2_G[0:1000,1]*nmax)
plt.plot(r2*Stall2_G[0:1000,0],Stall2_G[0:1000,1]*nmax)

plt.plot(Res2_KMC[0:1000,0],Res2_KMC[0:1000,1],ls='dashed')
plt.plot(Stall2_KMC[0:1000,0],Stall2_KMC[0:1000,1],ls='dashed')

plt.show()


plt.plot(r3*Res3_G[0:1500,0],Res3_G[0:1500,1]*nmax)
plt.plot(r3*Stall3_G[0:1500,0],Stall3_G[0:1500,1]*nmax)

plt.plot(Res3_KMC[0:1500,0],Res3_KMC[0:1500,1],ls='dashed')
plt.plot(Stall3_KMC[0:1500,0],Stall3_KMC[0:1500,1],ls='dashed')

plt.show()

#%%


###########################################################################################################

                                                #N=8#

###########################################################################################################


files0Res_KMC=['/media/luke/Elements/Simulations/Without _depletion/KMC/CooperativityKMC/Traces_Res_KMC_Coop_Berg_J=0.00_N=8.dat']
data0Res_KMC=[]
files0Stall_KMC=['/media/luke/Elements/Simulations/Without _depletion/KMC/CooperativityKMC/Traces_Stall_KMC_Coop_Berg_J=0.00_N=8.dat']
data0Stall_KMC=[]

for data_file in files0Res_KMC:
    data0Res_KMC.append(np.loadtxt(data_file))
for data_file in files0Stall_KMC:
    data0Stall_KMC.append(np.loadtxt(data_file))
    
Res0_KMC = data0Res_KMC[0]
Stall0_KMC = data0Stall_KMC[0]


files1Res_KMC=['/media/luke/Elements/Simulations/Without _depletion/KMC/CooperativityKMC/Traces_Res_KMC_Coop_Berg_J=1.00_N=8.dat']
data1Res_KMC=[]
files1Stall_KMC=['/media/luke/Elements/Simulations/Without _depletion/KMC/CooperativityKMC/Traces_Stall_KMC_Coop_Berg_J=1.00_N=8.dat']
data1Stall_KMC=[]

for data_file in files1Res_KMC:
    data1Res_KMC.append(np.loadtxt(data_file))
for data_file in files1Stall_KMC:
    data1Stall_KMC.append(np.loadtxt(data_file))
    
Res1_KMC = data1Res_KMC[0]
Stall1_KMC = data1Stall_KMC[0]


files2Res_KMC=['/media/luke/Elements/Simulations/Without _depletion/KMC/CooperativityKMC/Traces_Res_KMC_Coop_Berg_J=2.00_N=8.dat']
data2Res_KMC=[]
files2Stall_KMC=['/media/luke/Elements/Simulations/Without _depletion/KMC/CooperativityKMC/Traces_Stall_KMC_Coop_Berg_J=2.00_N=8.dat']
data2Stall_KMC=[]

for data_file in files2Res_KMC:
    data2Res_KMC.append(np.loadtxt(data_file))
for data_file in files2Stall_KMC:
    data2Stall_KMC.append(np.loadtxt(data_file))
    
Res2_KMC = data2Res_KMC[0]
Stall2_KMC = data2Stall_KMC[0]


files3Res_KMC=['/media/luke/Elements/Simulations/Without _depletion/KMC/CooperativityKMC/Traces_Res_KMC_Coop_Berg_J=3.00_N=8.dat']
data3Res_KMC=[]
files3Stall_KMC=['/media/luke/Elements/Simulations/Without _depletion/KMC/CooperativityKMC/Traces_Stall_KMC_Coop_Berg_J=3.00_N=8.dat']
data3Stall_KMC=[]

for data_file in files3Res_KMC:
    data3Res_KMC.append(np.loadtxt(data_file))
for data_file in files3Stall_KMC:
    data3Stall_KMC.append(np.loadtxt(data_file))
    
Res3_KMC = data3Res_KMC[0]
Stall3_KMC = data3Stall_KMC[0]


files4Res_KMC=['/media/luke/Elements/Simulations/Without _depletion/KMC/CooperativityKMC/Traces_Res_KMC_Coop_Berg_J=4.00_N=8.dat']
data4Res_KMC=[]
files4Stall_KMC=['/media/luke/Elements/Simulations/Without _depletion/KMC/CooperativityKMC/Traces_Stall_KMC_Coop_Berg_J=4.00_N=8.dat']
data4Stall_KMC=[]

for data_file in files4Res_KMC:
    data4Res_KMC.append(np.loadtxt(data_file))
for data_file in files4Stall_KMC:
    data4Stall_KMC.append(np.loadtxt(data_file))
    
Res4_KMC = data4Res_KMC[0]
Stall4_KMC = data4Stall_KMC[0]


files5Res_KMC=['/media/luke/Elements/Simulations/Without _depletion/KMC/CooperativityKMC/Traces_Res_KMC_Coop_Berg_J=5.00_N=8.dat']
data5Res_KMC=[]
files5Stall_KMC=['/media/luke/Elements/Simulations/Without _depletion/KMC/CooperativityKMC/Traces_Stall_KMC_Coop_Berg_J=5.00_N=8.dat']
data5Stall_KMC=[]

for data_file in files5Res_KMC:
    data5Res_KMC.append(np.loadtxt(data_file))
for data_file in files5Stall_KMC:
    data5Stall_KMC.append(np.loadtxt(data_file))
    
Res5_KMC = data5Res_KMC[0]
Stall5_KMC = data5Stall_KMC[0]




files0Res_G=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Traces_Res_J=0.00_N=4.dat'] #All resurrection simulations
data0Res_G=[]
files0Stall_G=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Traces_Stall_J=0.00_N=4.dat'] #All resurrection simulations
data0Stall_G=[]

for data_file in files0Res_G:
    data0Res_G.append(np.loadtxt(data_file))
for data_file in files0Stall_G:
    data0Stall_G.append(np.loadtxt(data_file))
    
Res0_G=data0Res_G[0]
Stall0_G=data0Stall_G[0]


files1Res_G=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Traces_Res_J=1.00_N=4.dat'] #All resurrection simulations
data1Res_G=[]
files1Stall_G=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Traces_Stall_J=1.00_N=4.dat'] #All resurrection simulations
data1Stall_G=[]

for data_file in files1Res_G:
    data1Res_G.append(np.loadtxt(data_file))
for data_file in files1Stall_G:
    data1Stall_G.append(np.loadtxt(data_file))
    
Res1_G=data1Res_G[0]
Stall1_G=data1Stall_G[0]


files2Res_G=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Traces_Res_J=2.00_N=4.dat'] #All resurrection simulations
data2Res_G=[]
files2Stall_G=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Traces_Stall_J=2.00_N=4.dat'] #All resurrection simulations
data2Stall_G=[]

for data_file in files2Res_G:
    data2Res_G.append(np.loadtxt(data_file))
for data_file in files2Stall_G:
    data2Stall_G.append(np.loadtxt(data_file))
    
Res2_G=data2Res_G[0]
Stall2_G=data2Stall_G[0]


files3Res_G=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Traces_Res_J=3.00_N=4.dat'] #All resurrection simulations
data3Res_G=[]
files3Stall_G=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Traces_Stall_J=3.00_N=4.dat'] #All resurrection simulations
data3Stall_G=[]

for data_file in files3Res_G:
    data3Res_G.append(np.loadtxt(data_file))
for data_file in files3Stall_G:
    data3Stall_G.append(np.loadtxt(data_file))
    
Res3_G=data3Res_G[0]
Stall3_G=data3Stall_G[0]


files4Res_G=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Traces_Res_J=4.00_N=8.dat'] #All resurrection simulations
data4Res_G=[]
files4Stall_G=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Traces_Stall_J=4.00_N=8.dat'] #All resurrection simulations
data4Stall_G=[]

for data_file in files4Res_G:
    data4Res_G.append(np.loadtxt(data_file))
for data_file in files4Stall_G:
    data4Stall_G.append(np.loadtxt(data_file))
    
Res4_G=data4Res_G[0]
Stall4_G=data4Stall_G[0]


files5Res_G=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Traces_Res_J=5.00_N=8.dat'] #All resurrection simulations
data5Res_G=[]
files5Stall_G=['/media/luke/Elements/Simulations/Without _depletion/Glauber/Traces_Stall_J=5.00_N=8.dat'] #All resurrection simulations
data5Stall_G=[]

for data_file in files5Res_G:
    data5Res_G.append(np.loadtxt(data_file))
for data_file in files5Stall_G:
    data5Stall_G.append(np.loadtxt(data_file))
    
Res5_G=data5Res_G[0]
Stall5_G=data5Stall_G[0]


#%%

nt=1000

r0= Res0_KMC[nt,0]/Res0_G[nt,0]
r1= Res1_KMC[nt,0]/Res1_G[nt,0]
r2= Res2_KMC[nt,0]/Res2_G[nt,0]
r3= Res3_KMC[nt,0]/Res3_G[nt,0]
r4= Res4_KMC[nt,0]/Res4_G[nt,0]
r5= Res5_KMC[nt,0]/Res5_G[nt,0]

r0

plt.plot(r0*Res0_G[0:500,0],Res0_G[0:500,1]*nmax)
plt.plot(r0*Stall0_G[0:500,0],Stall0_G[0:500,1]*nmax)

plt.plot(Res0_KMC[0:500,0],Res0_KMC[0:500,1],ls='dashed')
plt.plot(Stall0_KMC[0:500,0],Stall0_KMC[0:500,1],ls='dashed')

plt.show()


plt.plot(r1*Res1_G[0:500,0],Res1_G[0:500,1]*nmax)
plt.plot(r1*Stall1_G[0:500,0],Stall1_G[0:500,1]*nmax)

plt.plot(Res1_KMC[0:500,0],Res1_KMC[0:500,1],ls='dashed')
plt.plot(Stall1_KMC[0:500,0],Stall1_KMC[0:500,1],ls='dashed')

plt.show()


plt.plot(r2*Res2_G[0:1000,0],Res2_G[0:1000,1]*nmax)
plt.plot(r2*Stall2_G[0:1000,0],Stall2_G[0:1000,1]*nmax)

plt.plot(Res2_KMC[0:1000,0],Res2_KMC[0:1000,1],ls='dashed')
plt.plot(Stall2_KMC[0:1000,0],Stall2_KMC[0:1000,1],ls='dashed')

plt.show()


plt.plot(r3*Res3_G[0:1500,0],Res3_G[0:1500,1]*nmax)
plt.plot(r3*Stall3_G[0:1500,0],Stall3_G[0:1500,1]*nmax)

plt.plot(Res3_KMC[0:1500,0],Res3_KMC[0:1500,1],ls='dashed')
plt.plot(Stall3_KMC[0:1500,0],Stall3_KMC[0:1500,1],ls='dashed')

plt.show()


plt.plot(r4*Res4_G[0:2500,0],Res4_G[0:2500,1]*nmax)
plt.plot(r4*Stall4_G[0:2500,0],Stall4_G[0:2500,1]*nmax)

plt.plot(Res4_KMC[0:2500,0],Res4_KMC[0:2500,1],ls='dashed')
plt.plot(Stall4_KMC[0:2500,0],Stall4_KMC[0:2500,1],ls='dashed')

plt.show()


plt.plot(r5*Res5_G[0:5000,0],Res5_G[0:5000,1]*nmax)
plt.plot(r5*Stall5_G[0:5000,0],Stall5_G[0:5000,1]*nmax)

plt.plot(Res5_KMC[0:5000,0],Res5_KMC[0:5000,1],ls='dashed')
plt.plot(Stall5_KMC[0:5000,0],Stall5_KMC[0:5000,1],ls='dashed')

plt.show()




