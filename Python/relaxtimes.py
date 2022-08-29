#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 16:42:03 2022

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

l=50000

#%%

files030=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=0.00_a=30.44_N=4.dat']
data030=[]
files03=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=0.00_a= 3.06_N=4.dat']
data03=[]
files06=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=0.00_a= 6.15_N=4.dat']
data06=[]

for data_file in files030:
    data030.append(np.loadtxt(data_file))
for data_file in files03:
    data03.append(np.loadtxt(data_file))
for data_file in files06:
    data06.append(np.loadtxt(data_file))

Res030 = data030[0]
Res03 = data03[0]
Res06 = data06[0]


plt.plot(Res030[0:l,0],Res030[0:l,1])
plt.plot(Res03[0:l,0],Res03[0:l,1])
plt.plot(Res06[0:l,0],Res06[0:l,1])

files02530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=0.25_a=30.44_N=4.dat']
data02530=[]
files0253=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=0.25_a= 3.06_N=4.dat']
data0253=[]
files0256=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=0.25_a= 6.15_N=4.dat']
data0256=[]

for data_file in files02530:
    data02530.append(np.loadtxt(data_file))
for data_file in files0253:
    data0253.append(np.loadtxt(data_file))
for data_file in files0256:
    data0256.append(np.loadtxt(data_file))

Res02530 = data02530[0]
Res0253 = data0253[0]
Res0256 = data0256[0]


files05030=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=0.50_a=30.44_N=4.dat']
data05030=[]
files0503=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=0.50_a= 3.06_N=4.dat']
data0503=[]
files0506=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=0.50_a= 6.15_N=4.dat']
data0506=[]

for data_file in files05030:
    data05030.append(np.loadtxt(data_file))
for data_file in files0503:
    data0503.append(np.loadtxt(data_file))
for data_file in files0506:
    data0506.append(np.loadtxt(data_file))

Res05030 = data05030[0]
Res0503 = data0503[0]
Res0506 = data0506[0]


files07530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=0.75_a=30.44_N=4.dat']
data07530=[]
files0753=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=0.75_a= 3.06_N=4.dat']
data0753=[]
files0756=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=0.75_a= 6.15_N=4.dat']
data0756=[]

for data_file in files07530:
    data07530.append(np.loadtxt(data_file))
for data_file in files0753:
    data0753.append(np.loadtxt(data_file))
for data_file in files0756:
    data0756.append(np.loadtxt(data_file))

Res07530 = data07530[0]
Res0753 = data0753[0]
Res0756 = data0756[0]


#%%
files130=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=1.00_a=30.44_N=4.dat']
data130=[]
files13=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=1.00_a= 3.06_N=4.dat']
data13=[]
files16=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=1.00_a= 6.15_N=4.dat']
data16=[]

for data_file in files130:
    data130.append(np.loadtxt(data_file))
for data_file in files13:
    data13.append(np.loadtxt(data_file))
for data_file in files16:
    data16.append(np.loadtxt(data_file))

Res130 = data130[0]
Res13 = data13[0]
Res16 = data16[0]

plt.plot(Res130[0:l,0],Res130[0:l,1])
plt.plot(Res13[0:l,0],Res13[0:l,1])
plt.plot(Res16[0:l,0],Res16[0:l,1])


files12530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=1.25_a=30.44_N=4.dat']
data12530=[]
files1253=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=1.25_a= 3.06_N=4.dat']
data1253=[]
files1256=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=1.25_a= 6.15_N=4.dat']
data1256=[]

for data_file in files12530:
    data12530.append(np.loadtxt(data_file))
for data_file in files1253:
    data1253.append(np.loadtxt(data_file))
for data_file in files1256:
    data1256.append(np.loadtxt(data_file))

Res12530 = data12530[0]
Res1253 = data1253[0]
Res1256 = data1256[0]


files15030=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=1.50_a=30.44_N=4.dat']
data15030=[]
files1503=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=1.50_a= 3.06_N=4.dat']
data1503=[]
files1506=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=1.50_a= 6.15_N=4.dat']
data1506=[]

for data_file in files15030:
    data15030.append(np.loadtxt(data_file))
for data_file in files1503:
    data1503.append(np.loadtxt(data_file))
for data_file in files1506:
    data1506.append(np.loadtxt(data_file))

Res15030 = data15030[0]
Res1503 = data1503[0]
Res1506 = data1506[0]


files17530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=1.75_a=30.44_N=4.dat']
data17530=[]
files1753=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=1.75_a= 3.06_N=4.dat']
data1753=[]
files1756=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=1.75_a= 6.15_N=4.dat']
data1756=[]

for data_file in files17530:
    data17530.append(np.loadtxt(data_file))
for data_file in files1753:
    data1753.append(np.loadtxt(data_file))
for data_file in files1756:
    data1756.append(np.loadtxt(data_file))

Res17530 = data17530[0]
Res1753 = data1753[0]
Res1756 = data1756[0]

#%%

files230=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=2.00_a=30.44_N=4.dat']
data230=[]
files23=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=2.00_a= 3.06_N=4.dat']
data23=[]
files26=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=2.00_a= 6.15_N=4.dat']
data26=[]

for data_file in files230:
    data230.append(np.loadtxt(data_file))
for data_file in files23:
    data23.append(np.loadtxt(data_file))
for data_file in files26:
    data26.append(np.loadtxt(data_file))

Res230 = data230[0]
Res23 = data23[0]
Res26 = data26[0]

plt.plot(Res230[0:l,0],Res230[0:l,1])
plt.plot(Res23[0:l,0],Res23[0:l,1])
plt.plot(Res26[0:l,0],Res26[0:l,1])


files22530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=2.25_a=30.44_N=4.dat']
data22530=[]
files2253=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=2.25_a= 3.06_N=4.dat']
data2253=[]
files2256=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=2.25_a= 6.15_N=4.dat']
data2256=[]

for data_file in files22530:
    data22530.append(np.loadtxt(data_file))
for data_file in files2253:
    data2253.append(np.loadtxt(data_file))
for data_file in files2256:
    data2256.append(np.loadtxt(data_file))

Res22530 = data22530[0]
Res2253 = data2253[0]
Res2256 = data2256[0]


files25030=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=2.50_a=30.44_N=4.dat']
data25030=[]
files2503=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=2.50_a= 3.06_N=4.dat']
data2503=[]
files2506=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=2.50_a= 6.15_N=4.dat']
data2506=[]

for data_file in files25030:
    data25030.append(np.loadtxt(data_file))
for data_file in files2503:
    data2503.append(np.loadtxt(data_file))
for data_file in files2506:
    data2506.append(np.loadtxt(data_file))

Res25030 = data25030[0]
Res2503 = data2503[0]
Res2506 = data2506[0]


files27530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=2.75_a=30.44_N=4.dat']
data27530=[]
files2753=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=2.75_a= 3.06_N=4.dat']
data2753=[]
files2756=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=2.75_a= 6.15_N=4.dat']
data2756=[]

for data_file in files27530:
    data27530.append(np.loadtxt(data_file))
for data_file in files2753:
    data2753.append(np.loadtxt(data_file))
for data_file in files2756:
    data2756.append(np.loadtxt(data_file))

Res27530 = data27530[0]
Res2753 = data2753[0]
Res2756 = data2756[0]


#%%

files330=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=3.00_a=30.44_N=4.dat']
data330=[]
files33=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=3.00_a= 3.06_N=4.dat']
data33=[]
files36=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=3.00_a= 6.15_N=4.dat']
data36=[]

for data_file in files330:
    data330.append(np.loadtxt(data_file))
for data_file in files33:
    data33.append(np.loadtxt(data_file))
for data_file in files36:
    data36.append(np.loadtxt(data_file))

Res330 = data330[0]
Res33 = data33[0]
Res36 = data36[0]

plt.plot(Res330[0:l,0],Res330[0:l,1])
plt.plot(Res33[0:l,0],Res33[0:l,1])
plt.plot(Res36[0:l,0],Res36[0:l,1])



files32530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=3.25_a=30.44_N=4.dat']
data32530=[]
files3253=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=3.25_a= 3.06_N=4.dat']
data3253=[]
files3256=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=3.25_a= 6.15_N=4.dat']
data3256=[]

for data_file in files32530:
    data32530.append(np.loadtxt(data_file))
for data_file in files3253:
    data3253.append(np.loadtxt(data_file))
for data_file in files3256:
    data3256.append(np.loadtxt(data_file))

Res32530 = data32530[0]
Res3253 = data3253[0]
Res3256 = data3256[0]



files35030=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=3.50_a=30.44_N=4.dat']
data35030=[]
files3503=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=3.50_a= 3.06_N=4.dat']
data3503=[]
files3506=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=3.50_a= 6.15_N=4.dat']
data3506=[]

for data_file in files35030:
    data35030.append(np.loadtxt(data_file))
for data_file in files3503:
    data3503.append(np.loadtxt(data_file))
for data_file in files3506:
    data3506.append(np.loadtxt(data_file))

Res35030 = data35030[0]
Res3503 = data3503[0]
Res3506 = data3506[0]



files37530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=3.75_a=30.44_N=4.dat']
data37530=[]
files3753=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=3.75_a= 3.06_N=4.dat']
data3753=[]
files3756=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=3.75_a= 6.15_N=4.dat']
data3756=[]

for data_file in files37530:
    data37530.append(np.loadtxt(data_file))
for data_file in files3753:
    data3753.append(np.loadtxt(data_file))
for data_file in files3756:
    data3756.append(np.loadtxt(data_file))

Res37530 = data37530[0]
Res3753 = data3753[0]
Res3756 = data3756[0]

#%%

files430=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=4.00_a=30.44_N=4.dat']
data430=[]
files43=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=4.00_a= 3.06_N=4.dat']
data43=[]
files46=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=4.00_a= 6.15_N=4.dat']
data46=[]

for data_file in files430:
    data430.append(np.loadtxt(data_file))
for data_file in files43:
    data43.append(np.loadtxt(data_file))
for data_file in files46:
    data46.append(np.loadtxt(data_file))

Res430 = data430[0]
Res43 = data43[0]
Res46 = data46[0]



files42530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=4.25_a=30.44_N=4.dat']
data42530=[]
files4253=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=4.25_a= 3.06_N=4.dat']
data4253=[]
files4256=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=4.25_a= 6.15_N=4.dat']
data4256=[]

for data_file in files42530:
    data42530.append(np.loadtxt(data_file))
for data_file in files4253:
    data4253.append(np.loadtxt(data_file))
for data_file in files4256:
    data4256.append(np.loadtxt(data_file))

Res42530 = data42530[0]
Res4253 = data4253[0]
Res4256 = data4256[0]



files45030=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=4.50_a=30.44_N=4.dat']
data45030=[]
files4503=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=4.50_a= 3.06_N=4.dat']
data4503=[]
files4506=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=4.50_a= 6.15_N=4.dat']
data4506=[]

for data_file in files45030:
    data45030.append(np.loadtxt(data_file))
for data_file in files4503:
    data4503.append(np.loadtxt(data_file))
for data_file in files4506:
    data4506.append(np.loadtxt(data_file))

Res45030 = data45030[0]
Res4503 = data4503[0]
Res4506 = data4506[0]



files47530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=4.75_a=30.44_N=4.dat']
data47530=[]
files4753=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=4.75_a= 3.06_N=4.dat']
data4753=[]
files4756=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=4.75_a= 6.15_N=4.dat']
data4756=[]

for data_file in files47530:
    data47530.append(np.loadtxt(data_file))
for data_file in files4753:
    data4753.append(np.loadtxt(data_file))
for data_file in files4756:
    data4756.append(np.loadtxt(data_file))

Res47530 = data47530[0]
Res4753 = data4753[0]
Res4756 = data4756[0]



files530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=5.00_a=30.44_N=4.dat']
data530=[]
files53=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=5.00_a= 3.06_N=4.dat']
data53=[]
files56=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=5.00_a= 6.15_N=4.dat']
data56=[]

for data_file in files530:
    data530.append(np.loadtxt(data_file))
for data_file in files53:
    data53.append(np.loadtxt(data_file))
for data_file in files56:
    data56.append(np.loadtxt(data_file))

Res530 = data530[0]
Res53 = data53[0]
Res56 = data56[0]

plt.plot(Res530[0:l,0],Res530[0:l,1])
plt.plot(Res53[0:l,0],Res53[0:l,1])
plt.plot(Res56[0:l,0],Res56[0:l,1])


#%%

n0=0.0

tRes = np.zeros((3,21))

poRes030, pcRes030 = curve_fit(N,Res030[0:l,0],Res030[0:l,1])
poRes03, pcRes03 = curve_fit(N,Res03[0:l,0],Res03[0:l,1])
poRes06, pcRes06 = curve_fit(N,Res06[0:l,0],Res06[0:l,1])

poRes02530, pcRes02530 = curve_fit(N,Res02530[0:l,0],Res02530[0:l,1])
poRes0253, pcRes0253 = curve_fit(N,Res0253[0:l,0],Res0253[0:l,1])
poRes0256, pcRes0256 = curve_fit(N,Res0256[0:l,0],Res0256[0:l,1])

poRes05030, pcRes05030 = curve_fit(N,Res05030[0:l,0],Res05030[0:l,1])
poRes0503, pcRes0503 = curve_fit(N,Res0503[0:l,0],Res0503[0:l,1])
poRes0506, pcRes0506 = curve_fit(N,Res0506[0:l,0],Res0506[0:l,1])

poRes07530, pcRes07530 = curve_fit(N,Res07530[0:l,0],Res07530[0:l,1])
poRes0753, pcRes0753 = curve_fit(N,Res0753[0:l,0],Res0753[0:l,1])
poRes0756, pcRes0756 = curve_fit(N,Res0756[0:l,0],Res0756[0:l,1])



poRes130, pcRes130 = curve_fit(N,Res130[0:l,0],Res130[0:l,1])
poRes13, pcRes13 = curve_fit(N,Res13[0:l,0],Res13[0:l,1])
poRes16, pcRes16 = curve_fit(N,Res16[0:l,0],Res16[0:l,1])

poRes17530, pcRes17530 = curve_fit(N,Res17530[0:l,0],Res17530[0:l,1])
poRes1753, pcRes1753 = curve_fit(N,Res1753[0:l,0],Res1753[0:l,1])
poRes1756, pcRes1756 = curve_fit(N,Res1756[0:l,0],Res1756[0:l,1])

poRes12530, pcRes12530 = curve_fit(N,Res12530[0:l,0],Res12530[0:l,1])
poRes1253, pcRes1253 = curve_fit(N,Res1253[0:l,0],Res1253[0:l,1])
poRes1256, pcRes1256 = curve_fit(N,Res1256[0:l,0],Res1256[0:l,1])

poRes15030, pcRes15030 = curve_fit(N,Res15030[0:l,0],Res15030[0:l,1])
poRes1503, pcRes1503 = curve_fit(N,Res1503[0:l,0],Res1503[0:l,1])
poRes1506, pcRes1506 = curve_fit(N,Res1506[0:l,0],Res1506[0:l,1])




poRes230, pcRes230 = curve_fit(N,Res230[0:l,0],Res230[0:l,1])
poRes23, pcRes23 = curve_fit(N,Res23[0:l,0],Res23[0:l,1])
poRes26, pcRes26 = curve_fit(N,Res26[0:l,0],Res26[0:l,1])

poRes22530, pcRes230 = curve_fit(N,Res22530[0:l,0],Res22530[0:l,1])
poRes2253, pcRes2253 = curve_fit(N,Res2253[0:l,0],Res2253[0:l,1])
poRes2256, pcRes2256 = curve_fit(N,Res2256[0:l,0],Res2256[0:l,1])

poRes25030, pcRes25030 = curve_fit(N,Res25030[0:l,0],Res25030[0:l,1])
poRes2503, pcRes2503 = curve_fit(N,Res2503[0:l,0],Res2503[0:l,1])
poRes2506, pcRes2506 = curve_fit(N,Res2506[0:l,0],Res2506[0:l,1])

poRes27530, pcRes27530 = curve_fit(N,Res27530[0:l,0],Res27530[0:l,1])
poRes2753, pcRes2753 = curve_fit(N,Res2753[0:l,0],Res2753[0:l,1])
poRes2756, pcRes2756 = curve_fit(N,Res2756[0:l,0],Res2756[0:l,1])



poRes330, pcRes330 = curve_fit(N,Res330[0:l,0],Res330[0:l,1])
poRes33, pcRes33 = curve_fit(N,Res33[0:l,0],Res33[0:l,1])
poRes36, pcRes36 = curve_fit(N,Res36[0:l,0],Res36[0:l,1])

poRes32530, pcRes32530 = curve_fit(N,Res32530[0:l,0],Res32530[0:l,1])
poRes3253, pcRes3253 = curve_fit(N,Res3253[0:l,0],Res3253[0:l,1])
poRes3256, pcRes3256 = curve_fit(N,Res3256[0:l,0],Res3256[0:l,1])

poRes35030, pcRes35030 = curve_fit(N,Res35030[0:l,0],Res35030[0:l,1])
poRes3503, pcRes3503 = curve_fit(N,Res3503[0:l,0],Res3503[0:l,1])
poRes3506, pcRes3506 = curve_fit(N,Res3506[0:l,0],Res3506[0:l,1])

poRes37530, pcRes37530 = curve_fit(N,Res37530[0:l,0],Res37530[0:l,1])
poRes3753, pcRes3753 = curve_fit(N,Res3753[0:l,0],Res3753[0:l,1])
poRes3756, pcRes3756 = curve_fit(N,Res3756[0:l,0],Res3756[0:l,1])



poRes430, pcRes430 = curve_fit(N,Res430[0:l,0],Res430[0:l,1])
poRes43, pcRes43 = curve_fit(N,Res43[0:l,0],Res43[0:l,1])
poRes46, pcRes46 = curve_fit(N,Res46[0:l,0],Res46[0:l,1])

poRes42530, pcRes42530 = curve_fit(N,Res42530[0:l,0],Res42530[0:l,1])
poRes4253, pcRes4253 = curve_fit(N,Res4253[0:l,0],Res4253[0:l,1])
poRes4256, pcRes4256 = curve_fit(N,Res4256[0:l,0],Res4256[0:l,1])

poRes45030, pcRes45030 = curve_fit(N,Res45030[0:l,0],Res45030[0:l,1])
poRes4503, pcRes4503 = curve_fit(N,Res4503[0:l,0],Res4503[0:l,1])
poRes4506, pcRes4506 = curve_fit(N,Res4506[0:l,0],Res4506[0:l,1])

poRes47530, pcRes47530 = curve_fit(N,Res47530[0:l,0],Res47530[0:l,1])
poRes4753, pcRes4753 = curve_fit(N,Res4753[0:l,0],Res4753[0:l,1])
poRes4756, pcRes4756 = curve_fit(N,Res4756[0:l,0],Res4756[0:l,1])


poRes530, pcRes530 = curve_fit(N,Res530[0:l,0],Res530[0:l,1])
poRes53, pcRes53 = curve_fit(N,Res53[0:l,0],Res53[0:l,1])
poRes56, pcRes56 = curve_fit(N,Res56[0:l,0],Res56[0:l,1])

#%%

tRes[0,0] = poRes03[0]
tRes[1,0] = poRes06[0]
tRes[2,0] = poRes030[0]

tRes[0,1] = poRes0253[0]
tRes[1,1] = poRes0256[0]
tRes[2,1] = poRes02530[0]

tRes[0,2] = poRes0503[0]
tRes[1,2] = poRes0506[0]
tRes[2,2] = poRes05030[0]

tRes[0,3] = poRes0753[0]
tRes[1,3] = poRes0756[0]
tRes[2,3] = poRes07530[0]

tRes[0,4] = poRes13[0]
tRes[1,4] = poRes16[0]
tRes[2,4] = poRes130[0]

tRes[0,5] = poRes1253[0]
tRes[1,5] = poRes1256[0]
tRes[2,5] = poRes12530[0]

tRes[0,6] = poRes1503[0]
tRes[1,6] = poRes1506[0]
tRes[2,6] = poRes15030[0]

tRes[0,7] = poRes1753[0]
tRes[1,7] = poRes1756[0]
tRes[2,7] = poRes17530[0]

tRes[0,8] = poRes23[0]
tRes[1,8] = poRes26[0]
tRes[2,8] = poRes230[0]

tRes[0,9] = poRes2253[0]
tRes[1,9] = poRes2256[0]
tRes[2,9] = poRes22530[0]

tRes[0,10] = poRes2503[0]
tRes[1,10] = poRes2506[0]
tRes[2,10] = poRes25030[0]

tRes[0,11] = poRes2753[0]
tRes[1,11] = poRes2756[0]
tRes[2,11] = poRes27530[0]

tRes[0,12] = poRes33[0]
tRes[1,12] = poRes36[0]
tRes[2,12] = poRes330[0]

tRes[0,13] = poRes3253[0]
tRes[1,13] = poRes3256[0]
tRes[2,13] = poRes32530[0]

tRes[0,14] = poRes3503[0]
tRes[1,14] = poRes3506[0]
tRes[2,14] = poRes35030[0]

tRes[0,15] = poRes3753[0]
tRes[1,15] = poRes3756[0]
tRes[2,15] = poRes37530[0]

tRes[0,16] = poRes43[0]
tRes[1,16] = poRes46[0]
tRes[2,16] = poRes430[0]

tRes[0,17] = poRes4253[0]
tRes[1,17] = poRes4256[0]
tRes[2,17] = poRes42530[0]

tRes[0,18] = poRes4503[0]
tRes[1,18] = poRes4506[0]
tRes[2,18] = poRes45030[0]

tRes[0,19] = poRes4753[0]
tRes[1,19] = poRes4756[0]
tRes[2,19] = poRes47530[0]

tRes[0,20] = poRes53[0]
tRes[1,20] = poRes56[0]
tRes[2,20] = poRes530[0]


tRes

#%%

m=21

J=np.arange(0,5.25,0.25)


plt.rcParams["figure.figsize"] =[8.0,6.0]
plt.title('Resurrection relax times',fontsize=16)
plt.xlabel('J  ($k_BT$)', fontsize=16)
plt.ylabel(r'$\tau_{res}$ (s)',fontsize=16)
plt.grid()


plt.plot(J,tRes[0,:],marker='.',linestyle="",markersize=7,label=r'$\alpha=3.06$')
plt.plot(J,tRes[1,:],marker='.',linestyle="",markersize=7,label=r'$\alpha=6.15$')
plt.plot(J,tRes[2,:],marker='.',linestyle="",markersize=7,label=r'$\alpha=30.44$')

plt.legend(loc='upper left',prop={'size': 16})

#%%

filesS030=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=0.00_a=30.44_N=4.dat']
StallD030=[]
filesS03=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=0.00_a= 3.06_N=4.dat']
StallD03=[]
filesS06=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=0.00_a= 6.15_N=4.dat']
StallD06=[]

for data_file in filesS030:
    StallD030.append(np.loadtxt(data_file))
for data_file in filesS03:
    StallD03.append(np.loadtxt(data_file))
for data_file in filesS06:
    StallD06.append(np.loadtxt(data_file))

Stall030 = StallD030[0]
Stall03 = StallD03[0]
Stall06 = StallD06[0]


plt.plot(Stall030[0:l,0],Stall030[0:l,1])
plt.plot(Stall03[0:l,0],Stall03[0:l,1])
plt.plot(Stall06[0:l,0],Stall06[0:l,1])

filesS02530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=0.25_a=30.44_N=4.dat']
StallD02530=[]
filesS0253=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=0.25_a= 3.06_N=4.dat']
StallD0253=[]
filesS0256=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=0.25_a= 6.15_N=4.dat']
StallD0256=[]

for data_file in filesS02530:
    StallD02530.append(np.loadtxt(data_file))
for data_file in filesS0253:
    StallD0253.append(np.loadtxt(data_file))
for data_file in filesS0256:
    StallD0256.append(np.loadtxt(data_file))

Stall02530 = StallD02530[0]
Stall0253 = StallD0253[0]
Stall0256 = StallD0256[0]


filesS05030=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=0.50_a=30.44_N=4.dat']
StallD05030=[]
filesS0503=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=0.50_a= 3.06_N=4.dat']
StallD0503=[]
filesS0506=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=0.50_a= 6.15_N=4.dat']
StallD0506=[]

for data_file in filesS05030:
    StallD05030.append(np.loadtxt(data_file))
for data_file in filesS0503:
    StallD0503.append(np.loadtxt(data_file))
for data_file in filesS0506:
    StallD0506.append(np.loadtxt(data_file))

Stall05030 = StallD05030[0]
Stall0503 = StallD0503[0]
Stall0506 = StallD0506[0]


filesS07530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=0.75_a=30.44_N=4.dat']
StallD07530=[]
filesS0753=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=0.75_a= 3.06_N=4.dat']
StallD0753=[]
filesS0756=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=0.75_a= 6.15_N=4.dat']
StallD0756=[]

for data_file in filesS07530:
    StallD07530.append(np.loadtxt(data_file))
for data_file in filesS0753:
    StallD0753.append(np.loadtxt(data_file))
for data_file in filesS0756:
    StallD0756.append(np.loadtxt(data_file))

Stall07530 = StallD07530[0]
Stall0753 = StallD0753[0]
Stall0756 = StallD0756[0]


#%%
filesS130=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=1.00_a=30.44_N=4.dat']
StallD130=[]
filesS13=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=1.00_a= 3.06_N=4.dat']
StallD13=[]
filesS16=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=1.00_a= 6.15_N=4.dat']
StallD16=[]

for data_file in filesS130:
    StallD130.append(np.loadtxt(data_file))
for data_file in filesS13:
    StallD13.append(np.loadtxt(data_file))
for data_file in filesS16:
    StallD16.append(np.loadtxt(data_file))

Stall130 = StallD130[0]
Stall13 = StallD13[0]
Stall16 = StallD16[0]

plt.plot(Stall130[0:l,0],Stall130[0:l,1])
plt.plot(Stall13[0:l,0],Stall13[0:l,1])
plt.plot(Stall16[0:l,0],Stall16[0:l,1])


filesS12530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=1.25_a=30.44_N=4.dat']
StallD12530=[]
filesS1253=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=1.25_a= 3.06_N=4.dat']
StallD1253=[]
filesS1256=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=1.25_a= 6.15_N=4.dat']
StallD1256=[]

for data_file in filesS12530:
    StallD12530.append(np.loadtxt(data_file))
for data_file in filesS1253:
    StallD1253.append(np.loadtxt(data_file))
for data_file in filesS1256:
    StallD1256.append(np.loadtxt(data_file))

Stall12530 = StallD12530[0]
Stall1253 = StallD1253[0]
Stall1256 = StallD1256[0]


filesS15030=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=1.50_a=30.44_N=4.dat']
StallD15030=[]
filesS1503=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=1.50_a= 3.06_N=4.dat']
StallD1503=[]
filesS1506=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=1.50_a= 6.15_N=4.dat']
StallD1506=[]

for data_file in filesS15030:
    StallD15030.append(np.loadtxt(data_file))
for data_file in filesS1503:
    StallD1503.append(np.loadtxt(data_file))
for data_file in filesS1506:
    StallD1506.append(np.loadtxt(data_file))

Stall15030 = StallD15030[0]
Stall1503 = StallD1503[0]
Stall1506 = StallD1506[0]


filesS17530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=1.75_a=30.44_N=4.dat']
StallD17530=[]
filesS1753=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=1.75_a= 3.06_N=4.dat']
StallD1753=[]
filesS1756=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=1.75_a= 6.15_N=4.dat']
StallD1756=[]

for data_file in filesS17530:
    StallD17530.append(np.loadtxt(data_file))
for data_file in filesS1753:
    StallD1753.append(np.loadtxt(data_file))
for data_file in filesS1756:
    StallD1756.append(np.loadtxt(data_file))

Stall17530 = StallD17530[0]
Stall1753 = StallD1753[0]
Stall1756 = StallD1756[0]

#%%

filesS230=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=2.00_a=30.44_N=4.dat']
StallD230=[]
filesS23=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=2.00_a= 3.06_N=4.dat']
StallD23=[]
filesS26=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=2.00_a= 6.15_N=4.dat']
StallD26=[]

for data_file in filesS230:
    StallD230.append(np.loadtxt(data_file))
for data_file in filesS23:
    StallD23.append(np.loadtxt(data_file))
for data_file in filesS26:
    StallD26.append(np.loadtxt(data_file))

Stall230 = StallD230[0]
Stall23 = StallD23[0]
Stall26 = StallD26[0]

plt.plot(Stall230[0:l,0],Stall230[0:l,1])
plt.plot(Stall23[0:l,0],Stall23[0:l,1])
plt.plot(Stall26[0:l,0],Stall26[0:l,1])


filesS22530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=2.25_a=30.44_N=4.dat']
StallD22530=[]
filesS2253=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=2.25_a= 3.06_N=4.dat']
StallD2253=[]
filesS2256=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=2.25_a= 6.15_N=4.dat']
StallD2256=[]

for data_file in filesS22530:
    StallD22530.append(np.loadtxt(data_file))
for data_file in filesS2253:
    StallD2253.append(np.loadtxt(data_file))
for data_file in filesS2256:
    StallD2256.append(np.loadtxt(data_file))

Stall22530 = StallD22530[0]
Stall2253 = StallD2253[0]
Stall2256 = StallD2256[0]


filesS25030=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=2.50_a=30.44_N=4.dat']
StallD25030=[]
filesS2503=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=2.50_a= 3.06_N=4.dat']
StallD2503=[]
filesS2506=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=2.50_a= 6.15_N=4.dat']
StallD2506=[]

for data_file in filesS25030:
    StallD25030.append(np.loadtxt(data_file))
for data_file in filesS2503:
    StallD2503.append(np.loadtxt(data_file))
for data_file in filesS2506:
    StallD2506.append(np.loadtxt(data_file))

Stall25030 = StallD25030[0]
Stall2503 = StallD2503[0]
Stall2506 = StallD2506[0]


filesS27530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=2.75_a=30.44_N=4.dat']
StallD27530=[]
filesS2753=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=2.75_a= 3.06_N=4.dat']
StallD2753=[]
filesS2756=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=2.75_a= 6.15_N=4.dat']
StallD2756=[]

for data_file in filesS27530:
    StallD27530.append(np.loadtxt(data_file))
for data_file in filesS2753:
    StallD2753.append(np.loadtxt(data_file))
for data_file in filesS2756:
    StallD2756.append(np.loadtxt(data_file))

Stall27530 = StallD27530[0]
Stall2753 = StallD2753[0]
Stall2756 = StallD2756[0]


#%%

filesS330=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=3.00_a=30.44_N=4.dat']
StallD330=[]
filesS33=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=3.00_a= 3.06_N=4.dat']
StallD33=[]
filesS36=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=3.00_a= 6.15_N=4.dat']
StallD36=[]

for data_file in filesS330:
    StallD330.append(np.loadtxt(data_file))
for data_file in filesS33:
    StallD33.append(np.loadtxt(data_file))
for data_file in filesS36:
    StallD36.append(np.loadtxt(data_file))

Stall330 = StallD330[0]
Stall33 = StallD33[0]
Stall36 = StallD36[0]

plt.plot(Stall330[0:l,0],Stall330[0:l,1])
plt.plot(Stall33[0:l,0],Stall33[0:l,1])
plt.plot(Stall36[0:l,0],Stall36[0:l,1])



filesS32530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=3.25_a=30.44_N=4.dat']
StallD32530=[]
filesS3253=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=3.25_a= 3.06_N=4.dat']
StallD3253=[]
filesS3256=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=3.25_a= 6.15_N=4.dat']
StallD3256=[]

for data_file in filesS32530:
    StallD32530.append(np.loadtxt(data_file))
for data_file in filesS3253:
    StallD3253.append(np.loadtxt(data_file))
for data_file in filesS3256:
    StallD3256.append(np.loadtxt(data_file))

Stall32530 = StallD32530[0]
Stall3253 = StallD3253[0]
Stall3256 = StallD3256[0]



filesS35030=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=3.50_a=30.44_N=4.dat']
StallD35030=[]
filesS3503=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=3.50_a= 3.06_N=4.dat']
StallD3503=[]
filesS3506=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=3.50_a= 6.15_N=4.dat']
StallD3506=[]

for data_file in filesS35030:
    StallD35030.append(np.loadtxt(data_file))
for data_file in filesS3503:
    StallD3503.append(np.loadtxt(data_file))
for data_file in filesS3506:
    StallD3506.append(np.loadtxt(data_file))

Stall35030 = StallD35030[0]
Stall3503 = StallD3503[0]
Stall3506 = StallD3506[0]



filesS37530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=3.75_a=30.44_N=4.dat']
StallD37530=[]
filesS3753=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=3.75_a= 3.06_N=4.dat']
StallD3753=[]
filesS3756=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=3.75_a= 6.15_N=4.dat']
StallD3756=[]

for data_file in filesS37530:
    StallD37530.append(np.loadtxt(data_file))
for data_file in filesS3753:
    StallD3753.append(np.loadtxt(data_file))
for data_file in filesS3756:
    StallD3756.append(np.loadtxt(data_file))

Stall37530 = StallD37530[0]
Stall3753 = StallD3753[0]
Stall3756 = StallD3756[0]

#%%

filesS430=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=4.00_a=30.44_N=4.dat']
StallD430=[]
filesS43=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=4.00_a= 3.06_N=4.dat']
StallD43=[]
filesS46=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=4.00_a= 6.15_N=4.dat']
StallD46=[]

for data_file in filesS430:
    StallD430.append(np.loadtxt(data_file))
for data_file in filesS43:
    StallD43.append(np.loadtxt(data_file))
for data_file in filesS46:
    StallD46.append(np.loadtxt(data_file))

Stall430 = StallD430[0]
Stall43 = StallD43[0]
Stall46 = StallD46[0]




filesS42530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=4.25_a=30.44_N=4.dat']
StallD42530=[]
filesS4253=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=4.25_a= 3.06_N=4.dat']
StallD4253=[]
filesS4256=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=4.25_a= 6.15_N=4.dat']
StallD4256=[]

for data_file in filesS42530:
    StallD42530.append(np.loadtxt(data_file))
for data_file in filesS4253:
    StallD4253.append(np.loadtxt(data_file))
for data_file in filesS4256:
    StallD4256.append(np.loadtxt(data_file))

Stall42530 = StallD42530[0]
Stall4253 = StallD4253[0]
Stall4256 = StallD4256[0]



filesS45030=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=4.50_a=30.44_N=4.dat']
StallD45030=[]
filesS4503=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=4.50_a= 3.06_N=4.dat']
StallD4503=[]
filesS4506=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=4.50_a= 6.15_N=4.dat']
StallD4506=[]

for data_file in filesS45030:
    StallD45030.append(np.loadtxt(data_file))
for data_file in filesS4503:
    StallD4503.append(np.loadtxt(data_file))
for data_file in filesS4506:
    StallD4506.append(np.loadtxt(data_file))

Stall45030 = StallD45030[0]
Stall4503 = StallD4503[0]
Stall4506 = StallD4506[0]



filesS47530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=4.75_a=30.44_N=4.dat']
StallD47530=[]
filesS4753=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=4.75_a= 3.06_N=4.dat']
StallD4753=[]
filesS4756=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=4.75_a= 6.15_N=4.dat']
StallD4756=[]

for data_file in filesS47530:
    StallD47530.append(np.loadtxt(data_file))
for data_file in filesS4753:
    StallD4753.append(np.loadtxt(data_file))
for data_file in filesS4756:
    StallD4756.append(np.loadtxt(data_file))

Stall47530 = StallD47530[0]
Stall4753 = StallD4753[0]
Stall4756 = StallD4756[0]



filesS530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=5.00_a=30.44_N=4.dat']
StallD530=[]
filesS53=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=5.00_a= 3.06_N=4.dat']
StallD53=[]
filesS56=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=5.00_a= 6.15_N=4.dat']
StallD56=[]

for data_file in filesS530:
    StallD530.append(np.loadtxt(data_file))
for data_file in filesS53:
    StallD53.append(np.loadtxt(data_file))
for data_file in filesS56:
    StallD56.append(np.loadtxt(data_file))

Stall530 = StallD530[0]
Stall53 = StallD53[0]
Stall56 = StallD56[0]

plt.plot(Stall530[0:l,0],Stall530[0:l,1])
plt.plot(Stall53[0:l,0],Stall53[0:l,1])
plt.plot(Stall56[0:l,0],Stall56[0:l,1])

#%%

n0=6.0

tStall = np.zeros((3,21))

poStall030, pcStall030 = curve_fit(N,Stall030[0:l,0],Stall030[0:l,1])
poStall03, pcStall03 = curve_fit(N,Stall03[0:l,0],Stall03[0:l,1])
poStall06, pcStall06 = curve_fit(N,Stall06[0:l,0],Stall06[0:l,1])

poStall02530, pcStall02530 = curve_fit(N,Stall02530[0:l,0],Stall02530[0:l,1])
poStall0253, pcStall0253 = curve_fit(N,Stall0253[0:l,0],Stall0253[0:l,1])
poStall0256, pcStall0256 = curve_fit(N,Stall0256[0:l,0],Stall0256[0:l,1])

poStall05030, pcStall05030 = curve_fit(N,Stall05030[0:l,0],Stall05030[0:l,1])
poStall0503, pcStall0503 = curve_fit(N,Stall0503[0:l,0],Stall0503[0:l,1])
poStall0506, pcStall0506 = curve_fit(N,Stall0506[0:l,0],Stall0506[0:l,1])

poStall07530, pcStall07530 = curve_fit(N,Stall07530[0:l,0],Stall07530[0:l,1])
poStall0753, pcStall0753 = curve_fit(N,Stall0753[0:l,0],Stall0753[0:l,1])
poStall0756, pcStall0756 = curve_fit(N,Stall0756[0:l,0],Stall0756[0:l,1])



poStall130, pcStall130 = curve_fit(N,Stall130[0:l,0],Stall130[0:l,1])
poStall13, pcStall13 = curve_fit(N,Stall13[0:l,0],Stall13[0:l,1])
poStall16, pcStall16 = curve_fit(N,Stall16[0:l,0],Stall16[0:l,1])

poStall17530, pcStall17530 = curve_fit(N,Stall17530[0:l,0],Stall17530[0:l,1])
poStall1753, pcStall1753 = curve_fit(N,Stall1753[0:l,0],Stall1753[0:l,1])
poStall1756, pcStall1756 = curve_fit(N,Stall1756[0:l,0],Stall1756[0:l,1])

poStall12530, pcStall12530 = curve_fit(N,Stall12530[0:l,0],Stall12530[0:l,1])
poStall1253, pcStall1253 = curve_fit(N,Stall1253[0:l,0],Stall1253[0:l,1])
poStall1256, pcStall1256 = curve_fit(N,Stall1256[0:l,0],Stall1256[0:l,1])

poStall15030, pcStall15030 = curve_fit(N,Stall15030[0:l,0],Stall15030[0:l,1])
poStall1503, pcStall1503 = curve_fit(N,Stall1503[0:l,0],Stall1503[0:l,1])
poStall1506, pcStall1506 = curve_fit(N,Stall1506[0:l,0],Stall1506[0:l,1])




poStall230, pcStall230 = curve_fit(N,Stall230[0:l,0],Stall230[0:l,1])
poStall23, pcStall23 = curve_fit(N,Stall23[0:l,0],Stall23[0:l,1])
poStall26, pcStall26 = curve_fit(N,Stall26[0:l,0],Stall26[0:l,1])

poStall22530, pcStall230 = curve_fit(N,Stall22530[0:l,0],Stall22530[0:l,1])
poStall2253, pcStall2253 = curve_fit(N,Stall2253[0:l,0],Stall2253[0:l,1])
poStall2256, pcStall2256 = curve_fit(N,Stall2256[0:l,0],Stall2256[0:l,1])

poStall25030, pcStall25030 = curve_fit(N,Stall25030[0:l,0],Stall25030[0:l,1])
poStall2503, pcStall2503 = curve_fit(N,Stall2503[0:l,0],Stall2503[0:l,1])
poStall2506, pcStall2506 = curve_fit(N,Stall2506[0:l,0],Stall2506[0:l,1])

poStall27530, pcStall27530 = curve_fit(N,Stall27530[0:l,0],Stall27530[0:l,1])
poStall2753, pcStall2753 = curve_fit(N,Stall2753[0:l,0],Stall2753[0:l,1])
poStall2756, pcStall2756 = curve_fit(N,Stall2756[0:l,0],Stall2756[0:l,1])



poStall330, pcStall330 = curve_fit(N,Stall330[0:l,0],Stall330[0:l,1])
poStall33, pcStall33 = curve_fit(N,Stall33[0:l,0],Stall33[0:l,1])
poStall36, pcStall36 = curve_fit(N,Stall36[0:l,0],Stall36[0:l,1])

poStall32530, pcStall32530 = curve_fit(N,Stall32530[0:l,0],Stall32530[0:l,1])
poStall3253, pcStall3253 = curve_fit(N,Stall3253[0:l,0],Stall3253[0:l,1])
poStall3256, pcStall3256 = curve_fit(N,Stall3256[0:l,0],Stall3256[0:l,1])

poStall35030, pcStall35030 = curve_fit(N,Stall35030[0:l,0],Stall35030[0:l,1])
poStall3503, pcStall3503 = curve_fit(N,Stall3503[0:l,0],Stall3503[0:l,1])
poStall3506, pcStall3506 = curve_fit(N,Stall3506[0:l,0],Stall3506[0:l,1])

poStall37530, pcStall37530 = curve_fit(N,Stall37530[0:l,0],Stall37530[0:l,1])
poStall3753, pcStall3753 = curve_fit(N,Stall3753[0:l,0],Stall3753[0:l,1])
poStall3756, pcStall3756 = curve_fit(N,Stall3756[0:l,0],Stall3756[0:l,1])



poStall430, pcStall430 = curve_fit(N,Stall430[0:l,0],Stall430[0:l,1])
poStall43, pcStall43 = curve_fit(N,Stall43[0:l,0],Stall43[0:l,1])
poStall46, pcStall46 = curve_fit(N,Stall46[0:l,0],Stall46[0:l,1])

poStall42530, pcStall42530 = curve_fit(N,Stall42530[0:l,0],Stall42530[0:l,1])
poStall4253, pcStall4253 = curve_fit(N,Stall4253[0:l,0],Stall4253[0:l,1])
poStall4256, pcStall4256 = curve_fit(N,Stall4256[0:l,0],Stall4256[0:l,1])

poStall45030, pcStall45030 = curve_fit(N,Stall45030[0:l,0],Stall45030[0:l,1])
poStall4503, pcStall4503 = curve_fit(N,Stall4503[0:l,0],Stall4503[0:l,1])
poStall4506, pcStall4506 = curve_fit(N,Stall4506[0:l,0],Stall4506[0:l,1])

poStall47530, pcStall47530 = curve_fit(N,Stall47530[0:l,0],Stall47530[0:l,1])
poStall4753, pcStall4753 = curve_fit(N,Stall4753[0:l,0],Stall4753[0:l,1])
poStall4756, pcStall4756 = curve_fit(N,Stall4756[0:l,0],Stall4756[0:l,1])


poStall530, pcStall530 = curve_fit(N,Stall530[0:l,0],Stall530[0:l,1])
poStall53, pcStall53 = curve_fit(N,Stall53[0:l,0],Stall53[0:l,1])
poStall56, pcStall56 = curve_fit(N,Stall56[0:l,0],Stall56[0:l,1])

#%%

tStall[0,0] = poStall03[0]
tStall[1,0] = poStall06[0]
tStall[2,0] = poStall030[0]

tStall[0,1] = poStall0253[0]
tStall[1,1] = poStall0256[0]
tStall[2,1] = poStall02530[0]

tStall[0,2] = poStall0503[0]
tStall[1,2] = poStall0506[0]
tStall[2,2] = poStall05030[0]

tStall[0,3] = poStall0753[0]
tStall[1,3] = poStall0756[0]
tStall[2,3] = poStall07530[0]

tStall[0,4] = poStall13[0]
tStall[1,4] = poStall16[0]
tStall[2,4] = poStall130[0]

tStall[0,5] = poStall1253[0]
tStall[1,5] = poStall1256[0]
tStall[2,5] = poStall12530[0]

tStall[0,6] = poStall1503[0]
tStall[1,6] = poStall1506[0]
tStall[2,6] = poStall15030[0]

tStall[0,7] = poStall1753[0]
tStall[1,7] = poStall1756[0]
tStall[2,7] = poStall17530[0]

tStall[0,8] = poStall23[0]
tStall[1,8] = poStall26[0]
tStall[2,8] = poStall230[0]

tStall[0,9] = poStall2253[0]
tStall[1,9] = poStall2256[0]
tStall[2,9] = poStall22530[0]

tStall[0,10] = poStall2503[0]
tStall[1,10] = poStall2506[0]
tStall[2,10] = poStall25030[0]

tStall[0,11] = poStall2753[0]
tStall[1,11] = poStall2756[0]
tStall[2,11] = poStall27530[0]

tStall[0,12] = poStall33[0]
tStall[1,12] = poStall36[0]
tStall[2,12] = poStall330[0]

tStall[0,13] = poStall3253[0]
tStall[1,13] = poStall3256[0]
tStall[2,13] = poStall32530[0]

tStall[0,14] = poStall3503[0]
tStall[1,14] = poStall3506[0]
tStall[2,14] = poStall35030[0]

tStall[0,15] = poStall3753[0]
tStall[1,15] = poStall3756[0]
tStall[2,15] = poStall37530[0]

tStall[0,16] = poStall43[0]
tStall[1,16] = poStall46[0]
tStall[2,16] = poStall430[0]

tStall[0,17] = poStall4253[0]
tStall[1,17] = poStall4256[0]
tStall[2,17] = poStall42530[0]

tStall[0,18] = poStall4503[0]
tStall[1,18] = poStall4506[0]
tStall[2,18] = poStall45030[0]

tStall[0,19] = poStall4753[0]
tStall[1,19] = poStall4756[0]
tStall[2,19] = poStall47530[0]

tStall[0,20] = poStall53[0]
tStall[1,20] = poStall56[0]
tStall[2,20] = poStall530[0]


tStall

#%%

m=21

J=np.arange(0,5.25,0.25)


plt.rcParams["figure.figsize"] =[8.0,6.0]
plt.title('Resurrection relax times',fontsize=16)
plt.xlabel('J  ($k_BT$)', fontsize=16)
plt.ylabel(r'$\tau_{rel}$ (s)',fontsize=16)
plt.grid()


plt.plot(J,tStall[0,:],marker='.',linestyle="",markersize=7,label=r'$\alpha=3.06$')
plt.plot(J,tStall[1,:],marker='.',linestyle="",markersize=7,label=r'$\alpha=6.15$')
plt.plot(J,tStall[2,:],marker='.',linestyle="",markersize=7,label=r'$\alpha=30.44$')

plt.legend(loc='upper left',prop={'size': 16})

#%%

Dt03_N4= (tStall[0,:]-tRes[0,:])/(tStall[0,:]+tRes[0,:])
Dt06_N4= (tStall[1,:]-tRes[1,:])/(tStall[1,:]+tRes[1,:])
Dt30_N4= (tStall[2,:]-tRes[2,:])/(tStall[2,:]+tRes[2,:])

plt.rcParams["figure.figsize"] =[8.0,6.0]
plt.title('',fontsize=16)
plt.xlabel('J  ($k_BT$)', fontsize=16)
plt.ylabel(r'$\Delta$',fontsize=16)
plt.grid()

plt.plot(J,Dt03_N4,marker='.',linestyle="",markersize=7,label=r'$\alpha=3.06$')
plt.plot(J,Dt06_N4,marker='.',linestyle="",markersize=7,label=r'$\alpha=6.15$')
plt.plot(J,Dt30_N4,marker='.',linestyle="",markersize=7,label=r'$\alpha=30.44$')

plt.legend(loc='upper right',prop={'size': 16}, framealpha=0.3)




#%%

###############################################################################

########################## STEADY STATE = 8 ###################################

###############################################################################

files030=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=0.00_a=30.44_N=8.dat']
data030=[]
files03=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=0.00_a= 3.06_N=8.dat']
data03=[]
files06=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=0.00_a= 6.15_N=8.dat']
data06=[]

for data_file in files030:
    data030.append(np.loadtxt(data_file))
for data_file in files03:
    data03.append(np.loadtxt(data_file))
for data_file in files06:
    data06.append(np.loadtxt(data_file))

Res030 = data030[0]
Res03 = data03[0]
Res06 = data06[0]


plt.plot(Res030[0:l,0],Res030[0:l,1])
plt.plot(Res03[0:l,0],Res03[0:l,1])
plt.plot(Res06[0:l,0],Res06[0:l,1])

files02530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=0.25_a=30.44_N=8.dat']
data02530=[]
files0253=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=0.25_a= 3.06_N=8.dat']
data0253=[]
files0256=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=0.25_a= 6.15_N=8.dat']
data0256=[]

for data_file in files02530:
    data02530.append(np.loadtxt(data_file))
for data_file in files0253:
    data0253.append(np.loadtxt(data_file))
for data_file in files0256:
    data0256.append(np.loadtxt(data_file))

Res02530 = data02530[0]
Res0253 = data0253[0]
Res0256 = data0256[0]


files05030=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=0.50_a=30.44_N=8.dat']
data05030=[]
files0503=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=0.50_a= 3.06_N=8.dat']
data0503=[]
files0506=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=0.50_a= 6.15_N=8.dat']
data0506=[]

for data_file in files05030:
    data05030.append(np.loadtxt(data_file))
for data_file in files0503:
    data0503.append(np.loadtxt(data_file))
for data_file in files0506:
    data0506.append(np.loadtxt(data_file))

Res05030 = data05030[0]
Res0503 = data0503[0]
Res0506 = data0506[0]


files07530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=0.75_a=30.44_N=8.dat']
data07530=[]
files0753=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=0.75_a= 3.06_N=8.dat']
data0753=[]
files0756=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=0.75_a= 6.15_N=8.dat']
data0756=[]

for data_file in files07530:
    data07530.append(np.loadtxt(data_file))
for data_file in files0753:
    data0753.append(np.loadtxt(data_file))
for data_file in files0756:
    data0756.append(np.loadtxt(data_file))

Res07530 = data07530[0]
Res0753 = data0753[0]
Res0756 = data0756[0]


#%%
files130=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=1.00_a=30.44_N=8.dat']
data130=[]
files13=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=1.00_a= 3.06_N=8.dat']
data13=[]
files16=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=1.00_a= 6.15_N=8.dat']
data16=[]

for data_file in files130:
    data130.append(np.loadtxt(data_file))
for data_file in files13:
    data13.append(np.loadtxt(data_file))
for data_file in files16:
    data16.append(np.loadtxt(data_file))

Res130 = data130[0]
Res13 = data13[0]
Res16 = data16[0]

plt.plot(Res130[0:l,0],Res130[0:l,1])
plt.plot(Res13[0:l,0],Res13[0:l,1])
plt.plot(Res16[0:l,0],Res16[0:l,1])


files12530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=1.25_a=30.44_N=8.dat']
data12530=[]
files1253=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=1.25_a= 3.06_N=8.dat']
data1253=[]
files1256=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=1.25_a= 6.15_N=8.dat']
data1256=[]

for data_file in files12530:
    data12530.append(np.loadtxt(data_file))
for data_file in files1253:
    data1253.append(np.loadtxt(data_file))
for data_file in files1256:
    data1256.append(np.loadtxt(data_file))

Res12530 = data12530[0]
Res1253 = data1253[0]
Res1256 = data1256[0]


files15030=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=1.50_a=30.44_N=8.dat']
data15030=[]
files1503=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=1.50_a= 3.06_N=8.dat']
data1503=[]
files1506=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=1.50_a= 6.15_N=8.dat']
data1506=[]

for data_file in files15030:
    data15030.append(np.loadtxt(data_file))
for data_file in files1503:
    data1503.append(np.loadtxt(data_file))
for data_file in files1506:
    data1506.append(np.loadtxt(data_file))

Res15030 = data15030[0]
Res1503 = data1503[0]
Res1506 = data1506[0]


files17530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=1.75_a=30.44_N=8.dat']
data17530=[]
files1753=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=1.75_a= 3.06_N=8.dat']
data1753=[]
files1756=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=1.75_a= 6.15_N=8.dat']
data1756=[]

for data_file in files17530:
    data17530.append(np.loadtxt(data_file))
for data_file in files1753:
    data1753.append(np.loadtxt(data_file))
for data_file in files1756:
    data1756.append(np.loadtxt(data_file))

Res17530 = data17530[0]
Res1753 = data1753[0]
Res1756 = data1756[0]

#%%

files230=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=2.00_a=30.44_N=8.dat']
data230=[]
files23=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=2.00_a= 3.06_N=8.dat']
data23=[]
files26=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=2.00_a= 6.15_N=8.dat']
data26=[]

for data_file in files230:
    data230.append(np.loadtxt(data_file))
for data_file in files23:
    data23.append(np.loadtxt(data_file))
for data_file in files26:
    data26.append(np.loadtxt(data_file))

Res230 = data230[0]
Res23 = data23[0]
Res26 = data26[0]

plt.plot(Res230[0:l,0],Res230[0:l,1])
plt.plot(Res23[0:l,0],Res23[0:l,1])
plt.plot(Res26[0:l,0],Res26[0:l,1])


files22530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=2.25_a=30.44_N=8.dat']
data22530=[]
files2253=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=2.25_a= 3.06_N=8.dat']
data2253=[]
files2256=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=2.25_a= 6.15_N=8.dat']
data2256=[]

for data_file in files22530:
    data22530.append(np.loadtxt(data_file))
for data_file in files2253:
    data2253.append(np.loadtxt(data_file))
for data_file in files2256:
    data2256.append(np.loadtxt(data_file))

Res22530 = data22530[0]
Res2253 = data2253[0]
Res2256 = data2256[0]


files25030=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=2.50_a=30.44_N=8.dat']
data25030=[]
files2503=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=2.50_a= 3.06_N=8.dat']
data2503=[]
files2506=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=2.50_a= 6.15_N=8.dat']
data2506=[]

for data_file in files25030:
    data25030.append(np.loadtxt(data_file))
for data_file in files2503:
    data2503.append(np.loadtxt(data_file))
for data_file in files2506:
    data2506.append(np.loadtxt(data_file))

Res25030 = data25030[0]
Res2503 = data2503[0]
Res2506 = data2506[0]


files27530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=2.75_a=30.44_N=8.dat']
data27530=[]
files2753=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=2.75_a= 3.06_N=8.dat']
data2753=[]
files2756=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=2.75_a= 6.15_N=8.dat']
data2756=[]

for data_file in files27530:
    data27530.append(np.loadtxt(data_file))
for data_file in files2753:
    data2753.append(np.loadtxt(data_file))
for data_file in files2756:
    data2756.append(np.loadtxt(data_file))

Res27530 = data27530[0]
Res2753 = data2753[0]
Res2756 = data2756[0]


#%%

files330=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=3.00_a=30.44_N=8.dat']
data330=[]
files33=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=3.00_a= 3.06_N=8.dat']
data33=[]
files36=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=3.00_a= 6.15_N=8.dat']
data36=[]

for data_file in files330:
    data330.append(np.loadtxt(data_file))
for data_file in files33:
    data33.append(np.loadtxt(data_file))
for data_file in files36:
    data36.append(np.loadtxt(data_file))

Res330 = data330[0]
Res33 = data33[0]
Res36 = data36[0]

plt.plot(Res330[0:l,0],Res330[0:l,1])
plt.plot(Res33[0:l,0],Res33[0:l,1])
plt.plot(Res36[0:l,0],Res36[0:l,1])



files32530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=3.25_a=30.44_N=8.dat']
data32530=[]
files3253=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=3.25_a= 3.06_N=8.dat']
data3253=[]
files3256=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=3.25_a= 6.15_N=8.dat']
data3256=[]

for data_file in files32530:
    data32530.append(np.loadtxt(data_file))
for data_file in files3253:
    data3253.append(np.loadtxt(data_file))
for data_file in files3256:
    data3256.append(np.loadtxt(data_file))

Res32530 = data32530[0]
Res3253 = data3253[0]
Res3256 = data3256[0]



files35030=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=3.50_a=30.44_N=8.dat']
data35030=[]
files3503=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=3.50_a= 3.06_N=8.dat']
data3503=[]
files3506=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=3.50_a= 6.15_N=8.dat']
data3506=[]

for data_file in files35030:
    data35030.append(np.loadtxt(data_file))
for data_file in files3503:
    data3503.append(np.loadtxt(data_file))
for data_file in files3506:
    data3506.append(np.loadtxt(data_file))

Res35030 = data35030[0]
Res3503 = data3503[0]
Res3506 = data3506[0]



files37530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=3.75_a=30.44_N=8.dat']
data37530=[]
files3753=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=3.75_a= 3.06_N=8.dat']
data3753=[]
files3756=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=3.75_a= 6.15_N=8.dat']
data3756=[]

for data_file in files37530:
    data37530.append(np.loadtxt(data_file))
for data_file in files3753:
    data3753.append(np.loadtxt(data_file))
for data_file in files3756:
    data3756.append(np.loadtxt(data_file))

Res37530 = data37530[0]
Res3753 = data3753[0]
Res3756 = data3756[0]

#%%

files430=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=4.00_a=30.44_N=8.dat']
data430=[]
files43=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=4.00_a= 3.06_N=8.dat']
data43=[]
files46=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=4.00_a= 6.15_N=8.dat']
data46=[]

for data_file in files430:
    data430.append(np.loadtxt(data_file))
for data_file in files43:
    data43.append(np.loadtxt(data_file))
for data_file in files46:
    data46.append(np.loadtxt(data_file))

Res430 = data430[0]
Res43 = data43[0]
Res46 = data46[0]



files42530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=4.25_a=30.44_N=8.dat']
data42530=[]
files4253=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=4.25_a= 3.06_N=8.dat']
data4253=[]
files4256=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=4.25_a= 6.15_N=8.dat']
data4256=[]

for data_file in files42530:
    data42530.append(np.loadtxt(data_file))
for data_file in files4253:
    data4253.append(np.loadtxt(data_file))
for data_file in files4256:
    data4256.append(np.loadtxt(data_file))

Res42530 = data42530[0]
Res4253 = data4253[0]
Res4256 = data4256[0]



files45030=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=4.50_a=30.44_N=8.dat']
data45030=[]
files4503=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=4.50_a= 3.06_N=8.dat']
data4503=[]
files4506=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=4.50_a= 6.15_N=8.dat']
data4506=[]

for data_file in files45030:
    data45030.append(np.loadtxt(data_file))
for data_file in files4503:
    data4503.append(np.loadtxt(data_file))
for data_file in files4506:
    data4506.append(np.loadtxt(data_file))

Res45030 = data45030[0]
Res4503 = data4503[0]
Res4506 = data4506[0]



files47530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=4.75_a=30.44_N=8.dat']
data47530=[]
files4753=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=4.75_a= 3.06_N=8.dat']
data4753=[]
files4756=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=4.75_a= 6.15_N=8.dat']
data4756=[]

for data_file in files47530:
    data47530.append(np.loadtxt(data_file))
for data_file in files4753:
    data4753.append(np.loadtxt(data_file))
for data_file in files4756:
    data4756.append(np.loadtxt(data_file))

Res47530 = data47530[0]
Res4753 = data4753[0]
Res4756 = data4756[0]



files530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=5.00_a=30.44_N=8.dat']
data530=[]
files53=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=5.00_a= 3.06_N=8.dat']
data53=[]
files56=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg_J=5.00_a= 6.15_N=8.dat']
data56=[]

for data_file in files530:
    data530.append(np.loadtxt(data_file))
for data_file in files53:
    data53.append(np.loadtxt(data_file))
for data_file in files56:
    data56.append(np.loadtxt(data_file))

Res530 = data530[0]
Res53 = data53[0]
Res56 = data56[0]

plt.plot(Res530[0:l,0],Res530[0:l,1])
plt.plot(Res53[0:l,0],Res53[0:l,1])
plt.plot(Res56[0:l,0],Res56[0:l,1])

#%%

n0=0.0

tRes = np.zeros((3,21))

poRes030, pcRes030 = curve_fit(N,Res030[0:l,0],Res030[0:l,1])
poRes03, pcRes03 = curve_fit(N,Res03[0:l,0],Res03[0:l,1])
poRes06, pcRes06 = curve_fit(N,Res06[0:l,0],Res06[0:l,1])

poRes02530, pcRes02530 = curve_fit(N,Res02530[0:l,0],Res02530[0:l,1])
poRes0253, pcRes0253 = curve_fit(N,Res0253[0:l,0],Res0253[0:l,1])
poRes0256, pcRes0256 = curve_fit(N,Res0256[0:l,0],Res0256[0:l,1])

poRes05030, pcRes05030 = curve_fit(N,Res05030[0:l,0],Res05030[0:l,1])
poRes0503, pcRes0503 = curve_fit(N,Res0503[0:l,0],Res0503[0:l,1])
poRes0506, pcRes0506 = curve_fit(N,Res0506[0:l,0],Res0506[0:l,1])

poRes07530, pcRes07530 = curve_fit(N,Res07530[0:l,0],Res07530[0:l,1])
poRes0753, pcRes0753 = curve_fit(N,Res0753[0:l,0],Res0753[0:l,1])
poRes0756, pcRes0756 = curve_fit(N,Res0756[0:l,0],Res0756[0:l,1])



poRes130, pcRes130 = curve_fit(N,Res130[0:l,0],Res130[0:l,1])
poRes13, pcRes13 = curve_fit(N,Res13[0:l,0],Res13[0:l,1])
poRes16, pcRes16 = curve_fit(N,Res16[0:l,0],Res16[0:l,1])

poRes17530, pcRes17530 = curve_fit(N,Res17530[0:l,0],Res17530[0:l,1])
poRes1753, pcRes1753 = curve_fit(N,Res1753[0:l,0],Res1753[0:l,1])
poRes1756, pcRes1756 = curve_fit(N,Res1756[0:l,0],Res1756[0:l,1])

poRes12530, pcRes12530 = curve_fit(N,Res12530[0:l,0],Res12530[0:l,1])
poRes1253, pcRes1253 = curve_fit(N,Res1253[0:l,0],Res1253[0:l,1])
poRes1256, pcRes1256 = curve_fit(N,Res1256[0:l,0],Res1256[0:l,1])

poRes15030, pcRes15030 = curve_fit(N,Res15030[0:l,0],Res15030[0:l,1])
poRes1503, pcRes1503 = curve_fit(N,Res1503[0:l,0],Res1503[0:l,1])
poRes1506, pcRes1506 = curve_fit(N,Res1506[0:l,0],Res1506[0:l,1])




poRes230, pcRes230 = curve_fit(N,Res230[0:l,0],Res230[0:l,1])
poRes23, pcRes23 = curve_fit(N,Res23[0:l,0],Res23[0:l,1])
poRes26, pcRes26 = curve_fit(N,Res26[0:l,0],Res26[0:l,1])

poRes22530, pcRes230 = curve_fit(N,Res22530[0:l,0],Res22530[0:l,1])
poRes2253, pcRes2253 = curve_fit(N,Res2253[0:l,0],Res2253[0:l,1])
poRes2256, pcRes2256 = curve_fit(N,Res2256[0:l,0],Res2256[0:l,1])

poRes25030, pcRes25030 = curve_fit(N,Res25030[0:l,0],Res25030[0:l,1])
poRes2503, pcRes2503 = curve_fit(N,Res2503[0:l,0],Res2503[0:l,1])
poRes2506, pcRes2506 = curve_fit(N,Res2506[0:l,0],Res2506[0:l,1])

poRes27530, pcRes27530 = curve_fit(N,Res27530[0:l,0],Res27530[0:l,1])
poRes2753, pcRes2753 = curve_fit(N,Res2753[0:l,0],Res2753[0:l,1])
poRes2756, pcRes2756 = curve_fit(N,Res2756[0:l,0],Res2756[0:l,1])



poRes330, pcRes330 = curve_fit(N,Res330[0:l,0],Res330[0:l,1])
poRes33, pcRes33 = curve_fit(N,Res33[0:l,0],Res33[0:l,1])
poRes36, pcRes36 = curve_fit(N,Res36[0:l,0],Res36[0:l,1])

poRes32530, pcRes32530 = curve_fit(N,Res32530[0:l,0],Res32530[0:l,1])
poRes3253, pcRes3253 = curve_fit(N,Res3253[0:l,0],Res3253[0:l,1])
poRes3256, pcRes3256 = curve_fit(N,Res3256[0:l,0],Res3256[0:l,1])

poRes35030, pcRes35030 = curve_fit(N,Res35030[0:l,0],Res35030[0:l,1])
poRes3503, pcRes3503 = curve_fit(N,Res3503[0:l,0],Res3503[0:l,1])
poRes3506, pcRes3506 = curve_fit(N,Res3506[0:l,0],Res3506[0:l,1])

poRes37530, pcRes37530 = curve_fit(N,Res37530[0:l,0],Res37530[0:l,1])
poRes3753, pcRes3753 = curve_fit(N,Res3753[0:l,0],Res3753[0:l,1])
poRes3756, pcRes3756 = curve_fit(N,Res3756[0:l,0],Res3756[0:l,1])



poRes430, pcRes430 = curve_fit(N,Res430[0:l,0],Res430[0:l,1])
poRes43, pcRes43 = curve_fit(N,Res43[0:l,0],Res43[0:l,1])
poRes46, pcRes46 = curve_fit(N,Res46[0:l,0],Res46[0:l,1])

poRes42530, pcRes42530 = curve_fit(N,Res42530[0:l,0],Res42530[0:l,1])
poRes4253, pcRes4253 = curve_fit(N,Res4253[0:l,0],Res4253[0:l,1])
poRes4256, pcRes4256 = curve_fit(N,Res4256[0:l,0],Res4256[0:l,1])

poRes45030, pcRes45030 = curve_fit(N,Res45030[0:l,0],Res45030[0:l,1])
poRes4503, pcRes4503 = curve_fit(N,Res4503[0:l,0],Res4503[0:l,1])
poRes4506, pcRes4506 = curve_fit(N,Res4506[0:l,0],Res4506[0:l,1])

poRes47530, pcRes47530 = curve_fit(N,Res47530[0:l,0],Res47530[0:l,1])
poRes4753, pcRes4753 = curve_fit(N,Res4753[0:l,0],Res4753[0:l,1])
poRes4756, pcRes4756 = curve_fit(N,Res4756[0:l,0],Res4756[0:l,1])


poRes530, pcRes530 = curve_fit(N,Res530[0:l,0],Res530[0:l,1])
poRes53, pcRes53 = curve_fit(N,Res53[0:l,0],Res53[0:l,1])
poRes56, pcRes56 = curve_fit(N,Res56[0:l,0],Res56[0:l,1])

#%%

tRes[0,0] = poRes03[0]
tRes[1,0] = poRes06[0]
tRes[2,0] = poRes030[0]

tRes[0,1] = poRes0253[0]
tRes[1,1] = poRes0256[0]
tRes[2,1] = poRes02530[0]

tRes[0,2] = poRes0503[0]
tRes[1,2] = poRes0506[0]
tRes[2,2] = poRes05030[0]

tRes[0,3] = poRes0753[0]
tRes[1,3] = poRes0756[0]
tRes[2,3] = poRes07530[0]

tRes[0,4] = poRes13[0]
tRes[1,4] = poRes16[0]
tRes[2,4] = poRes130[0]

tRes[0,5] = poRes1253[0]
tRes[1,5] = poRes1256[0]
tRes[2,5] = poRes12530[0]

tRes[0,6] = poRes1503[0]
tRes[1,6] = poRes1506[0]
tRes[2,6] = poRes15030[0]

tRes[0,7] = poRes1753[0]
tRes[1,7] = poRes1756[0]
tRes[2,7] = poRes17530[0]

tRes[0,8] = poRes23[0]
tRes[1,8] = poRes26[0]
tRes[2,8] = poRes230[0]

tRes[0,9] = poRes2253[0]
tRes[1,9] = poRes2256[0]
tRes[2,9] = poRes22530[0]

tRes[0,10] = poRes2503[0]
tRes[1,10] = poRes2506[0]
tRes[2,10] = poRes25030[0]

tRes[0,11] = poRes2753[0]
tRes[1,11] = poRes2756[0]
tRes[2,11] = poRes27530[0]

tRes[0,12] = poRes33[0]
tRes[1,12] = poRes36[0]
tRes[2,12] = poRes330[0]

tRes[0,13] = poRes3253[0]
tRes[1,13] = poRes3256[0]
tRes[2,13] = poRes32530[0]

tRes[0,14] = poRes3503[0]
tRes[1,14] = poRes3506[0]
tRes[2,14] = poRes35030[0]

tRes[0,15] = poRes3753[0]
tRes[1,15] = poRes3756[0]
tRes[2,15] = poRes37530[0]

tRes[0,16] = poRes43[0]
tRes[1,16] = poRes46[0]
tRes[2,16] = poRes430[0]

tRes[0,17] = poRes4253[0]
tRes[1,17] = poRes4256[0]
tRes[2,17] = poRes42530[0]

tRes[0,18] = poRes4503[0]
tRes[1,18] = poRes4506[0]
tRes[2,18] = poRes45030[0]

tRes[0,19] = poRes4753[0]
tRes[1,19] = poRes4756[0]
tRes[2,19] = poRes47530[0]

tRes[0,20] = poRes53[0]
tRes[1,20] = poRes56[0]
tRes[2,20] = poRes530[0]


tRes

#%%

m=21

J=np.arange(0,5.25,0.25)


plt.rcParams["figure.figsize"] =[8.0,6.0]
plt.title('Resurrection relax times',fontsize=16)
plt.xlabel('J  ($k_BT$)', fontsize=16)
plt.ylabel(r'$\tau_{res}$ (s)',fontsize=16)
plt.grid()


plt.plot(J,tRes[0,:],marker='.',linestyle="",markersize=7,label=r'$\alpha=3.06$')
plt.plot(J,tRes[1,:],marker='.',linestyle="",markersize=7,label=r'$\alpha=6.15$')
plt.plot(J,tRes[2,:],marker='.',linestyle="",markersize=7,label=r'$\alpha=30.44$')

plt.legend(loc='upper left',prop={'size': 16})


#%%

filesS030=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=0.00_a=30.44_N=8.dat']
StallD030=[]
filesS03=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=0.00_a= 3.06_N=8.dat']
StallD03=[]
filesS06=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=0.00_a= 6.15_N=8.dat']
StallD06=[]

for data_file in filesS030:
    StallD030.append(np.loadtxt(data_file))
for data_file in filesS03:
    StallD03.append(np.loadtxt(data_file))
for data_file in filesS06:
    StallD06.append(np.loadtxt(data_file))

Stall030 = StallD030[0]
Stall03 = StallD03[0]
Stall06 = StallD06[0]


plt.plot(Stall030[0:l,0],Stall030[0:l,1])
plt.plot(Stall03[0:l,0],Stall03[0:l,1])
plt.plot(Stall06[0:l,0],Stall06[0:l,1])

filesS02530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=0.25_a=30.44_N=8.dat']
StallD02530=[]
filesS0253=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=0.25_a= 3.06_N=8.dat']
StallD0253=[]
filesS0256=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=0.25_a= 6.15_N=8.dat']
StallD0256=[]

for data_file in filesS02530:
    StallD02530.append(np.loadtxt(data_file))
for data_file in filesS0253:
    StallD0253.append(np.loadtxt(data_file))
for data_file in filesS0256:
    StallD0256.append(np.loadtxt(data_file))

Stall02530 = StallD02530[0]
Stall0253 = StallD0253[0]
Stall0256 = StallD0256[0]


filesS05030=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=0.50_a=30.44_N=8.dat']
StallD05030=[]
filesS0503=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=0.50_a= 3.06_N=8.dat']
StallD0503=[]
filesS0506=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=0.50_a= 6.15_N=8.dat']
StallD0506=[]

for data_file in filesS05030:
    StallD05030.append(np.loadtxt(data_file))
for data_file in filesS0503:
    StallD0503.append(np.loadtxt(data_file))
for data_file in filesS0506:
    StallD0506.append(np.loadtxt(data_file))

Stall05030 = StallD05030[0]
Stall0503 = StallD0503[0]
Stall0506 = StallD0506[0]


filesS07530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=0.75_a=30.44_N=8.dat']
StallD07530=[]
filesS0753=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=0.75_a= 3.06_N=8.dat']
StallD0753=[]
filesS0756=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=0.75_a= 6.15_N=8.dat']
StallD0756=[]

for data_file in filesS07530:
    StallD07530.append(np.loadtxt(data_file))
for data_file in filesS0753:
    StallD0753.append(np.loadtxt(data_file))
for data_file in filesS0756:
    StallD0756.append(np.loadtxt(data_file))

Stall07530 = StallD07530[0]
Stall0753 = StallD0753[0]
Stall0756 = StallD0756[0]


#%%
filesS130=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=1.00_a=30.44_N=8.dat']
StallD130=[]
filesS13=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=1.00_a= 3.06_N=8.dat']
StallD13=[]
filesS16=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=1.00_a= 6.15_N=8.dat']
StallD16=[]

for data_file in filesS130:
    StallD130.append(np.loadtxt(data_file))
for data_file in filesS13:
    StallD13.append(np.loadtxt(data_file))
for data_file in filesS16:
    StallD16.append(np.loadtxt(data_file))

Stall130 = StallD130[0]
Stall13 = StallD13[0]
Stall16 = StallD16[0]

plt.plot(Stall130[0:l,0],Stall130[0:l,1])
plt.plot(Stall13[0:l,0],Stall13[0:l,1])
plt.plot(Stall16[0:l,0],Stall16[0:l,1])


filesS12530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=1.25_a=30.44_N=8.dat']
StallD12530=[]
filesS1253=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=1.25_a= 3.06_N=8.dat']
StallD1253=[]
filesS1256=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=1.25_a= 6.15_N=8.dat']
StallD1256=[]

for data_file in filesS12530:
    StallD12530.append(np.loadtxt(data_file))
for data_file in filesS1253:
    StallD1253.append(np.loadtxt(data_file))
for data_file in filesS1256:
    StallD1256.append(np.loadtxt(data_file))

Stall12530 = StallD12530[0]
Stall1253 = StallD1253[0]
Stall1256 = StallD1256[0]


filesS15030=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=1.50_a=30.44_N=8.dat']
StallD15030=[]
filesS1503=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=1.50_a= 3.06_N=8.dat']
StallD1503=[]
filesS1506=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=1.50_a= 6.15_N=8.dat']
StallD1506=[]

for data_file in filesS15030:
    StallD15030.append(np.loadtxt(data_file))
for data_file in filesS1503:
    StallD1503.append(np.loadtxt(data_file))
for data_file in filesS1506:
    StallD1506.append(np.loadtxt(data_file))

Stall15030 = StallD15030[0]
Stall1503 = StallD1503[0]
Stall1506 = StallD1506[0]


filesS17530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=1.75_a=30.44_N=8.dat']
StallD17530=[]
filesS1753=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=1.75_a= 3.06_N=8.dat']
StallD1753=[]
filesS1756=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=1.75_a= 6.15_N=8.dat']
StallD1756=[]

for data_file in filesS17530:
    StallD17530.append(np.loadtxt(data_file))
for data_file in filesS1753:
    StallD1753.append(np.loadtxt(data_file))
for data_file in filesS1756:
    StallD1756.append(np.loadtxt(data_file))

Stall17530 = StallD17530[0]
Stall1753 = StallD1753[0]
Stall1756 = StallD1756[0]

#%%

filesS230=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=2.00_a=30.44_N=8.dat']
StallD230=[]
filesS23=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=2.00_a= 3.06_N=8.dat']
StallD23=[]
filesS26=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=2.00_a= 6.15_N=8.dat']
StallD26=[]

for data_file in filesS230:
    StallD230.append(np.loadtxt(data_file))
for data_file in filesS23:
    StallD23.append(np.loadtxt(data_file))
for data_file in filesS26:
    StallD26.append(np.loadtxt(data_file))

Stall230 = StallD230[0]
Stall23 = StallD23[0]
Stall26 = StallD26[0]

plt.plot(Stall230[0:l,0],Stall230[0:l,1])
plt.plot(Stall23[0:l,0],Stall23[0:l,1])
plt.plot(Stall26[0:l,0],Stall26[0:l,1])


filesS22530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=2.25_a=30.44_N=8.dat']
StallD22530=[]
filesS2253=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=2.25_a= 3.06_N=8.dat']
StallD2253=[]
filesS2256=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=2.25_a= 6.15_N=8.dat']
StallD2256=[]

for data_file in filesS22530:
    StallD22530.append(np.loadtxt(data_file))
for data_file in filesS2253:
    StallD2253.append(np.loadtxt(data_file))
for data_file in filesS2256:
    StallD2256.append(np.loadtxt(data_file))

Stall22530 = StallD22530[0]
Stall2253 = StallD2253[0]
Stall2256 = StallD2256[0]


filesS25030=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=2.50_a=30.44_N=8.dat']
StallD25030=[]
filesS2503=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=2.50_a= 3.06_N=8.dat']
StallD2503=[]
filesS2506=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=2.50_a= 6.15_N=8.dat']
StallD2506=[]

for data_file in filesS25030:
    StallD25030.append(np.loadtxt(data_file))
for data_file in filesS2503:
    StallD2503.append(np.loadtxt(data_file))
for data_file in filesS2506:
    StallD2506.append(np.loadtxt(data_file))

Stall25030 = StallD25030[0]
Stall2503 = StallD2503[0]
Stall2506 = StallD2506[0]


filesS27530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=2.75_a=30.44_N=8.dat']
StallD27530=[]
filesS2753=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=2.75_a= 3.06_N=8.dat']
StallD2753=[]
filesS2756=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=2.75_a= 6.15_N=8.dat']
StallD2756=[]

for data_file in filesS27530:
    StallD27530.append(np.loadtxt(data_file))
for data_file in filesS2753:
    StallD2753.append(np.loadtxt(data_file))
for data_file in filesS2756:
    StallD2756.append(np.loadtxt(data_file))

Stall27530 = StallD27530[0]
Stall2753 = StallD2753[0]
Stall2756 = StallD2756[0]


#%%

filesS330=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=3.00_a=30.44_N=8.dat']
StallD330=[]
filesS33=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=3.00_a= 3.06_N=8.dat']
StallD33=[]
filesS36=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=3.00_a= 6.15_N=8.dat']
StallD36=[]

for data_file in filesS330:
    StallD330.append(np.loadtxt(data_file))
for data_file in filesS33:
    StallD33.append(np.loadtxt(data_file))
for data_file in filesS36:
    StallD36.append(np.loadtxt(data_file))

Stall330 = StallD330[0]
Stall33 = StallD33[0]
Stall36 = StallD36[0]

plt.plot(Stall330[0:l,0],Stall330[0:l,1])
plt.plot(Stall33[0:l,0],Stall33[0:l,1])
plt.plot(Stall36[0:l,0],Stall36[0:l,1])



filesS32530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=3.25_a=30.44_N=8.dat']
StallD32530=[]
filesS3253=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=3.25_a= 3.06_N=8.dat']
StallD3253=[]
filesS3256=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=3.25_a= 6.15_N=8.dat']
StallD3256=[]

for data_file in filesS32530:
    StallD32530.append(np.loadtxt(data_file))
for data_file in filesS3253:
    StallD3253.append(np.loadtxt(data_file))
for data_file in filesS3256:
    StallD3256.append(np.loadtxt(data_file))

Stall32530 = StallD32530[0]
Stall3253 = StallD3253[0]
Stall3256 = StallD3256[0]



filesS35030=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=3.50_a=30.44_N=8.dat']
StallD35030=[]
filesS3503=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=3.50_a= 3.06_N=8.dat']
StallD3503=[]
filesS3506=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=3.50_a= 6.15_N=8.dat']
StallD3506=[]

for data_file in filesS35030:
    StallD35030.append(np.loadtxt(data_file))
for data_file in filesS3503:
    StallD3503.append(np.loadtxt(data_file))
for data_file in filesS3506:
    StallD3506.append(np.loadtxt(data_file))

Stall35030 = StallD35030[0]
Stall3503 = StallD3503[0]
Stall3506 = StallD3506[0]



filesS37530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=3.75_a=30.44_N=8.dat']
StallD37530=[]
filesS3753=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=3.75_a= 3.06_N=8.dat']
StallD3753=[]
filesS3756=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=3.75_a= 6.15_N=8.dat']
StallD3756=[]

for data_file in filesS37530:
    StallD37530.append(np.loadtxt(data_file))
for data_file in filesS3753:
    StallD3753.append(np.loadtxt(data_file))
for data_file in filesS3756:
    StallD3756.append(np.loadtxt(data_file))

Stall37530 = StallD37530[0]
Stall3753 = StallD3753[0]
Stall3756 = StallD3756[0]

#%%

filesS430=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=4.00_a=30.44_N=8.dat']
StallD430=[]
filesS43=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=4.00_a= 3.06_N=8.dat']
StallD43=[]
filesS46=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=4.00_a= 6.15_N=8.dat']
StallD46=[]

for data_file in filesS430:
    StallD430.append(np.loadtxt(data_file))
for data_file in filesS43:
    StallD43.append(np.loadtxt(data_file))
for data_file in filesS46:
    StallD46.append(np.loadtxt(data_file))

Stall430 = StallD430[0]
Stall43 = StallD43[0]
Stall46 = StallD46[0]




filesS42530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=4.25_a=30.44_N=8.dat']
StallD42530=[]
filesS4253=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=4.25_a= 3.06_N=8.dat']
StallD4253=[]
filesS4256=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=4.25_a= 6.15_N=8.dat']
StallD4256=[]

for data_file in filesS42530:
    StallD42530.append(np.loadtxt(data_file))
for data_file in filesS4253:
    StallD4253.append(np.loadtxt(data_file))
for data_file in filesS4256:
    StallD4256.append(np.loadtxt(data_file))

Stall42530 = StallD42530[0]
Stall4253 = StallD4253[0]
Stall4256 = StallD4256[0]



filesS45030=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=4.50_a=30.44_N=8.dat']
StallD45030=[]
filesS4503=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=4.50_a= 3.06_N=8.dat']
StallD4503=[]
filesS4506=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=4.50_a= 6.15_N=8.dat']
StallD4506=[]

for data_file in filesS45030:
    StallD45030.append(np.loadtxt(data_file))
for data_file in filesS4503:
    StallD4503.append(np.loadtxt(data_file))
for data_file in filesS4506:
    StallD4506.append(np.loadtxt(data_file))

Stall45030 = StallD45030[0]
Stall4503 = StallD4503[0]
Stall4506 = StallD4506[0]



filesS47530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=4.75_a=30.44_N=8.dat']
StallD47530=[]
filesS4753=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=4.75_a= 3.06_N=8.dat']
StallD4753=[]
filesS4756=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=4.75_a= 6.15_N=8.dat']
StallD4756=[]

for data_file in filesS47530:
    StallD47530.append(np.loadtxt(data_file))
for data_file in filesS4753:
    StallD4753.append(np.loadtxt(data_file))
for data_file in filesS4756:
    StallD4756.append(np.loadtxt(data_file))

Stall47530 = StallD47530[0]
Stall4753 = StallD4753[0]
Stall4756 = StallD4756[0]



filesS530=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=5.00_a=30.44_N=8.dat']
StallD530=[]
filesS53=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=5.00_a= 3.06_N=8.dat']
StallD53=[]
filesS56=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg_J=5.00_a= 6.15_N=8.dat']
StallD56=[]

for data_file in filesS530:
    StallD530.append(np.loadtxt(data_file))
for data_file in filesS53:
    StallD53.append(np.loadtxt(data_file))
for data_file in filesS56:
    StallD56.append(np.loadtxt(data_file))

Stall530 = StallD530[0]
Stall53 = StallD53[0]
Stall56 = StallD56[0]

plt.plot(Stall530[0:l,0],Stall530[0:l,1])
plt.plot(Stall53[0:l,0],Stall53[0:l,1])
plt.plot(Stall56[0:l,0],Stall56[0:l,1])

#%%

n0=9.0

tStall = np.zeros((3,21))

poStall030, pcStall030 = curve_fit(N,Stall030[0:l,0],Stall030[0:l,1])
poStall03, pcStall03 = curve_fit(N,Stall03[0:l,0],Stall03[0:l,1])
poStall06, pcStall06 = curve_fit(N,Stall06[0:l,0],Stall06[0:l,1])

poStall02530, pcStall02530 = curve_fit(N,Stall02530[0:l,0],Stall02530[0:l,1])
poStall0253, pcStall0253 = curve_fit(N,Stall0253[0:l,0],Stall0253[0:l,1])
poStall0256, pcStall0256 = curve_fit(N,Stall0256[0:l,0],Stall0256[0:l,1])

poStall05030, pcStall05030 = curve_fit(N,Stall05030[0:l,0],Stall05030[0:l,1])
poStall0503, pcStall0503 = curve_fit(N,Stall0503[0:l,0],Stall0503[0:l,1])
poStall0506, pcStall0506 = curve_fit(N,Stall0506[0:l,0],Stall0506[0:l,1])

poStall07530, pcStall07530 = curve_fit(N,Stall07530[0:l,0],Stall07530[0:l,1])
poStall0753, pcStall0753 = curve_fit(N,Stall0753[0:l,0],Stall0753[0:l,1])
poStall0756, pcStall0756 = curve_fit(N,Stall0756[0:l,0],Stall0756[0:l,1])



poStall130, pcStall130 = curve_fit(N,Stall130[0:l,0],Stall130[0:l,1])
poStall13, pcStall13 = curve_fit(N,Stall13[0:l,0],Stall13[0:l,1])
poStall16, pcStall16 = curve_fit(N,Stall16[0:l,0],Stall16[0:l,1])

poStall17530, pcStall17530 = curve_fit(N,Stall17530[0:l,0],Stall17530[0:l,1])
poStall1753, pcStall1753 = curve_fit(N,Stall1753[0:l,0],Stall1753[0:l,1])
poStall1756, pcStall1756 = curve_fit(N,Stall1756[0:l,0],Stall1756[0:l,1])

poStall12530, pcStall12530 = curve_fit(N,Stall12530[0:l,0],Stall12530[0:l,1])
poStall1253, pcStall1253 = curve_fit(N,Stall1253[0:l,0],Stall1253[0:l,1])
poStall1256, pcStall1256 = curve_fit(N,Stall1256[0:l,0],Stall1256[0:l,1])

poStall15030, pcStall15030 = curve_fit(N,Stall15030[0:l,0],Stall15030[0:l,1])
poStall1503, pcStall1503 = curve_fit(N,Stall1503[0:l,0],Stall1503[0:l,1])
poStall1506, pcStall1506 = curve_fit(N,Stall1506[0:l,0],Stall1506[0:l,1])




poStall230, pcStall230 = curve_fit(N,Stall230[0:l,0],Stall230[0:l,1])
poStall23, pcStall23 = curve_fit(N,Stall23[0:l,0],Stall23[0:l,1])
poStall26, pcStall26 = curve_fit(N,Stall26[0:l,0],Stall26[0:l,1])

poStall22530, pcStall230 = curve_fit(N,Stall22530[0:l,0],Stall22530[0:l,1])
poStall2253, pcStall2253 = curve_fit(N,Stall2253[0:l,0],Stall2253[0:l,1])
poStall2256, pcStall2256 = curve_fit(N,Stall2256[0:l,0],Stall2256[0:l,1])

poStall25030, pcStall25030 = curve_fit(N,Stall25030[0:l,0],Stall25030[0:l,1])
poStall2503, pcStall2503 = curve_fit(N,Stall2503[0:l,0],Stall2503[0:l,1])
poStall2506, pcStall2506 = curve_fit(N,Stall2506[0:l,0],Stall2506[0:l,1])

poStall27530, pcStall27530 = curve_fit(N,Stall27530[0:l,0],Stall27530[0:l,1])
poStall2753, pcStall2753 = curve_fit(N,Stall2753[0:l,0],Stall2753[0:l,1])
poStall2756, pcStall2756 = curve_fit(N,Stall2756[0:l,0],Stall2756[0:l,1])



poStall330, pcStall330 = curve_fit(N,Stall330[0:l,0],Stall330[0:l,1])
poStall33, pcStall33 = curve_fit(N,Stall33[0:l,0],Stall33[0:l,1])
poStall36, pcStall36 = curve_fit(N,Stall36[0:l,0],Stall36[0:l,1])

poStall32530, pcStall32530 = curve_fit(N,Stall32530[0:l,0],Stall32530[0:l,1])
poStall3253, pcStall3253 = curve_fit(N,Stall3253[0:l,0],Stall3253[0:l,1])
poStall3256, pcStall3256 = curve_fit(N,Stall3256[0:l,0],Stall3256[0:l,1])

poStall35030, pcStall35030 = curve_fit(N,Stall35030[0:l,0],Stall35030[0:l,1])
poStall3503, pcStall3503 = curve_fit(N,Stall3503[0:l,0],Stall3503[0:l,1])
poStall3506, pcStall3506 = curve_fit(N,Stall3506[0:l,0],Stall3506[0:l,1])

poStall37530, pcStall37530 = curve_fit(N,Stall37530[0:l,0],Stall37530[0:l,1])
poStall3753, pcStall3753 = curve_fit(N,Stall3753[0:l,0],Stall3753[0:l,1])
poStall3756, pcStall3756 = curve_fit(N,Stall3756[0:l,0],Stall3756[0:l,1])



poStall430, pcStall430 = curve_fit(N,Stall430[0:l,0],Stall430[0:l,1])
poStall43, pcStall43 = curve_fit(N,Stall43[0:l,0],Stall43[0:l,1])
poStall46, pcStall46 = curve_fit(N,Stall46[0:l,0],Stall46[0:l,1])

poStall42530, pcStall42530 = curve_fit(N,Stall42530[0:l,0],Stall42530[0:l,1])
poStall4253, pcStall4253 = curve_fit(N,Stall4253[0:l,0],Stall4253[0:l,1])
poStall4256, pcStall4256 = curve_fit(N,Stall4256[0:l,0],Stall4256[0:l,1])

poStall45030, pcStall45030 = curve_fit(N,Stall45030[0:l,0],Stall45030[0:l,1])
poStall4503, pcStall4503 = curve_fit(N,Stall4503[0:l,0],Stall4503[0:l,1])
poStall4506, pcStall4506 = curve_fit(N,Stall4506[0:l,0],Stall4506[0:l,1])

poStall47530, pcStall47530 = curve_fit(N,Stall47530[0:l,0],Stall47530[0:l,1])
poStall4753, pcStall4753 = curve_fit(N,Stall4753[0:l,0],Stall4753[0:l,1])
poStall4756, pcStall4756 = curve_fit(N,Stall4756[0:l,0],Stall4756[0:l,1])


poStall530, pcStall530 = curve_fit(N,Stall530[0:l,0],Stall530[0:l,1])
poStall53, pcStall53 = curve_fit(N,Stall53[0:l,0],Stall53[0:l,1])
poStall56, pcStall56 = curve_fit(N,Stall56[0:l,0],Stall56[0:l,1])

#%%

tStall[0,0] = poStall03[0]
tStall[1,0] = poStall06[0]
tStall[2,0] = poStall030[0]

tStall[0,1] = poStall0253[0]
tStall[1,1] = poStall0256[0]
tStall[2,1] = poStall02530[0]

tStall[0,2] = poStall0503[0]
tStall[1,2] = poStall0506[0]
tStall[2,2] = poStall05030[0]

tStall[0,3] = poStall0753[0]
tStall[1,3] = poStall0756[0]
tStall[2,3] = poStall07530[0]

tStall[0,4] = poStall13[0]
tStall[1,4] = poStall16[0]
tStall[2,4] = poStall130[0]

tStall[0,5] = poStall1253[0]
tStall[1,5] = poStall1256[0]
tStall[2,5] = poStall12530[0]

tStall[0,6] = poStall1503[0]
tStall[1,6] = poStall1506[0]
tStall[2,6] = poStall15030[0]

tStall[0,7] = poStall1753[0]
tStall[1,7] = poStall1756[0]
tStall[2,7] = poStall17530[0]

tStall[0,8] = poStall23[0]
tStall[1,8] = poStall26[0]
tStall[2,8] = poStall230[0]

tStall[0,9] = poStall2253[0]
tStall[1,9] = poStall2256[0]
tStall[2,9] = poStall22530[0]

tStall[0,10] = poStall2503[0]
tStall[1,10] = poStall2506[0]
tStall[2,10] = poStall25030[0]

tStall[0,11] = poStall2753[0]
tStall[1,11] = poStall2756[0]
tStall[2,11] = poStall27530[0]

tStall[0,12] = poStall33[0]
tStall[1,12] = poStall36[0]
tStall[2,12] = poStall330[0]

tStall[0,13] = poStall3253[0]
tStall[1,13] = poStall3256[0]
tStall[2,13] = poStall32530[0]

tStall[0,14] = poStall3503[0]
tStall[1,14] = poStall3506[0]
tStall[2,14] = poStall35030[0]

tStall[0,15] = poStall3753[0]
tStall[1,15] = poStall3756[0]
tStall[2,15] = poStall37530[0]

tStall[0,16] = poStall43[0]
tStall[1,16] = poStall46[0]
tStall[2,16] = poStall430[0]

tStall[0,17] = poStall4253[0]
tStall[1,17] = poStall4256[0]
tStall[2,17] = poStall42530[0]

tStall[0,18] = poStall4503[0]
tStall[1,18] = poStall4506[0]
tStall[2,18] = poStall45030[0]

tStall[0,19] = poStall4753[0]
tStall[1,19] = poStall4756[0]
tStall[2,19] = poStall47530[0]

tStall[0,20] = poStall53[0]
tStall[1,20] = poStall56[0]
tStall[2,20] = poStall530[0]


tStall

#%%

m=21

J=np.arange(0,5.25,0.25)


plt.rcParams["figure.figsize"] =[8.0,6.0]
plt.title('Resurrection relax times',fontsize=16)
plt.xlabel('J  ($k_BT$)', fontsize=16)
plt.ylabel(r'$\tau_{rel}$ (s)',fontsize=16)
plt.grid()


plt.plot(J,tStall[0,:],marker='.',linestyle="",markersize=7,label=r'$\alpha=3.06$')
plt.plot(J,tStall[1,:],marker='.',linestyle="",markersize=7,label=r'$\alpha=6.15$')
plt.plot(J,tStall[2,:],marker='.',linestyle="",markersize=7,label=r'$\alpha=30.44$')

plt.legend(loc='upper right',prop={'size': 16})

#%%

Dt03_N8= (tStall[0,:]-tRes[0,:])/(tStall[0,:]+tRes[0,:])
Dt06_N8= (tStall[1,:]-tRes[1,:])/(tStall[1,:]+tRes[1,:])
Dt30_N8= (tStall[2,:]-tRes[2,:])/(tStall[2,:]+tRes[2,:])

plt.rcParams["figure.figsize"] =[8.0,6.0]
plt.title('',fontsize=16)
plt.xlabel('J  ($k_BT$)', fontsize=16)
plt.ylabel(r'$\Delta$',fontsize=16)
plt.grid()

plt.plot(J,Dt03_N8,marker='.',linestyle="",markersize=7,label=r'$\alpha=3.06$')
plt.plot(J,Dt06_N8,marker='.',linestyle="",markersize=7,label=r'$\alpha=6.15$')
plt.plot(J,Dt30_N8,marker='.',linestyle="",markersize=7,label=r'$\alpha=30.44$')

plt.legend(loc='upper right',prop={'size': 16}, framealpha=0.3)

#%%

plt.rcParams["figure.figsize"] =[8.0,6.0]
plt.title('',fontsize=16)
plt.xlabel('J  ($k_BT$)', fontsize=16)
plt.ylabel(r'$\Delta$',fontsize=16)
plt.grid()

plt.plot(J,Dt03_N8,marker='^',linestyle="",markersize=5, label=r'$\alpha=3.06$',color='purple')
plt.plot(J,Dt06_N8,marker='^',linestyle="",markersize=5,label=r'$\alpha=6.15$',color='orange')
plt.plot(J,Dt30_N8,marker='^',linestyle="",markersize=5,label=r'$\alpha=30.44$',color='royalblue')

plt.plot(J,Dt03_N4,marker='o',linestyle="",markersize=5,label=r'$\alpha=3.06$',color='purple')
plt.plot(J,Dt06_N4,marker='o',linestyle="",markersize=5,label=r'$\alpha=6.15$',color='orange')
plt.plot(J,Dt30_N4,marker='o',linestyle="",markersize=5,label=r'$\alpha=30.44$',color='royalblue')


a306= mlines.Line2D([],[],color="purple",marker="s",linestyle="",label=r'$\alpha=3.06$')
a615= mlines.Line2D([],[],color="orange",marker="s",linestyle="",label=r'$\alpha=6.15$')
a3044= mlines.Line2D([],[],color="royalblue",marker="s",linestyle="",label=r'$\alpha=30.44$')

n4= mlines.Line2D([],[],color="black",marker="o",linestyle="",label=r'$<N>=4$')
n8= mlines.Line2D([],[],color="black",marker="^",linestyle="",label=r'$<N>=8$')

plt.legend(loc='upper right',prop={'size': 16}, framealpha=0.3,bbox_to_anchor=(1.35,1),handles=[a306,
           a615,a3044,n4,n8])
    
