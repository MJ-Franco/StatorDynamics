#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 10:58:11 2022

@author: MJ Franco
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
plt.rcParams["figure.figsize"] = [10.0,8.0]

"""Data reading"""

# resurrection data
FileR=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/NewBergKMC/Traces_Res_KMC_NBerg.dat']
DataR=[]

for data_file in FileR:
    DataR.append(np.loadtxt(data_file))
    
Res=DataR[0]

# release data
FileRel=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/NewBergKMC/Traces_Stall_KMC_NBerg.dat']
DataRel=[]

for data_file in FileRel:
    DataRel.append(np.loadtxt(data_file))
    
Rel=DataRel[0]

"""Useful parameters"""

ntime = 5000

"""Plot traces"""

plt.title('Berg model (parameters from article)',fontsize=16)
plt.xlabel('Time / s',fontsize = 14)
plt.ylabel('Stator number',fontsize= 14)
plt.grid()
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)

plt.plot(Res[0:ntime,0],Res[0:ntime,1])
plt.plot(Rel[0:ntime,0],Rel[0:ntime,1])


""""Relaxation times"""

