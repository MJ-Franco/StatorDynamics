# -*- coding: utf-8 -*-
"""
Created on Fri Sep 17 11:03:06 2021

@author: mjfra
"""

#Libraries

import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import curve_fit
import matplotlib.lines as mlines

#%% Exponential fit

def N_exp(t,tau,nss):
    return nss + (n0-nss)*np.exp(-t/tau)

#%% Ruben's Data load

with open('/media/mfranco/Elements/Simulations/Data/Speed-rate/friction300.pkl', 'rb') as f:
    data300 = pickle.load(f)
    
with open('/media/mfranco/Elements/Simulations/Data/Speed-rate/friction500.pkl', 'rb') as f:
    data500 = pickle.load(f)

with open('/media/mfranco/Elements/Simulations/Data/Speed-rate/friction1300.pkl', 'rb') as f:
    data1300 = pickle.load(f)   
    

#%% Organisation of Ruben's data

res300_dat = data300['resurrection']
rel300_dat = data300['release']

res500_dat = data500['resurrection']
rel500_dat = data500['release']

res1300_dat = data1300['resurrection']
rel1300_dat = data1300['release']

#%% Data from KMC simulations

files1=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg.dat']
data1=[]

files2=['/media/mfranco/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg.dat']
data2=[]

for data_file in files1:
    data1.append(np.loadtxt(data_file))

for data_file in files2:
    data2.append(np.loadtxt(data_file))


Res=data1[0]
Rel=data2[0]

l=10000
nmax=13



#%% Mean field 

def model(N,t):
    a_g = 30.44
    k0 = 0.00575
    mu = 1.19
    Nm = 13
    dNdt = k0*(1. - np.exp(-a_g/N))*(Nm - N*(np.exp(-mu) + 1.))
    
    return dNdt

t = np.linspace(0,800,10000)

#%% Solution of mean field

N0_res = 0
N_res = odeint(model,N0_res,t)


N0_rel = 9.965
N_rel = odeint(model,N0_rel,t)


#%% Fit resurrection

n0 = 0

p0=[100,4]


poRes_R, pcRes_R = curve_fit(N_exp,res1300_dat['t(s)'],res300_dat['N'])
poRes_MF, pcRes_MF = curve_fit(N_exp,t,N_res[:,0])
poRes_sim, pcRes_sim = curve_fit(N_exp,Res[0:l,0],Res[0:l,1],p0)

#%% Fit release

n0 = 9.25


p0=[100,4]


poRel_R, pcRel_R = curve_fit(N_exp,rel1300_dat['t(s)'],rel300_dat['N'])
poRel_MF, pcRel_MF = curve_fit(N_exp,t,N_rel[:,0])
poRel_sim, pcRel_sim = curve_fit(N_exp,Rel[0:l,0],Rel[0:l,1],p0)

#%% Comparison Ruben vs mean field

plt.rcParams["figure.figsize"] =[8.0,6.0]
plt.title('1300 nm bead',fontsize=14)
plt.xlabel('Time (s)', fontsize=16)
plt.ylabel('N',fontsize=16)
plt.grid()

plt.plot(res1300_dat['t(s)'],res1300_dat['N'],color='blue',label=r'$\tau$={:1.2f}'.format(poRes_R[0]))
plt.plot(rel1300_dat['t(s)'],rel1300_dat['N'],color='orange',label=r'$\tau$={:1.2f}'.format(poRel_R[0]))

plt.plot(t,N_res,color='blue',label=r'$\tau$={:1.2f}'.format(poRes_MF[0]),linestyle=':')
plt.plot(t,N_rel,color='orange',label=r'$\tau$={:1.2f}'.format(poRel_MF[0]),linestyle=':')

plt.legend()

black_line1 = mlines.Line2D([],[],linestyle='', label=r'$\Delta \tau_{R}=0.03$')
black_line2 = mlines.Line2D([],[],linestyle='', label=r'$\Delta \tau_{MF}=0.02$')
Res_R = mlines.Line2D([],[],linestyle='-', color='blue',label=r'$\tau$={:1.1f}'.format(poRes_R[0]))
Rel_R = mlines.Line2D([],[],linestyle='-', color='orange',label=r'$\tau$={:1.1f}'.format(poRel_R[0]))
Res_MF = mlines.Line2D([],[],linestyle=':', color='blue',label=r'$\tau$={:1.1f}'.format(poRes_MF[0]))
Rel_MF = mlines.Line2D([],[],linestyle=':', color='orange',label=r'$\tau$={:1.1f}'.format(poRel_MF[0]))
title_1 = mlines.Line2D([],[],linestyle='', label='Master equation')
title_2 = mlines.Line2D([],[],linestyle='', label='Mean field')


plt.legend(loc='lower right',handles=[title_1,Res_R,Rel_R,black_line1,title_2,Res_MF,Rel_MF,black_line2])

#%% Comparison Ruben vs sim

plt.rcParams["figure.figsize"] =[8.0,6.0]
plt.title('1300 nm bead',fontsize=16)
plt.xlabel('Time (s)', fontsize=16)
plt.ylabel('N',fontsize=16)
plt.grid()

plt.plot(res1300_dat['t(s)'],res1300_dat['N'],color='black',label=r'$\tau$={:1.2f}'.format(poRes_R[0]))
plt.plot(rel1300_dat['t(s)'],rel1300_dat['N'],color='black',label=r'$\tau$={:1.2f}'.format(poRel_R[0]))

plt.plot(Res[0:l,0],Res[0:l,1],linestyle='--',color='red',label=r'$\tau$={:1.2f}'.format(poRes_sim[0]))
plt.plot(Rel[0:l,0],Rel[0:l,1],linestyle='--',color='red',label=r'$\tau$={:1.2f}'.format(poRel_sim[0]))

plt.legend()

black_line1 = mlines.Line2D([],[],linestyle='', label=r'$\Delta \tau_{R}=0.28$')
black_line2 = mlines.Line2D([],[],linestyle='', label=r'$\Delta \tau_{MC}=0.30$')
Res_R = mlines.Line2D([],[],linestyle='-', color='blue',label=r'$\tau$={:1.1f}'.format(poRes_R[0]))
Rel_R = mlines.Line2D([],[],linestyle='-', color='orange',label=r'$\tau$={:1.1f}'.format(poRel_R[0]))
Res_MF = mlines.Line2D([],[],linestyle='--', color='blue',label=r'$\tau$={:1.1f}'.format(poRes_sim[0]))
Rel_MF = mlines.Line2D([],[],linestyle='--', color='orange',label=r'$\tau$={:1.1f}'.format(poRel_sim[0]))
title_1 = mlines.Line2D([],[],linestyle='-', color='black', label='Master equation')
title_2 = mlines.Line2D([],[],linestyle='--', color='red', label='KMC')



#plt.legend(loc='lower right',handles=[title_1,Res_R,Rel_R,black_line1,title_2,Res_MF,Rel_MF,black_line2])

plt.legend(loc='lower right',handles=[title_1,title_2],prop={'size': 16})

#%%


plt.plot(Res[0:l,0],Res[0:l,1],linestyle='--',color='red',label=r'$\tau$={:1.2f}'.format(poRes_sim[0]))
