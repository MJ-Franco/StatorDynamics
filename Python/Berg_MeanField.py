# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.lines as mlines

#%% Differential equation

def model(N,t):
    a_g = 3.06
    k0 = 0.00575
    mu = -0.846
    Nm = 13
    dNdt = k0*(1. - np.exp(-a_g/N))*(Nm - N*(np.exp(-mu) + 1.))
    
    return dNdt

t = np.linspace(0,1000,10000)

#%% Resurrection

N0_res = 0

N_res = odeint(model,N0_res,t)

plt.plot(t,N_res)

#%% Release

N0_rel = 6.2

N_rel = odeint(model,N0_rel,t)

plt.plot(t,N_rel)

#%% Exponential fit

def N_exp(t,tau,nss):
    return nss + (n0-nss)*np.exp(-t/tau)


#%% Fit resurrection

n0 = N0_res

poRes, pcRes = curve_fit(N_exp,t,N_res[:,0])

exp_Res = N_exp(t,*poRes)

#%%
 
n0 = N0_rel

poRel, pcRel = curve_fit(N_exp,t,N_rel[:,0])

exp_Rel = N_exp(t,*poRel)

#%%

plt.rcParams["figure.figsize"] = [8.0,6.0]
plt.title(r'$\mu$ = -0.846, $k_0$ = 0.00575, $\alpha_\gamma = 3.06$',fontsize=14)
plt.xlabel('t / s',fontsize=16)
plt.ylabel('N',fontsize=16)
plt.grid()

plt.plot(t,N_res,label=r'$\tau$={:1.2f}'.format(poRes[0]))
plt.plot(t,N_rel,label=r'$\tau$={:1.2f}'.format(poRel[0]))

plt.legend()

#%%

files1=['/media/luke/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Res_KMC_Coop_Berg.dat']
data1=[]

files2=['/media/luke/Elements/Simulations/Without _depletion/KMC/Berg_CoopKMC/Traces_Stall_KMC_Coop_Berg.dat']
data2=[]

for data_file in files1:
    data1.append(np.loadtxt(data_file))

for data_file in files2:
    data2.append(np.loadtxt(data_file))


Res=data1[0]
Rel=data2[0]

l=10000
nmax=13

#%%

n0 = N0_res

p0=[100,4]

poRes_sim, pcRes_sim = curve_fit(N_exp,Res[0:l,1],Res[0:l,2]*nmax,p0)

exp_Res_sim = N_exp(Res[0:l,1],*poRes_sim)

#%%

p0=[100,4]

n0 = N0_rel

poRel_sim, pcRel_sim = curve_fit(N_exp,Rel[0:l,1],Rel[0:l,2]*nmax,p0)

exp_Rel_sim = N_exp(Rel[0:l,1],*poRel_sim)


#%%

plt.title('300 nm',fontsize=14)
plt.xlabel('Time',fontsize=16)
plt.ylabel('N',fontsize=16)
plt.rcParams["figure.figsize"] = [8.0,6.0]
plt.grid()

plt.plot(t,N_res,color='blue',linestyle=':')
plt.plot(t,N_rel,color='orange',linestyle=':')
plt.plot(Res[0:l,0],Res[0:l,1],linestyle='--',color='blue')
plt.plot(Rel[0:l,0],Rel[0:l,1],linestyle='--',color='orange')

plt.legend()

black_line1 = mlines.Line2D([],[],linestyle='', label=r'$\Delta \tau_{MF}=0.31$')
black_line2 = mlines.Line2D([],[],linestyle='', label=r'$\Delta \tau_{MC}=0.23$')
Res_R = mlines.Line2D([],[],linestyle=':', color='blue',label=r'$\tau$={:1.2f}'.format(poRes[0]))
Rel_R = mlines.Line2D([],[],linestyle=':', color='orange',label=r'$\tau$={:1.1f}'.format(poRel[0]))
Res_MF = mlines.Line2D([],[],linestyle='--', color='blue',label=r'$\tau$={:1.2f}'.format(poRes_sim[0]))
Rel_MF = mlines.Line2D([],[],linestyle='--', color='orange',label=r'$\tau$={:1.2f}'.format(poRel_sim[0]))
title_1 = mlines.Line2D([],[],linestyle='', label='Mean field')
title_2 = mlines.Line2D([],[],linestyle='', label='Monte Carlo')


plt.legend(loc='lower right',handles=[title_1,Res_R,Rel_R,black_line1,title_2,Res_MF,Rel_MF,black_line2])
