#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 14:20:23 2021

@author: mariajose
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.patches as mpatches
import matplotlib.lines as mlines


#%%

nmax=12
l=80000
#Function used for fits
def N(x,tau,nss):
    return nss + (n0-nss)*np.exp(-x/tau) 

#%% Files for resurrection

files1=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Res_J=0.00  .dat'] #All resurrection simulations
data1=[]
files11=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Res_J=0.25  .dat'] #All resurrection simulations
data11=[]
files12=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Res_J=0.50  .dat'] #All resurrection simulations
data12=[]
files13=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Res_J=0.75  .dat'] #All resurrection simulations
data13=[]

for data_file in files1:
    data1.append(np.loadtxt(data_file))
    
for data_file in files11:
    data11.append(np.loadtxt(data_file))
    
for data_file in files12:
    data12.append(np.loadtxt(data_file))

for data_file in files13:
    data13.append(np.loadtxt(data_file))
    
Res_J_0=data1[0]
Res_J_02=data11[0]
Res_J_05=data12[0]
Res_J_07=data12[0]

files2=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Res_J=1.00  .dat'] #Average of resurrection simulations
data2=[]
files21=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Res_J=1.25  .dat'] #Average of resurrection simulations
data21=[]
files22=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Res_J=1.50  .dat'] #Average of resurrection simulations
data22=[]
files23=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Res_J=1.75  .dat'] #Average of resurrection simulations
data23=[]

for data_file in files2:
    data2.append(np.loadtxt(data_file))
    
for data_file in files21:
    data21.append(np.loadtxt(data_file))
    
for data_file in files22:
    data22.append(np.loadtxt(data_file))

for data_file in files23:
    data23.append(np.loadtxt(data_file))
    
Res_J_1=data2[0]
Res_J_12=data21[0]
Res_J_15=data22[0]
Res_J_17=data23[0]

files3=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Res_J=2.00  .dat'] #Reservoir in each simulation
data3=[]
files31=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Res_J=2.25  .dat'] #Reservoir in each simulation
data31=[]
files32=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Res_J=2.50  .dat'] #Reservoir in each simulation
data32=[]
files33=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Res_J=2.75  .dat'] #Reservoir in each simulation
data33=[]

for data_file in files3:
    data3.append(np.loadtxt(data_file)) 
    
for data_file in files31:
    data31.append(np.loadtxt(data_file)) 
    
for data_file in files32:
    data32.append(np.loadtxt(data_file))
    
for data_file in files33:
    data33.append(np.loadtxt(data_file)) 
    
Res_J_2=data3[0]
Res_J_22=data31[0]
Res_J_25=data32[0]
Res_J_27=data33[0]

files4=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Res_J=3.00  .dat'] #Average of reservoir in resurrection simulations
data4=[]
files41=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Res_J=3.25  .dat'] #Average of reservoir in resurrection simulations
data41=[]
files42=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Res_J=3.50  .dat'] #Average of reservoir in resurrection simulations
data42=[]
files43=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Res_J=3.75  .dat'] #Average of reservoir in resurrection simulations
data43=[]

for data_file in files4:
    data4.append(np.loadtxt(data_file))
    
for data_file in files41:
    data41.append(np.loadtxt(data_file))
    
for data_file in files42:
    data42.append(np.loadtxt(data_file))
    
for data_file in files43:
    data43.append(np.loadtxt(data_file))
    
Res_J_3=data4[0]
Res_J_32=data41[0]
Res_J_35=data42[0]
Res_J_37=data43[0]

files5=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Res_J=4.00  .dat'] #Average of reservoir in resurrection simulations
data5=[]

for data_file in files5:
    data5.append(np.loadtxt(data_file))

Res_J_4=data5[0]


#%% Files for stall

files6=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Stall_J=0.00  .dat'] #All resurrection simulations
data6=[]
files61=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Stall_J=0.25  .dat'] #All resurrection simulations
data61=[]
files62=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Stall_J=0.50  .dat'] #All resurrection simulations
data62=[]
files63=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Stall_J=0.75  .dat'] #All resurrection simulations
data63=[]

for data_file in files6:
    data6.append(np.loadtxt(data_file))
    
for data_file in files61:
    data61.append(np.loadtxt(data_file))
    
for data_file in files62:
    data62.append(np.loadtxt(data_file))
    
for data_file in files63:
    data63.append(np.loadtxt(data_file))
    
    
Stall_J_0=data6[0]
Stall_J_02=data61[0]
Stall_J_05=data62[0]
Stall_J_07=data63[0]

files7=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Stall_J=1.00  .dat'] #Average of resurrection simulations
data7=[]
files71=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Stall_J=1.25  .dat'] #Average of resurrection simulations
data71=[]
files72=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Stall_J=1.50  .dat'] #Average of resurrection simulations
data72=[]
files73=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Stall_J=1.75  .dat'] #Average of resurrection simulations
data73=[]

for data_file in files7:
    data7.append(np.loadtxt(data_file))
    
for data_file in files71:
    data71.append(np.loadtxt(data_file))
    
for data_file in files72:
    data72.append(np.loadtxt(data_file))

for data_file in files73:
    data73.append(np.loadtxt(data_file))
    
    
Stall_J_1=data7[0]
Stall_J_12=data71[0]
Stall_J_15=data72[0]
Stall_J_17=data73[0]

files8=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Stall_J=2.00  .dat'] #Reservoir in each simulation
data8=[]
files81=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Stall_J=2.25  .dat'] #Reservoir in each simulation
data81=[]
files82=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Stall_J=2.50  .dat'] #Reservoir in each simulation
data82=[]
files83=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Stall_J=2.75  .dat'] #Reservoir in each simulation
data83=[]

for data_file in files8:
    data8.append(np.loadtxt(data_file))
    
for data_file in files81:
    data81.append(np.loadtxt(data_file))
    
for data_file in files82:
    data82.append(np.loadtxt(data_file))
    
for data_file in files83:
    data83.append(np.loadtxt(data_file))
    
    
Stall_J_2=data8[0]
Stall_J_22=data81[0]
Stall_J_25=data82[0]
Stall_J_27=data83[0]

files9=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Stall_J=3.00  .dat'] #Average of reservoir in resurrection simulations
data9=[]
files91=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Stall_J=3.25  .dat'] #Average of reservoir in resurrection simulations
data91=[]
files92=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Stall_J=3.50  .dat'] #Average of reservoir in resurrection simulations
data92=[]
files93=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Stall_J=3.75  .dat'] #Average of reservoir in resurrection simulations
data93=[]

for data_file in files9:
    data9.append(np.loadtxt(data_file))
    
for data_file in files91:
    data91.append(np.loadtxt(data_file))
    
for data_file in files92:
    data92.append(np.loadtxt(data_file))
    
for data_file in files93:
    data93.append(np.loadtxt(data_file))
    
Stall_J_3=data9[0]
Stall_J_32=data91[0]
Stall_J_35=data92[0]
Stall_J_37=data93[0]

files10=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Stall_J=4.00  .dat'] #Average of reservoir in resurrection simulations
data10=[]
    
for data_file in files10:
    data10.append(np.loadtxt(data_file))

Stall_J_4=data10[0]

#%%

n0=0

poRes0, pcRes0 = curve_fit(N,Res_J_0[0:l,0],Res_J_0[0:l,1]*nmax)
tau_Res0=poRes0[0]
err_Res0=np.sqrt(pcRes0[0,0])

poRes02, pcRes02 = curve_fit(N,Res_J_02[0:l,0],Res_J_02[0:l,1]*nmax)
tau_Res02=poRes02[0]
err_Res02=np.sqrt(pcRes02[0,0])

poRes05, pcRes05 = curve_fit(N,Res_J_05[0:l,0],Res_J_05[0:l,1]*nmax)
tau_Res05=poRes05[0]
err_Res05=np.sqrt(pcRes05[0,0])

poRes07, pcRes07 = curve_fit(N,Res_J_07[0:l,0],Res_J_07[0:l,1]*nmax)
tau_Res07=poRes07[0]
err_Res07=np.sqrt(pcRes07[0,0])



poRes1, pcRes1 = curve_fit(N,Res_J_1[0:l,0],Res_J_1[0:l,1]*nmax)
tau_Res1=poRes1[0]
err_Res1=np.sqrt(pcRes1[0,0])

poRes12, pcRes12 = curve_fit(N,Res_J_12[0:l,0],Res_J_12[0:l,1]*nmax)
tau_Res12=poRes12[0]
err_Res12=np.sqrt(pcRes12[0,0])

poRes15, pcRes15 = curve_fit(N,Res_J_15[0:l,0],Res_J_15[0:l,1]*nmax)
tau_Res15=poRes15[0]
err_Res15=np.sqrt(pcRes15[0,0])

poRes17, pcRes17 = curve_fit(N,Res_J_17[0:l,0],Res_J_17[0:l,1]*nmax)
tau_Res17=poRes17[0]
err_Res17=np.sqrt(pcRes17[0,0])



poRes2, pcRes2 = curve_fit(N,Res_J_2[0:l,0],Res_J_2[0:l,1]*nmax)
tau_Res2=poRes2[0]
err_Res2=np.sqrt(pcRes2[0,0])

poRes22, pcRes22 = curve_fit(N,Res_J_22[0:l,0],Res_J_22[0:l,1]*nmax)
tau_Res22=poRes22[0]
err_Res22=np.sqrt(pcRes22[0,0])

poRes25, pcRes25 = curve_fit(N,Res_J_25[0:l,0],Res_J_25[0:l,1]*nmax)
tau_Res25=poRes25[0]
err_Res25=np.sqrt(pcRes25[0,0])

poRes27, pcRes27 = curve_fit(N,Res_J_27[0:l,0],Res_J_27[0:l,1]*nmax)
tau_Res27=poRes27[0]
err_Res27=np.sqrt(pcRes27[0,0])



poRes3, pcRes3 = curve_fit(N,Res_J_3[0:l,0],Res_J_3[0:l,1]*nmax)
tau_Res3=poRes3[0]
err_Res3=np.sqrt(pcRes3[0,0])

poRes32, pcRes32 = curve_fit(N,Res_J_32[0:l,0],Res_J_32[0:l,1]*nmax)
tau_Res32=poRes32[0]
err_Res32=np.sqrt(pcRes32[0,0])

poRes35, pcRes35 = curve_fit(N,Res_J_35[0:l,0],Res_J_35[0:l,1]*nmax)
tau_Res35=poRes35[0]
err_Res35=np.sqrt(pcRes35[0,0])

poRes37, pcRes37 = curve_fit(N,Res_J_37[0:l,0],Res_J_37[0:l,1]*nmax)
tau_Res37=poRes37[0]
err_Res37=np.sqrt(pcRes37[0,0])



poRes4, pcRes4 = curve_fit(N,Res_J_4[0:l,0],Res_J_4[0:l,1]*nmax)
tau_Res4=poRes4[0]
err_Res4=np.sqrt(pcRes4[0,0])

#%%

plt.title('Resurrection without depletion ($N_{ss}=8$)',fontsize=17)
plt.xlabel('Simulation time',fontsize=14)
plt.ylabel('Stator number',fontsize=14)
plt.rcParams["figure.figsize"] = [8.0,6.0]

plt.plot(Res_J_0[0:l,0],Res_J_0[0:l,1]*nmax,label=r'J=0 $k_BT$, $\tau$={:1.2f}'.format(tau_Res0))
plt.plot(Res_J_1[0:l,0],Res_J_1[0:l,1]*nmax,label=r'J=1 $k_BT$, $\tau$={:1.2f}'.format(tau_Res1))
plt.plot(Res_J_2[0:l,0],Res_J_2[0:l,1]*nmax,label=r'J=2 $k_BT$, $\tau$={:1.2f}'.format(tau_Res2))
plt.plot(Res_J_3[0:l,0],Res_J_3[0:l,1]*nmax,label=r'J=3 $k_BT$, $\tau$={:1.2f}'.format(tau_Res3))
plt.plot(Res_J_4[0:l,0],Res_J_4[0:l,1]*nmax,label=r'J=4 $k_BT$, $\tau$={:1.2f}'.format(tau_Res4))      

plt.legend(loc='lower right', prop={'size': 15})

#%% Fit for stall without depletion

n0=12

poStall0, pcStall0 = curve_fit(N,Stall_J_0[0:l,0],Stall_J_0[0:l,1]*nmax)
tau_Stall0=poStall0[0]
err_Stall0=np.sqrt(pcStall0[0,0])

poStall02, pcStall02= curve_fit(N,Stall_J_02[0:l,0],Stall_J_02[0:l,1]*nmax)
tau_Stall02=poStall02[0]
err_Stall02=np.sqrt(pcStall02[0,0])

poStall05, pcStall05 = curve_fit(N,Stall_J_05[0:l,0],Stall_J_05[0:l,1]*nmax)
tau_Stall05=poStall05[0]
err_Stall05=np.sqrt(pcStall05[0,0])

poStall07, pcStall07 = curve_fit(N,Stall_J_07[0:l,0],Stall_J_07[0:l,1]*nmax)
tau_Stall07=poStall07[0]
err_Stall07=np.sqrt(pcStall07[0,0])



poStall1, pcStall1 = curve_fit(N,Stall_J_1[0:l,0],Stall_J_1[0:l,1]*nmax)
tau_Stall1=poStall1[0]
err_Stall1=np.sqrt(pcStall1[0,0])

poStall12, pcStall12 = curve_fit(N,Stall_J_12[0:l,0],Stall_J_12[0:l,1]*nmax)
tau_Stall12=poStall12[0]
err_Stall12=np.sqrt(pcStall12[0,0])

poStall15, pcStall15 = curve_fit(N,Stall_J_15[0:l,0],Stall_J_15[0:l,1]*nmax)
tau_Stall15=poStall15[0]
err_Stall15=np.sqrt(pcStall15[0,0])

poStall17, pcStall17 = curve_fit(N,Stall_J_17[0:l,0],Stall_J_17[0:l,1]*nmax)
tau_Stall17=poStall17[0]
err_Stall17=np.sqrt(pcStall17[0,0])



poStall2, pcStall2 = curve_fit(N,Stall_J_2[0:l,0],Stall_J_2[0:l,1]*nmax)
tau_Stall2=poStall2[0]
err_Stall2=np.sqrt(pcStall2[0,0])

poStall22, pcStall22 = curve_fit(N,Stall_J_22[0:l,0],Stall_J_22[0:l,1]*nmax)
tau_Stall22=poStall22[0]
err_Stall22=np.sqrt(pcStall22[0,0])

poStall25, pcStall25 = curve_fit(N,Stall_J_25[0:l,0],Stall_J_25[0:l,1]*nmax)
tau_Stall25=poStall25[0]
err_Stall25=np.sqrt(pcStall25[0,0])

poStall27, pcStall27 = curve_fit(N,Stall_J_27[0:l,0],Stall_J_27[0:l,1]*nmax)
tau_Stall27=poStall27[0]
err_Stall27=np.sqrt(pcStall27[0,0])



poStall3, pcStall3 = curve_fit(N,Stall_J_3[0:l,0],Stall_J_3[0:l,1]*nmax)
tau_Stall3=poStall3[0]
err_Stall3=np.sqrt(pcStall3[0,0])

poStall32, pcStall32 = curve_fit(N,Stall_J_32[0:l,0],Stall_J_32[0:l,1]*nmax)
tau_Stall32=poStall32[0]
err_Stall32=np.sqrt(pcStall32[0,0])

poStall35, pcStall35 = curve_fit(N,Stall_J_35[0:l,0],Stall_J_35[0:l,1]*nmax)
tau_Stall35=poStall35[0]
err_Stall35=np.sqrt(pcStall35[0,0])

poStall37, pcStall37 = curve_fit(N,Stall_J_37[0:l,0],Stall_J_37[0:l,1]*nmax)
tau_Stall37=poStall37[0]
err_Stall37=np.sqrt(pcStall37[0,0])



poStall4, pcStall4 = curve_fit(N,Stall_J_4[0:l,0],Stall_J_4[0:l,1]*nmax)
tau_Stall4=poStall4[0]
err_Stall4=np.sqrt(pcStall4[0,0])

#%%

plt.title('Stall without depletion ($N_{ss}=8$)',fontsize=17)
plt.xlabel('Simulation time',fontsize=14)
plt.ylabel('Stator number',fontsize=14)
plt.rcParams["figure.figsize"] = [8.0,6.0]

plt.plot(Stall_J_0[0:l,0],Stall_J_0[0:l,1]*nmax,label=r'J=0 $k_BT$, $\tau$={:1.2f}'.format(tau_Stall0))
plt.plot(Stall_J_1[0:l,0],Stall_J_1[0:l,1]*nmax,label=r'J=1 $k_BT$, $\tau$={:1.2f}'.format(tau_Stall1))
plt.plot(Stall_J_2[0:l,0],Stall_J_2[0:l,1]*nmax,label=r'J=2 $k_BT$, $\tau$={:1.2f}'.format(tau_Stall2))
plt.plot(Stall_J_3[0:l,0],Stall_J_3[0:l,1]*nmax,label=r'J=3 $k_BT$, $\tau$={:1.2f}'.format(tau_Stall3))
plt.plot(Stall_J_4[0:l,0],Stall_J_4[0:l,1]*nmax,label=r'J=4 $k_BT$, $\tau$={:1.2f}'.format(tau_Stall4))

#plt.hlines(8.05,0,200, linestyles='dashed',alpha=0.4)    

plt.legend(loc='right', prop={'size': 15})

#%%No depletion

J=np.arange(0,4.25,0.25)


Relax_t_Res = np.array([tau_Res0,tau_Res02,tau_Res05,tau_Res07,tau_Res1,tau_Res12,tau_Res15,tau_Res17,tau_Res2,tau_Res22,
                        tau_Res25,tau_Res27,tau_Res3,tau_Res32,tau_Res35,tau_Res37,tau_Res4])

avg_Relax_t_Res = sum(Relax_t_Res)/len(Relax_t_Res)

#std_Relax_t_Res = np.sqrt(sum((Relax_t_Res-avg_Relax_t_Res)**2)/len(Relax_t_Res))

std_Relax_t_Res=np.array([err_Res0,err_Res02,err_Res05,err_Res07,err_Res1,err_Res12,err_Res15,err_Res17,err_Res2,err_Res22,
                          err_Res25,err_Res27,err_Res3,err_Res32,err_Res35,err_Res37,err_Res4])

Relax_t_Stall = np.array([tau_Stall0,tau_Stall02,tau_Stall05,tau_Stall07,tau_Stall1,tau_Stall12,tau_Stall15,tau_Stall17,
                          tau_Stall2,tau_Stall22,tau_Stall25,tau_Stall27,tau_Stall3,tau_Stall32,tau_Stall35,tau_Stall37,tau_Stall4])

avg_Relax_t_Stall = sum(Relax_t_Stall)/len(Relax_t_Stall)

#std_Relax_t_Stall = np.sqrt(sum((Relax_t_Stall - avg_Relax_t_Stall)**2)/len(Relax_t_Stall))

std_Relax_t_Stall = np.array([err_Stall0,err_Stall02,err_Stall05,err_Stall07,err_Stall1,err_Stall12,err_Stall15,err_Stall17,
                              err_Stall2,err_Stall22,err_Stall25,err_Stall27,err_Stall3,err_Stall32,err_Stall35,err_Stall37,err_Stall4])

#%%

#np.savetxt('/home/mariajose/Escritorio/Simulations/Without depletion/Relaxation times/ResTime_6.dat',(np.c_[J,Relax_t_Res,std_Relax_t_Res]))
#np.savetxt('/home/mariajose/Escritorio/Simulations/Without depletion/Relaxation times/StallTime_6.dat',(np.c_[J,Relax_t_Stall,std_Relax_t_Stall]))

#%% Metropolis without depletion

M_tau_Res_03 = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without depletion/Relaxation times/ResTime_4.dat')
M_tau_Stall_03 = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without depletion/Relaxation times/StallTime_4.dat')

M_tau_Res_05 = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without depletion/Relaxation times/ResTime_6.dat')
M_tau_Stall_05 = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without depletion/Relaxation times/StallTime_6.dat')

M_tau_Res_07 = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without depletion/Relaxation times/ResTime_8.dat')
M_tau_Stall_07 = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without depletion/Relaxation times/StallTime_8.dat')


#%% Metropolis with depletion



#%% Glauber without depletion

G_tau_Res_05 = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/RelaxTimes/ResTime_J=-mu.dat')
G_tau_Stall_05 = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/RelaxTimes/StallTime_J=-mu.dat')

G_tau_Res_03 = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/RelaxTimes/ResTime_4.dat')
G_tau_Stall_03 = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/RelaxTimes/StallTime_4.dat')

G_tau_Res_07 = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/RelaxTimes/ResTime_8.dat')
G_tau_Stall_07 = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/RelaxTimes/StallTime_8.dat')

#%% Glauber with depletion

G_wd_tau_Res_03 = np.loadtxt('/home/mariajose/Escritorio/Simulations/With depletion/Glauber/RelaxTimes/ResTime_4.dat')
G_wd_tau_Stall_03 = np.loadtxt('/home/mariajose/Escritorio/Simulations/With depletion/Glauber/RelaxTimes/StallTime_4.dat')

G_wd_tau_Res_05 = np.loadtxt('/home/mariajose/Escritorio/Simulations/With depletion/Glauber/RelaxTimes/ResTime_6.dat')
G_wd_tau_Stall_05 = np.loadtxt('/home/mariajose/Escritorio/Simulations/With depletion/Glauber/RelaxTimes/StallTime_6.dat')

G_wd_tau_Res_07 = np.loadtxt('/home/mariajose/Escritorio/Simulations/With depletion/Glauber/RelaxTimes/ResTime_8.dat')
G_wd_tau_Stall_07 = np.loadtxt('/home/mariajose/Escritorio/Simulations/With depletion/Glauber/RelaxTimes/StallTime_8.dat')

#%% Fred witohut depletion

F_tau_Res_03 = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without depletion/Fred Dynamics/Relaxation times/ResTime_4.dat')
F_tau_Stall_03 = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without depletion/Fred Dynamics/Relaxation times/StallTime_4.dat')

F_tau_Res_05 = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without depletion/Fred Dynamics/Relaxation times/ResTime_6.dat')
F_tau_Stall_05 = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without depletion/Fred Dynamics/Relaxation times/StallTime_6.dat')

F_tau_Res_07 = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without depletion/Fred Dynamics/Relaxation times/ResTime_8.dat')
F_tau_Stall_07 = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without depletion/Fred Dynamics/Relaxation times/StallTime_8.dat')

#%% Fred with depletion


#%%

ax = plt.subplot(111)
#ax.set_xscale("log", nonposx='clip')
#ax.set_yscale("log", nonposy='clip')

ax.set_title('Comparison Metropolis/Glauber/Fred ($<\Phi> = 1/2$)',fontsize=17)
ax.set_xlabel(r'J / $k_BT$',fontsize=14)
ax.set_ylabel('Effective relaxation time',fontsize=14)
#plt.rcParams["figure.figsize"] = [8.0,6.0]

ax.errorbar(J,M_tau_Res_05[:,1],M_tau_Res_05[:,2],capsize=1,fmt='s',label='Resurrection',color='blue')
ax.errorbar(J,M_tau_Stall_05[:,1],M_tau_Stall_05[:,2],capsize=1,fmt='^',label='Stall',color='blue')
ax.errorbar(J,G_tau_Res_05[:,1],G_tau_Res_05[:,2],capsize=1,fmt='s',label='Resurrection w depletion',color='orange')
ax.errorbar(J,G_tau_Stall_05[:,1],G_tau_Stall_05[:,2],capsize=1,fmt='^',label='Stall',color='orange')
ax.errorbar(J,F_tau_Res_05[:,1],F_tau_Res_05[:,2],capsize=1,fmt='s',label='Resurrection w depletion',color='green')
ax.errorbar(J,F_tau_Stall_05[:,1],F_tau_Stall_05[:,2],capsize=1,fmt='^',label='Stall',color='green')

black_line1 = mlines.Line2D([],[],linestyle='', color='black', marker='s', label='Resurrection')
black_line2 = mlines.Line2D([],[],linestyle='', color='black', marker='^', label='Stall')
#black_line3 = mlines.Line2D([],[],linestyle='', color='black', marker='.', label='Fred')
blue_patch = mpatches.Patch(color='blue', label='Metropolis')
orange_patch = mpatches.Patch(color='orange', label='Glauber')
green_patch = mpatches.Patch(color='green', label='Fred')
plt.legend(loc='upper left', prop={'size': 15},handles=[blue_patch,orange_patch,green_patch,black_line1,black_line2])

#%%

ax = plt.subplot(111)
#ax.set_xscale("log", nonposx='clip')
#ax.set_yscale("log", nonposy='clip')

ax.set_title('Fred dynamics',fontsize=17)
ax.set_xlabel(r'J / $k_BT$',fontsize=14)
ax.set_ylabel('Effective relaxation time',fontsize=14)
#plt.rcParams["figure.figsize"] = [8.0,6.0]

ax.errorbar(J,F_tau_Res_03[:,1],F_tau_Res_03[:,2],capsize=1,fmt='s',label='Resurrection',color='blue')
ax.errorbar(J,F_tau_Stall_03[:,1],F_tau_Stall_03[:,2],capsize=1,fmt='^',label='Stall',color='blue')
ax.errorbar(J,F_tau_Res_05[:,1],F_tau_Res_05[:,2],capsize=1,fmt='s',label='Resurrection w depletion',color='orange')
ax.errorbar(J,F_tau_Stall_05[:,1],F_tau_Stall_03[:,2],capsize=1,fmt='^',label='Stall',color='orange')
ax.errorbar(J,F_tau_Res_07[:,1],F_tau_Res_07[:,2],capsize=1,fmt='s',label='Resurrection w depletion',color='green')
ax.errorbar(J,F_tau_Stall_07[:,1],F_tau_Stall_07[:,2],capsize=1,fmt='^',label='Stall',color='green')

black_line1 = mlines.Line2D([],[],linestyle='', color='black', marker='s', label='Resurrection')
black_line2 = mlines.Line2D([],[],linestyle='', color='black', marker='^', label='Stall')
#black_line3 = mlines.Line2D([],[],linestyle='', color='black', marker='.', label='Fred')
blue_patch = mpatches.Patch(color='blue', label='$\Phi=1/3$')
orange_patch = mpatches.Patch(color='orange', label='$\Phi=1/2$')
green_patch = mpatches.Patch(color='green', label='$\Phi=2/3$')
plt.legend(loc='upper left', prop={'size': 15},handles=[blue_patch,orange_patch,green_patch,black_line1,black_line2])