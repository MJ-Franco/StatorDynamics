#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 15:42:52 2020

@author: mariajose
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit



#%% Useful variables and functions

T=10 #Number of trajectories we want to plot
S=3000 #Total number of steps of the simulation
l=250 #This numbers limits the number of steps we represent
m=S-l #number used for representation
nmax=13

#Function used for fits
def N(x,tau,nss):
    return nss + (n0-nss)*np.exp(-x/tau) 

plt.rcParams["figure.figsize"] = [8.0,6.0]

#%% Files for resurrection

files1=['/media/mfranco/Elements/Simulations/Without _depletion/Glauber/Traces_Res_J=0.00  _mu= 0.47  .dat'] #All resurrection simulations
data1=[]
files2=['/media/mfranco/Elements/Simulations/Without _depletion/Glauber/Traces_Stall_J=0.00  _mu= 0.47  _.dat'] #Average of resurrection simulations
data2=[]
#files3=['/home/mariajose/Escritorio/Simulations/With depletion/Reservoir_Res.dat'] #Reservoir in each simulation
#data3=[]
#files4=['/home/mariajose/Escritorio/Simulations/With depletion/Traces_Res_Avg_Reservoir.dat'] #Average of reservoir in resurrection simulations
#data4=[]

for data_file in files1:
    data1.append(np.loadtxt(data_file))

for data_file in files2:
    data2.append(np.loadtxt(data_file))
    
#for data_file in files3:
#    data3.append(np.loadtxt(data_file)) 
    
#for data_file in files4:
#    data4.append(np.loadtxt(data_file))


Stall=data2[0]
Res=data1[0]

#All_Res_Reservoir=data3[0]
#Res_Reservoir=data4[0]


#%% Files for stall
"""
files5=['/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/GTraces_Stall.dat'] #All stall simulations
data5=[]
files6=['/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/Traces_Stall_J=4.00  .dat'] #Average of stall simulations
data6=[]
#files7=['/home/mariajose/Escritorio/Simulations/With depletion/Reservoir_Stall.dat'] #Reservoir in each stall simulation
#data7=[]
#files8 = ['/home/mariajose/Escritorio/Simulations/With depletion/Traces_Stall_Avg_Reservoir.dat'] #Average of reservoir in stall simulations
#data8=[] 

for data_file in files5:
    data5.append(np.loadtxt(data_file))

for data_file in files6:
    data6.append(np.loadtxt(data_file))

#for data_file in files7:
#    data7.append(np.loadtxt(data_file))
    
#for data_file in files8:
#    data8.append(np.loadtxt(data_file))  

  
All_Stall=data5[0]
Stall=data6[0]

#All_Stall_Reservoir=data7[0]
#Stall_Reservoir=data8[0]

#%% Resurrection representation. Fit

#n0=0


#poRes, pcRes = curve_fit(N,Res[0:l,0],Res[0:l,1]*nmax)
#modelRes=N(Res[0:l,0],*poRes)

#plt.title(r'J=3 $k_BT$, $\mu$=-3 $k_BT$',fontsize=17)
plt.xlabel('Simulation time',fontsize=14)
plt.ylabel('Stator number',fontsize=14)
plt.rcParams["figure.figsize"] = [8.0,6.0]

#Representation of T trajectories
for i in range(0,T):
    plt.plot(All_Res[(S*i):(S*(i+1)-m),1],All_Res[(S*i):(S*(i+1)-m),2]*nmax, alpha=0.2)

#tau_Res=poRes[0]

plt.plot(Res[0:l,0],Res[0:l,1]*nmax)    
#plt.plot(Res[0:l,0],modelRes, label=r'$\tau$={:1.2f}'.format(tau_Res))

#plt.legend(loc='right', prop={'size': 15})


#tau_seq_Res[3] = poRes[0]
#ns_seq_Res[3] = poRes[1]

#tau_seq_Res_er[3]=np.sqrt(pcRes[0,0])

#%% Stall representation. Fit

n0=8

poStall, pcStall = curve_fit(N,Stall[0:l,0],Stall[0:l,1]*nmax)
modelStall=N(Stall[0:l,0],*poStall)

len(modelStall)

plt.title(r'J=0 $k_BT$, $\mu$=-0.69 $k_BT$',fontsize=17)
plt.xlabel('Simulation time',fontsize=14)
plt.ylabel('Stator number',fontsize=14)
plt.rcParams["figure.figsize"] = [8.0,6.0]

#Representation of T trajectories
for i in range(0,T):
    plt.plot(All_Stall[(S*i):(S*(i+1)-m),1],All_Stall[(S*i):(S*(i+1)-m),2]*nmax, alpha=0.2)

plt.plot(Stall[0:l,0],Stall[0:l,1]*nmax)    
plt.plot(Stall[0:l,0],modelStall, label=r'$\tau$={:1.2f}'.format(poStall[0]))
plt.legend(loc='lower right', prop={'size': 15})

#print(poStall[0])

#tau_seq_Stall[3] = poStall[0]
#ns_seq_Stall[3] = poStall[1]

#tau_seq_Stall_er[3]=np.sqrt(pcStall[0,0])
"""
#%%
plt.plot(Res[0:l,0],Res[0:l,1]*nmax)    
plt.plot(Stall[0:l,0],Stall[0:l,1]*nmax)

#%%

def Nl(x,kon,koff,n0):
    nss = nmax/(1+koff/kon)
    return nss + (n0-nss)*np.exp(-x*(kon+koff)) 

#%%
    
kon = 2.47e-3
koff= 1.57e-3

t = np.linspace(0,4000,10000)

plt.grid()

plt.plot(t,Nl(t,kon,koff,0))
plt.plot(t,Nl(t,kon,koff,9.25))

nmax/(1+koff/kon)

#%%

n0=0


poRes, pcRes = curve_fit(N,Res[0:l,0],Res[0:l,1]*nmax)
modelRes=N(Res[0:l,0],*poRes)
  

#%%

n0=9.25

poStall, pcStall = curve_fit(N,Stall[0:l,0],Stall[0:l,1]*nmax)
modelStall=N(Stall[0:l,0],*poStall)

#%%

plt.title('500 nm bead',fontsize=14)
plt.xlabel('Time',fontsize=16)
plt.ylabel('Stators',fontsize=16)
plt.rcParams["figure.figsize"] = [8.0,6.0]
plt.grid()


plt.plot(Res[0:l,0],Res[0:l,1]*nmax,label='Glauber',linestyle='--',color='black')    
plt.plot(Stall[0:l,0],Stall[0:l,1]*nmax,linestyle='--',color='black')
plt.plot(t,Nl(t,kon,koff,0),color='gray',label='Langmuir')
plt.plot(t,Nl(t,kon,koff,9.25),color='gray')
#plt.plot(Stall[0:l], modelStall)
#plt.plot(Res[0:l], modelRes)


plt.legend()

#%%
tau_res = np.zeros(6)
tau_stall = np.zeros(6)
#%%

tau_res[4] = poRes[0]

tau_stall[4]= poStall[0]

#%%

J=np.arange(0,6,1)

plt.title('$<\phi>=1/2$, $\lambda = 0.046$, $\Gamma = 300$')
plt.xlabel('J ($k_B T$)')
plt.ylabel('$t_c$')
plt.plot(J,tau_res,linestyle='',marker='o',label='Resurrection')
plt.plot(J,tau_stall,linestyle='',marker='o',label='Release')

plt.legend()

#%%

files9=['/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/Initial_conditions_J=4.00  .dat'] #Average of stall simulations
data9=[]

for data_file in files9:
    data9.append(np.loadtxt(data_file))
    
IC=data9[0]*nmax

files10=['/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/Statistics_Initial_contidions_J=0.00  .dat'] #Average of stall simulations
data10=[]

files101=['/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/Statistics_Initial_contidions_J=0.25  .dat'] #Average of stall simulations
data101=[]
files102=['/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/Statistics_Initial_contidions_J=0.50  .dat'] #Average of stall simulations
data102=[]
files103=['/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/Statistics_Initial_contidions_J=0.75  .dat'] #Average of stall simulations
data103=[]

files11=['/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/Statistics_Initial_contidions_J=1.00  .dat'] #Average of stall simulations
data11=[]

files111=['/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/Statistics_Initial_contidions_J=1.25  .dat'] #Average of stall simulations
data111=[]
files112=['/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/Statistics_Initial_contidions_J=1.50  .dat'] #Average of stall simulations
data112=[]
files113=['/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/Statistics_Initial_contidions_J=1.75  .dat'] #Average of stall simulations
data113=[]

files12=['/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/Statistics_Initial_contidions_J=2.00  .dat'] #Average of stall simulations
data12=[]

files121=['/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/Statistics_Initial_contidions_J=2.25  .dat'] #Average of stall simulations
data121=[]
files122=['/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/Statistics_Initial_contidions_J=2.50  .dat'] #Average of stall simulations
data122=[]
files123=['/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/Statistics_Initial_contidions_J=2.75  .dat'] #Average of stall simulations
data123=[]

files13=['/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/Statistics_Initial_contidions_J=3.00  .dat'] #Average of stall simulations
data13=[]

files131=['/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/Statistics_Initial_contidions_J=3.25  .dat'] #Average of stall simulations
data131=[]
files132=['/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/Statistics_Initial_contidions_J=3.50  .dat'] #Average of stall simulations
data132=[]
files133=['/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/Statistics_Initial_contidions_J=3.75  .dat'] #Average of stall simulations
data133=[]

files14=['/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/Statistics_Initial_contidions_J=4.00  .dat'] #Average of stall simulations
data14=[]
files141=['/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/Statistics_Initial_contidions_J=4.25  .dat'] #Average of stall simulations
data141=[]
files142=['/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/Statistics_Initial_contidions_J=4.50  .dat'] #Average of stall simulations
data142=[]
files143=['/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/Statistics_Initial_contidions_J=4.75  .dat'] #Average of stall simulations
data143=[]

files15=['/home/mariajose/Escritorio/Simulations/Without depletion/Glauber/Statistics_Initial_contidions_J=5.00  .dat'] #Average of stall simulations
data15=[]

for data_file in files10:
    data10.append(np.loadtxt(data_file))
for data_file in files101:
    data101.append(np.loadtxt(data_file))    
for data_file in files102:
    data102.append(np.loadtxt(data_file))
for data_file in files103:
    data103.append(np.loadtxt(data_file))

for data_file in files11:
    data11.append(np.loadtxt(data_file))
for data_file in files111:
    data111.append(np.loadtxt(data_file))
for data_file in files112:
    data112.append(np.loadtxt(data_file))
for data_file in files113:
    data113.append(np.loadtxt(data_file))

for data_file in files12:
    data12.append(np.loadtxt(data_file))
for data_file in files121:
    data121.append(np.loadtxt(data_file))
for data_file in files122:
    data122.append(np.loadtxt(data_file))
for data_file in files123:
    data123.append(np.loadtxt(data_file))
    
for data_file in files13:
    data13.append(np.loadtxt(data_file))
for data_file in files131:
    data131.append(np.loadtxt(data_file))
for data_file in files132:
    data132.append(np.loadtxt(data_file))
for data_file in files133:
    data133.append(np.loadtxt(data_file))
    
for data_file in files14:
    data14.append(np.loadtxt(data_file))
for data_file in files141:
    data141.append(np.loadtxt(data_file))
for data_file in files142:
    data142.append(np.loadtxt(data_file))
for data_file in files143:
    data143.append(np.loadtxt(data_file))

for data_file in files15:
    data15.append(np.loadtxt(data_file))



stats_0 = data10[0]
stats_025 = data101[0]
stats_050 = data102[0]
stats_075 = data103[0]

stats_1 = data11[0]
stats_125 = data111[0]
stats_150 = data112[0]
stats_175 = data113[0]

stats_2 = data12[0]
stats_225 = data121[0]
stats_250 = data122[0]
stats_275 = data123[0]

stats_3 = data13[0]
stats_325 = data131[0]
stats_350 = data132[0]
stats_375 = data133[0]

stats_4 = data14[0]
stats_425 = data141[0]
stats_450 = data142[0]
stats_475 = data143[0]

stats_5 = data15[0]

#%%

u_IC, co_IC = np.unique(IC,return_counts=True)

plt.title('$<N_0> = $, $J = 5 \ k_B T$',size=17)
plt.rcParams["figure.figsize"] = [8.0,6.0]
#plt.xlabel()
plt.bar(u_IC,co_IC/len(IC),width=1)
plt.vlines(9, 0, 0.26, linestyles='dashed')


#%%

J=np.arange(0,5.25,0.25)

skew=np.array([stats_0[2],stats_025[2],stats_050[2],stats_075[2],stats_1[2],stats_125[2],stats_150[2],stats_175[2],
              stats_2[2],stats_225[2],stats_250[2],stats_275[2],stats_3[2],stats_325[2],stats_350[2],stats_375[2],
              stats_4[2],stats_425[2],stats_450[2],stats_475[2],stats_5[2]])

kurt=np.array([stats_0[3],stats_025[3],stats_050[3],stats_075[3],stats_1[3],stats_125[3],stats_150[3],stats_175[3],
              stats_2[3],stats_225[3],stats_250[3],stats_275[3],stats_3[3],stats_325[3],stats_350[3],stats_375[3],
              stats_4[3],stats_425[3],stats_450[3],stats_475[3],stats_5[3]])

plt.title('Initial conditions $<n_0>=4$',size=17)
plt.xlabel('J / $k_B T$',size=14)
plt.plot(J,skew,linestyle='',marker='.',label='Skewness')
plt.plot(J,kurt,linestyle='',marker='.',label='Kurtosis')

plt.legend(prop={'size': 15})

#%%

J=np.linspace(0, 4,5)

plt.title('$\mu$= 0 $k_BT$',fontsize=17)
plt.xlabel('J / kT',fontsize=14)
plt.ylabel('Resurrection time',fontsize=14)
plt.errorbar(J,tau_seq_Res,tau_seq_Res_er,marker='.',label='Resurrection',capsize=5)
plt.errorbar(J,tau_seq_Stall,tau_seq_Stall_er,marker='.',label='Stall',capsize=5)

plt.legend(prop={'size': 15})


np.savetxt('/home/mariajose/Escritorio/Simulations/Without depletion/Relaxation times/ResTime_J=-mu.dat',(np.c_[J,tau_seq_Res,tau_seq_Res_er]))
np.savetxt('/home/mariajose/Escritorio/Simulations/Without depletion/Relaxation times/StallTime_J=-mu.dat',(np.c_[J,tau_seq_Stall,tau_seq_Stall_er]))

#%% Representation of stall and resurrection togheter

plt.title(r'J=3 $k_BT$, $\mu$=-3 $k_BT$',fontsize=17)
plt.xlabel('Simulation time',fontsize=14)
plt.ylabel('Stator number',fontsize=14)
plt.rcParams["figure.figsize"] = [8.0,6.0]

plt.plot(Res[0:l,0],Res[0:l,1]*nmax,label=r'$\tau=137.64$')
plt.plot(Stall[0:l,0],Stall[0:l,1]*nmax,label=r'$\tau=138.78$')

#plt.hlines(8,0,l,alpha=0.3,linestyles='dashed', label=r'$N_{ss}$=8')

plt.legend(loc='lower right', prop={'size': 15})

plt.show()

#%% Representation of standard deviation

plt.title('J=0 kT, $\mu=1$ kT')
plt.xlabel('Simulation time')
plt.ylabel('$(<\phi^2>-<\phi>^2)^{1/2}$')
plt.plot(Res[0:l,0],Res[0:l,2],label='Resurrection')
plt.plot(Stall[0:l,0],Stall[0:l,2],label='Stall')
plt.hlines(0.13,0,200,label='Relaxation value=0.13',linestyles='dashed')

plt.legend(loc='upper right')