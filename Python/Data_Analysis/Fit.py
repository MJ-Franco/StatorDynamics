#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 15:45:24 2020

@author: mariajose

Study of the experimental fluctuations
"""
import numpy as np
import matplotlib.pyplot as plt
FPS=1000.
n_steps=2000


#%% Ashley program

def LoadExpData():
    '''Loads the Experimental Data.  Change the pathnames below to match your comp'''
    
    D300 = np.load('/home/mariajose/Escritorio/Master PIV/Stage/Data/NilsD300_V2.p',allow_pickle=True)
    D500 = np.load('/home/mariajose/Escritorio/Master PIV/Stage/Data/NilsD500.p',allow_pickle=True)
    D1300 = np.load('/home/mariajose/Escritorio/Master PIV/Stage/Data/NilsD1300.p',allow_pickle=True)
    return D300, D500, D1300



def CalculateAvgs(D300, D500, D1300):
    '''
    Calculates the averages and standard deviations of the data for each load size. 'Before' 
    refers to before stall, i.e. steady state.  'After' refers to immediately after the motor
    is released from stall.  'Resurrection' refers to immediately after fresh motility buffer 
    has been added, allowing for the stators to return.
    Returns a dictionary.
    '''
    
    _,_,_,_,_,_,N300_before,N300_after,N300_resurrection = Dict2DataMatrix(D300)
    N300_before[N300_before==0],N300_after[N300_after==0] = np.nan, np.nan
    N300_before_avg, N300_after_avg, N300_resurrection_avg = np.array(np.nanmean(N300_before,0)),np.array(np.nanmean(N300_after,0)),np.array(np.nanmean(N300_resurrection,0))
    N300_before_std, N300_after_std, N300_resurrection_std = np.array(np.nanstd(N300_before,0)),np.array(np.nanstd(N300_after,0)),np.array(np.nanstd(N300_resurrection,0))
    _,_,_,_,_,_,N500_before,N500_after,N500_resurrection = Dict2DataMatrix(D500)
    N500_before[N500_before==0],N500_after[N500_after==0] = np.nan, np.nan
    N500_before_avg, N500_after_avg, N500_resurrection_avg = np.array(np.nanmean(N500_before,0)),np.array(np.nanmean(N500_after,0)),np.array(np.nanmean(N500_resurrection,0))
    N500_before_std, N500_after_std, N500_resurrection_std = np.array(np.nanstd(N500_before,0)),np.array(np.nanstd(N500_after,0)),np.array(np.nanstd(N500_resurrection,0))
    _,_,_,_,_,_,N1300_before,N1300_after,N1300_resurrection = Dict2DataMatrix(D1300)
    N1300_before[N1300_before==0],N1300_after[N1300_after==0] = np.nan, np.nan
    N1300_before_avg, N1300_after_avg, N1300_resurrection_avg = np.array(np.nanmean(N1300_before,0)),np.array(np.nanmean(N1300_after,0)),np.array(np.nanmean(N1300_resurrection,0))
    N1300_before_std, N1300_after_std, N1300_resurrection_std = np.array(np.nanstd(N1300_before,0)),np.array(np.nanstd(N1300_after,0)),np.array(np.nanstd(N1300_resurrection,0))
    
    AvgDict = {'N300_before_avg':N300_before_avg, 'N300_after_avg':N300_after_avg, \
               'N300_resurrection_avg':N300_resurrection_avg, 'N300_before_std':N300_before_std, \
               'N300_after_std':N300_after_std, 'N300_resurrection_std':N300_resurrection_std, \
               'N500_before_avg':N500_before_avg, 'N500_after_avg':N500_after_avg, \
               'N500_resurrection_avg':N500_resurrection_avg, 'N500_before_std':N500_before_std, \
               'N500_after_std':N500_after_std, 'N500_resurrection_std':N500_resurrection_std, \
               'N1300_before_avg':N1300_before_avg, 'N1300_after_avg':N1300_after_avg, \
               'N1300_resurrection_avg':N1300_resurrection_avg, 'N1300_before_std':N1300_before_std, \
               'N1300_after_std':N1300_after_std, 'N1300_resurrection_std':N1300_resurrection_std}
    
    return AvgDict



def PlotExpData(AvgDict):
    #cut resurrections and after stalls at 800 s
    time_lim = int(800 *FPS)
    
    t = np.arange(len(AvgDict['N300_before_avg']))/FPS
    plt.plot(t,AvgDict['N300_before_avg'],'k',label='300nm bead steady-state')
    t = np.arange(len(AvgDict['N300_resurrection_avg']))/FPS
    plt.plot(t[:time_lim],AvgDict['N300_resurrection_avg'][:time_lim],'r',label='300nm bead resurrecion')
    t = np.arange(len(AvgDict['N300_after_avg']))/FPS
    plt.plot(t[:time_lim],AvgDict['N300_after_avg'][:time_lim],'r',label='300nm bead after stall')
    
    t = np.arange(len(AvgDict['N500_before_avg']))/FPS
    plt.plot(t,AvgDict['N500_before_avg'],'k',label='500nm bead steady-state')
    t = np.arange(len(AvgDict['N500_resurrection_avg']))/FPS
    plt.plot(t[:time_lim],AvgDict['N500_resurrection_avg'][:time_lim],'g',label='500nm bead resurrecion')
    t = np.arange(len(AvgDict['N500_after_avg']))/FPS
    plt.plot(t[:time_lim],AvgDict['N500_after_avg'][:time_lim],'g',label='500nm bead after stall')
    
    t = np.arange(len(AvgDict['N1300_before_avg']))/FPS
    plt.plot(t,AvgDict['N1300_before_avg'],'k',label='1300nm bead steady-state')
    t = np.arange(len(AvgDict['N1300_resurrection_avg']))/FPS
    plt.plot(t[:time_lim],AvgDict['N1300_resurrection_avg'][:time_lim],'b',label='1300nm bead resurrecion')
    t = np.arange(len(AvgDict['N1300_after_avg']))/FPS
    plt.plot(t[:time_lim],AvgDict['N1300_after_avg'][:time_lim],'b',label='1300nm bead after stall')
    
    plt.autoscale(enable=True, axis='both', tight=True)
    plt.xlabel('time (s)')
    plt.ylabel('# stators')
    plt.legend()
    
    return

  
def Dict2DataMatrix(D):
    
    max_before = np.max([len(D[k]['torque_before_stall_pNnm']) for k in D.keys()])
    max_after = np.max([len(D[k]['torque_after_release_pNnm']) for k in D.keys()])
    max_resurrection = np.max([len(D[k]['torque_resurrection_pNnm']) for k in D.keys() if 'torque_resurrection_pNnm' in D[k].keys()])
        
    T_before = np.ones((np.max([*D.keys()])+1,max_before)) * np.nan
    T_after = np.ones((np.max([*D.keys()])+1,max_after)) * np.nan
    T_resurrection = np.ones((np.max([*D.keys()])+1,max_resurrection)) * np.nan
    S_before = np.ones((np.max([*D.keys()])+1,max_before)) * np.nan
    S_after = np.ones((np.max([*D.keys()])+1,max_after)) * np.nan
    S_resurrection = np.ones((np.max([*D.keys()])+1,max_resurrection)) * np.nan
    N_before = np.ones((np.max([*D.keys()])+1,max_before)) * np.nan
    N_after = np.ones((np.max([*D.keys()])+1,max_after)) * np.nan
    N_resurrection = np.ones((np.max([*D.keys()])+1,max_resurrection)) * np.nan
    
    for k in D.keys():
        drag = D[k]['drag_Nms']
        sign_before = np.sign(np.mean(D[k]['torque_before_stall_pNnm']))
        sign_after = np.sign(np.mean(D[k]['torque_after_release_pNnm']))
        T_before[k,max_before-len(D[k]['torque_before_stall_pNnm']):] = D[k]['torque_before_stall_pNnm'] * sign_before
        T_after[k,:len(D[k]['torque_after_release_pNnm'])] = D[k]['torque_after_release_pNnm'] * sign_after
        S_before[k,max_before-len(D[k]['torque_before_stall_pNnm']):] = D[k]['torque_before_stall_pNnm'] / (2*np.pi*drag*1e21) * sign_before
        S_after[k,:len(D[k]['torque_after_release_pNnm'])] = D[k]['torque_after_release_pNnm']/ (2*np.pi*drag*1e21) * sign_after
        N_before[k,max_before-len(D[k]['torque_before_stall_pNnm']):] = D[k]['statnum_before_stall']
        N_after[k,:len(D[k]['torque_after_release_pNnm'])] = D[k]['statnum_after_release']
        if 'fit_speed_resurrection_Hz' in D[k].keys():
            sign_resurrection = np.sign(np.mean(D[k]['torque_resurrection_pNnm']))
            T_resurrection[k,:len(D[k]['torque_resurrection_pNnm'])] = D[k]['torque_resurrection_pNnm'] * sign_resurrection
            S_resurrection[k,:len(D[k]['torque_resurrection_pNnm'])] = D[k]['torque_resurrection_pNnm']/ (2*np.pi*drag*1e21) * sign_resurrection
            N_resurrection[k,:len(D[k]['torque_resurrection_pNnm'])] = D[k]['statnum_resurrection'] 
        
    return T_before,T_after,T_resurrection,S_before,S_after,S_resurrection,N_before,N_after,N_resurrection

#Load the experimental data:
D300, D500, D1300 = LoadExpData()
#   Calculate the averages, return a dictionary:
AvgDict = CalculateAvgs(D300, D500, D1300)

#%% Lecture of the simulation data for resurrection

files1=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Res.dat'] #All resurrection simulations
data1=[]
files2=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Res_Avg.dat'] #Average of resurrection simulations
data2=[]
files3=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Stall.dat'] #All resurrection simulations
data3=[]
files4=['/home/mariajose/Escritorio/Simulations/Without depletion/Traces/Traces_Stall_Avg.dat'] #Average of resurrection simulations
data4=[]


for data_file in files1:
    data1.append(np.loadtxt(data_file))

for data_file in files2:
    data2.append(np.loadtxt(data_file))

for data_file in files3:
    data3.append(np.loadtxt(data_file))

for data_file in files4:
    data4.append(np.loadtxt(data_file))
    
All_Res=data1[0]
Res=data2[0]

All_Stall=data3[0]
Stall=data4[0]

 #simulation step number
#%% Lecture and representation of the data from experiments of resurrection 300 nm

Nmax=12

Res300=AvgDict['N300_resurrection_avg']/Nmax

Res300_std = AvgDict['N300_resurrection_std']/Nmax

data_cut = int(1000*FPS)

Res300_cut=np.zeros(data_cut)
Res300_std_cut=np.zeros(data_cut)

for i in range(data_cut):
    Res300_cut[i] = Res300[i]
    Res300_std_cut[i] = Res300_std[i]

t=np.arange(len((Res300_cut)))/FPS

time_300 = t[len(Res300_cut)-1]

Res300_avg=sum(Res300_cut)/len(Res300_cut)

plt.title('Resurrection 300 nm bead ($N_{max}=12$)')
plt.plot(t,Res300_cut)
plt.hlines(Res300_avg,0,time_300,label='$<\phi>$={:1.2f}'.format(Res300_avg),alpha=0.3,linestyles='dashed')

plt.legend()

#%% Constantes y vectores para separar en ventanas los datos 300nm

#t_windows=time_300/n_steps #Original time window

t_windows=10

data_windows= int(t_windows*FPS) #number of data in each window

n_windows=int(time_300/t_windows) #number of windows

t = np.linspace(0,time_300,n_windows)  #time

#%% Calculation of matrices for separation in windows of the data


Res300_windows = np.zeros((data_windows,n_windows))
Res300_std_windows = np.zeros((data_windows,n_windows))

for i in range(n_windows):
    for j in range(data_windows):

        Res300_windows[j,i] = Res300_cut[j+data_windows*i]
        Res300_std_windows[j,i] = Res300_std_cut[j+data_windows*i]

    
Res300_std_sum = 0

AVbeg=int(100*FPS)

for i in range(AVbeg,len(Res300_std_cut)):
    Res300_std_sum = Res300_std_sum + Res300_std_cut[i]

Res_std_avg = Res300_std_sum/len(Res300_std_cut)

Res_std_avg


#%% Representation of the evolution of the standard deviation

plt.title('Resurrection, $<\phi>=0.30$, $N_{max}=12$')
plt.xlabel('t / s')
plt.ylabel('$(<\phi^2> - <\phi>^2)^{1/2}$')

for i in range(n_windows):   
    plt.plot(t[i*data_windows:(i+1)*data_windows], Res300_std_windows[i,i*data_windows:(i+1)*data_windows],color='blue')
    
    if(i==1):
        plt.hlines(Res_std_avg, 0, time_300,alpha=0.3,linestyles='dashed',label='Average={:1.2f}'.format(Res_std_avg))
        plt.legend()
        

#%% Lecture and representation of the data of resurrection 500 nm bead

Nmax=12

Res500=AvgDict['N500_resurrection_avg']/Nmax
Res500_std = AvgDict['N500_resurrection_std']/Nmax

data_cut = int(1000*FPS)

Res500_cut=np.zeros(data_cut)
Res500_std_cut=np.zeros(data_cut)

for i in range(data_cut):
    Res500_cut[i] = Res500[i]
    Res500_std_cut[i] = Res500_std[i]


t=np.arange(len((Res500_cut)))/FPS

time_500 = t[len(Res500_cut)-1]

Res500_avg=sum(Res500_cut)/len(Res500_cut)

plt.title('Resurrection 500 nm bead ($N_{max}=12$)')
plt.xlabel('t / s')
plt.ylabel('$<\phi>$')
plt.plot(t,Res500_cut)
plt.hlines(Res500_avg,0,len(t)/FPS,alpha=0.3,linestyles='dashed',label='$<\phi>$={:1.2f}'.format(Res500_avg))

plt.legend()

#%%500 nm bead

#t_windows= time_500/n_steps

t_windows=10

data_windows= int(t_windows*FPS) #number of data in each window

n_windows=int(time_500/t_windows) #number of windows

t = np.linspace(0,time_500,n_windows)  #time
#%% Calculation of matrices for separation in windows of the data 500 nm


Res500_windows = np.zeros((data_windows,n_windows))
Res500_std_windows = np.zeros((data_windows,n_windows))

for i in range(n_windows):
    for j in range(data_windows):

        Res500_windows[j,i] = Res500_cut[j+data_windows*i]
        Res500_std_windows[j,i] = Res500_std_cut[j+data_windows*i]

AVbeg=int(100*FPS)

Res500_std_sum=0
for i in range(AVbeg,len(Res500_std_cut)):
    Res500_std_sum = Res500_std_sum + Res500_std_cut[i]

Res500_std_avg = Res500_std_sum/len(Res500_std_cut)


#%% Representation of the evolution of the standard deviation


plt.title('Resurrection, $<\phi>=0.52$, $N_{max}=12$')
plt.xlabel('t / s')
plt.ylabel('$(<\phi^2> - <\phi>^2)^{1/2}$')

for i in range(n_windows):
    
    plt.plot(t[i*data_windows:(i+1)*data_windows], Res500_std_windows[i,i*data_windows:(i+1)*data_windows],color='blue')
    
    if(i==1):
        plt.hlines(Res500_std_avg, 0, time_500,alpha=0.3,linestyles='dashed',label='Average={:1.2f}'.format(Res500_std_avg))
        plt.legend()
  

  
#%% 1300 nm  resurrection

Nmax=12

Res1300=AvgDict['N1300_resurrection_avg']/Nmax
Res1300_std = AvgDict['N1300_resurrection_std']/Nmax

data_cut = int(1000*FPS)

Res1300_cut=np.zeros(data_cut)
Res1300_std_cut=np.zeros(data_cut)

for i in range(data_cut):
    Res1300_cut[i] = Res1300[i]
    Res1300_std_cut[i] = Res1300_std[i]

t=np.arange(len((Res1300_cut)))/FPS

time_1300 = t[len(Res1300_cut)-1]

Res1300_avg=sum(Res1300_cut)/len(Res1300_cut)


plt.title('Resurrection 1300 nm bead ($N_{max}=12$)')
plt.xlabel('t / s')
plt.ylabel('$<\phi>$')
plt.plot(t,Res1300_cut)
plt.hlines(Res1300_avg,0,len(t)/FPS,alpha=0.3,linestyles='dashed',label='$<\phi>$={:1.2f}'.format(Res1300_avg))

plt.legend()

#%%1300 nm bead resurrection

#t_windows= time_1300/n_steps

t_windows=10

data_windows= int(t_windows*FPS) #number of data in each window

n_windows=int(time_1300/t_windows) #number of windows

t = np.linspace(0,time_1300,n_windows)  #time

#%% Calculation of matrices for separation in windows of the data 500 nm

Res1300_windows = np.zeros((data_windows,n_windows))
Res1300_std_windows = np.zeros((data_windows,n_windows))

for i in range(n_windows):
    for j in range(data_windows):

        Res1300_windows[j,i] = Res1300_cut[j+data_windows*i]
        Res1300_std_windows[j,i] = Res1300_std_cut[j+data_windows*i]


AVbeg=int(100*FPS)

Res1300_std_sum=0
for i in range(AVbeg,len(Res500_std_cut)):
    Res1300_std_sum = Res1300_std_sum + Res1300_std_cut[i]

Res1300_std_avg = Res1300_std_sum/len(Res1300_std_cut)


#%% Representation of the evolution of the standard deviation

plt.title('Resurrection, $<\phi>=0.74$, $N_{max}=12$')
plt.xlabel('t / s')
plt.ylabel('$(<\phi^2> - <\phi>^2)^{1/2}$')

for i in range(n_windows):
    
    plt.plot(t[i*data_windows:(i+1)*data_windows], Res1300_std_windows[i,i*data_windows:(i+1)*data_windows],color='blue')
    
    if(i==1):
        plt.hlines(Res1300_std_avg, 0, time_1300,alpha=0.3,linestyles='dashed',label='Average={:1.2f}'.format(Res1300_std_avg))
        plt.legend()
  
 
plt.legend()
plt.show()

#%% 300nm stall

Nmax=12

Stall300=AvgDict['N300_after_avg']/Nmax
Stall300_std = AvgDict['N300_after_std']/Nmax

t=np.arange(len((Stall300)))/FPS

time_300 = t[len(Stall300)-1]

Stall300_avg=sum(Stall300)/len(Stall300)

plt.plot(t,Stall300)
plt.hlines(Stall300_avg,0,time_300,label=Stall300_avg,alpha=0.3,linestyles='dashed')

plt.legend()


#%% Constantes y vectores para separar en ventanas los datos 300nm

#t_windows=time_300/n_steps #time window

t_windows=10

data_windows= int(t_windows*FPS) #number of data in each window

n_windows=int(time_300/t_windows) #number of windows

t = np.linspace(0,time_300,n_windows)  #time
#%% Calculation of matrices for separation in windows of the data


Stall300_windows = np.zeros((data_windows,n_windows))
Stall300_std_windows = np.zeros((data_windows,n_windows))

for i in range(n_windows):
    for j in range(data_windows):

        Stall300_windows[j,i] = Stall300[j+data_windows*i]
        Stall300_std_windows[j,i] = Stall300_std[j+data_windows*i]


Stall300_std_avg = sum(Stall300_std)/len(Stall300_std)

Stall300_std_avg


#%% Representation of the evolution of the standard deviation

plt.title('Stall $<\phi>=0.3$, $N_{max}=12$')
plt.xlabel('t / s')
plt.ylabel('$(<\phi^2> - <\phi>^2)^{1/2}$')

for i in range(n_windows):
    
    plt.plot(t[i*data_windows:(i+1)*data_windows], Stall300_std_windows[i,i*data_windows:(i+1)*data_windows],color='blue')
    
    if(i==1):
        plt.hlines(Stall300_std_avg, 0, time_300,alpha=0.3,linestyles='dashed',label='Average={:1.2f}'.format(Stall300_std_avg))
        plt.legend()

#%% 500 nm stall

Nmax=12

Stall500=AvgDict['N500_after_avg']
Stall500_std = AvgDict['N500_after_std']/Nmax

t=np.arange(len((Stall500)))/FPS

time_500 = t[len(Stall500)-1]

Stall500_avg=sum(Stall500)/len(Stall500)

plt.plot(t,Stall500)
plt.hlines(Stall500_avg,0,time_500,label=Stall500_avg,alpha=0.3,linestyles='dashed')

plt.legend()


#%% Constantes y vectores para separar en ventanas los datos 300nm

#t_windows=time_500/n_steps #time window

t_windows=10

data_windows= int(t_windows*FPS) #number of data in each window

n_windows=int(time_500/t_windows) #number of windows

t = np.linspace(0,time_500,n_windows)  #time

#%% Calculation of matrices for separation in windows of the data


Stall500_windows = np.zeros((data_windows,n_windows))
Stall500_std_windows = np.zeros((data_windows,n_windows))

for i in range(n_windows):
    for j in range(data_windows):

        Stall500_windows[j,i] = Stall500[j+data_windows*i]
        Stall500_std_windows[j,i] = Stall500_std[j+data_windows*i]


Stall500_std_avg = sum(Stall500_std)/len(Stall500_std)

Stall500_std_avg


#%% Representation of the evolution of the standard deviation

plt.title('Stall, $<\phi>=0.52$, $N_{max}=12$')
plt.xlabel('t / s')
plt.ylabel('$(<\phi^2> - <\phi>^2)^{1/2}$')

for i in range(n_windows):
    
    plt.plot(t[i*data_windows:(i+1)*data_windows], Stall500_std_windows[i,i*data_windows:(i+1)*data_windows],color='blue')
    
    if(i==1):
        plt.hlines(Stall500_std_avg, 0, time_500,alpha=0.3,linestyles='dashed',label='Average={:1.2f}'.format(Stall500_std_avg))
        plt.legend()

#%% 1300 nm stall

Nmax=12

Stall1300=AvgDict['N1300_after_avg']

Stall1300_std = AvgDict['N1300_after_std']/Nmax

t=np.arange(len((Stall1300)))/FPS

time_1300 = t[len(Stall1300)-1]

Stall1300_avg=sum(Stall1300)/len(Stall1300)

plt.plot(t,Stall1300)
plt.hlines(Stall1300_avg,0,time_1300,label=Stall1300_avg,alpha=0.3,linestyles='dashed')

plt.legend()


#%% Constantes y vectores para separar en ventanas los datos 300nm

#t_windows=time_1300/n_steps #time window

t_windows=10

data_windows= int(t_windows*FPS) #number of data in each window

n_windows=int(time_1300/t_windows) #number of windows

t = np.linspace(0,time_1300,n_windows)  #time

#%% Calculation of matrices for separation in windows of the data


Stall1300_windows = np.zeros((data_windows,n_windows))
Stall1300_std_windows = np.zeros((data_windows,n_windows))

for i in range(n_windows):
    for j in range(data_windows):

        Stall1300_windows[j,i] = Stall1300[j+data_windows*i]
        Stall1300_std_windows[j,i] = Stall1300_std[j+data_windows*i]


Stall1300_std_avg = sum(Stall1300_std)/len(Stall1300_std)

Stall1300_std_avg


#%% Representation of the evolution of the standard deviation

plt.title('Stall, $<\phi>=0.74$, $N_{max}=12$')
plt.xlabel('t / s')
plt.ylabel('$(<\phi^2> - <\phi>^2)^{1/2}$')

for i in range(n_windows):
    
    plt.plot(t[i*data_windows:(i+1)*data_windows], Stall1300_std_windows[i,i*data_windows:(i+1)*data_windows],color='blue')
    
    if(i==1):
        plt.hlines(Stall1300_std_avg, 0, time_1300,alpha=0.3,linestyles='dashed',label='Average={:1.2f}'.format(Stall1300_std_avg))
        plt.legend()



