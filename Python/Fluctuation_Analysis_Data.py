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
#   Plot the averages:
#PlotExpData(AvgDict)



#%% Data resurrection 300 nm bead

Res300=AvgDict['N300_resurrection_avg']

Res300_std=AvgDict['N300_resurrection_std']

AVbeg = int(200*FPS)
AVend = int(1000*FPS)

SS_data = AVend - AVbeg

Nsum=0
for i in range(AVbeg,AVend):
    Nsum = Nsum + Res300[i]

Res_300_avg=Nsum/SS_data

stdsum=0
for i in range(AVbeg,AVend):
    stdsum= stdsum + Res300_std[i]

Res300_std_avg=stdsum/SS_data

Res300_var_avg = Res300_std_avg**2

errorsum=0
avgsum=0
for i in range(AVbeg,AVend):
    errorsum= errorsum + (Res300_std[i]-Res300_std_avg)**2
    avgsum = avgsum + (Res300[i]-Res_300_avg)**2
    
Res_300_avg_error = np.sqrt(avgsum/SS_data)   
Res_300_std_error = np.sqrt(errorsum/SS_data)

#%% Data stall 300 nm bead

Stall300=AvgDict['N300_after_avg']

Stall300_std=AvgDict['N300_after_std']

Nsum=0
for i in range(len(Stall300)):
    Nsum = Nsum + Stall300[i]
    
Stall_300_avg = Nsum/len(Stall300)

stdsum=0
for i in range(len(Stall300)):
    stdsum = stdsum + Stall300_std[i]
    
Stall_300_std_avg = stdsum/len(Stall300)

Stall_300_var_std = Stall_300_std_avg**2

errorsum=0
avgsum=0
for i in range(len(Stall300)):
    errorsum= errorsum + (Stall300_std[i]-Stall_300_std_avg)**2
    avgsum = avgsum + (Stall300[i]-Stall_300_avg)**2
    
Stall_300_avg_error = np.sqrt(avgsum/len(Stall300))   
Stall_300_std_error = np.sqrt(errorsum/len(Stall300))

#%% Data resurrection 500 nm bead

Res500=AvgDict['N500_resurrection_avg']

Res500_std=AvgDict['N500_resurrection_std']

AVbeg = int(200*FPS)
AVend = int(1000*FPS)

SS_data=AVend - AVbeg

Nsum=0
for i in range(AVbeg,AVend):
    Nsum = Nsum + Res500[i]

Res_500_avg=Nsum/SS_data

stdsum=0
for i in range(AVbeg,AVend):
    stdsum= stdsum + Res500_std[i]

Res500_std_avg=stdsum/SS_data

Res500_var_avg = Res500_std_avg**2

errorsum=0
avgsum=0
for i in range(AVbeg,AVend):
    errorsum= errorsum + (Res500_std[i]-Res500_std_avg)**2
    avgsum = avgsum + (Res500[i]-Res_500_avg)**2
    
Res_500_avg_error = np.sqrt(avgsum/SS_data)   
Res_500_std_error = np.sqrt(errorsum/SS_data)


#%% Data stall 500 nm bead

Stall500=AvgDict['N500_after_avg']

Stall500_std=AvgDict['N500_after_std']

Nsum=0
for i in range(len(Stall500)):
    Nsum = Nsum + Stall500[i]
    
Stall_500_avg = Nsum/len(Stall500)

stdsum=0
for i in range(len(Stall500)):
    stdsum = stdsum + Stall500_std[i]
    
Stall_500_std_avg = stdsum/len(Stall500)

Stall_500_var_std = Stall_500_std_avg**2

errorsum=0
avgsum=0
for i in range(len(Stall500)):
    errorsum= errorsum + (Stall500_std[i]-Stall_500_std_avg)**2
    avgsum = avgsum + (Stall500[i]-Stall_500_avg)**2
    
Stall_500_avg_error = np.sqrt(avgsum/len(Stall500))   
Stall_500_std_error = np.sqrt(errorsum/len(Stall500))


#%% Data resurrection 1300 nm bead

Res1300=AvgDict['N1300_resurrection_avg']

Res1300_std=AvgDict['N1300_resurrection_std']

AVbeg = int(200*FPS)
AVend = int(1000*FPS)

SS_data=AVend - AVbeg

Nsum=0
for i in range(AVbeg,AVend):
    Nsum = Nsum + Res1300[i]

Res_1300_avg=Nsum/SS_data

stdsum=0
for i in range(AVbeg,AVend):
    stdsum= stdsum + Res1300_std[i]

Res1300_std_avg=stdsum/SS_data

Res1300_var_avg = Res1300_std_avg**2

errorsum=0
avgsum=0
for i in range(AVbeg,AVend):
    errorsum= errorsum + (Res1300_std[i]-Res1300_std_avg)**2
    avgsum = avgsum + (Res1300[i]-Res_1300_avg)**2
    
Res_1300_avg_error = np.sqrt(avgsum/SS_data)   
Res_1300_std_error = np.sqrt(errorsum/SS_data)

#%% Data stall 1300 nm bead

Stall1300=AvgDict['N1300_after_avg']

Stall1300_std=AvgDict['N1300_after_std']

Nsum=0
for i in range(len(Stall1300)):
    Nsum = Nsum + Stall1300[i]
    
Stall_1300_avg = Nsum/len(Stall1300)

stdsum=0
for i in range(len(Stall1300)):
    stdsum = stdsum + Stall1300_std[i]
    
Stall_1300_std_avg = stdsum/len(Stall1300)

Stall_1300_var_std = Stall_1300_std_avg**2

errorsum=0
avgsum=0
for i in range(len(Stall1300)):
    errorsum= errorsum + (Stall1300_std[i]-Stall_1300_std_avg)**2
    avgsum = avgsum + (Stall1300[i]-Stall_1300_avg)**2
    
Stall_1300_avg_error = np.sqrt(avgsum/len(Stall1300))   
Stall_1300_std_error = np.sqrt(errorsum/len(Stall1300))

#%% Arrays for representations

Nmax=13

SS_Res = np.array([Res_300_avg, Res_500_avg, Res_1300_avg])/Nmax

SS_Stall = np.array([Stall_300_avg,Stall_500_avg,Stall_1300_avg])/Nmax

Std_Res = np.array([Res300_std_avg,Res500_std_avg,Res1300_std_avg])/Nmax

Std_Stall = np.array([Stall_300_std_avg,Stall_500_std_avg,Stall_1300_std_avg])/Nmax

Error_avg_Res = np.array([Res_300_avg_error,Res_500_avg_error,Res_1300_avg_error])/Nmax

Error_std_Res = np.array([Res_300_std_error,Res_500_std_error,Res_1300_std_error])/Nmax

Error_avg_Stall = np.array([Stall_300_avg_error,Stall_500_avg_error,Stall_1300_avg_error])/Nmax

Error_std_Stall = np.array([Stall_300_std_error,Stall_500_std_error,Stall_1300_std_error])/Nmax



Var_Res = np.array([Res300_var_avg,Res500_var_avg,Res1300_var_avg])/Nmax

Var_Stall = np.array([Stall_300_var_std,Stall_500_var_std,Stall_1300_var_std])/Nmax

beads = np.array([300,500,1300])

#%% Representation of the steady states

plt.title('Steady state vs bead diameter ($N_{max}=13$)')
plt.xlabel('Bead diameter / nm')
plt.ylabel('$<\phi>$')

plt.errorbar(beads,SS_Res,Std_Res,marker='^',label='Resurretion',capsize=3)
plt.errorbar(beads,SS_Stall,Std_Stall,marker='o',label='Stall',capsize=3)

plt.legend(loc='lower right')

plt.show()


#%%Data 300 before

T_300 = AvgDict['N300_before_avg']

T_std_300 = AvgDict['N300_before_std']

Nsum=0
for i in range(len(T_300)):
    Nsum = Nsum + T_300[i]
    
T_300_avg = Nsum/len(T_300)

stdsum=0
for i in range(len(T_300)):
    stdsum = stdsum + T_std_300[i]
    
T_300_std_avg = stdsum/len(T_300)

T_300_var_std = T_300_std_avg**2

errorsum=0
avgsum=0
for i in range(len(T_300)):
    errorsum= errorsum + (T_std_300[i]-T_300_std_avg)**2
    avgsum = avgsum + (T_300[i]-T_300_avg)**2
    
T_300_std_error = np.sqrt(errorsum/len(T_300))
T_300_avg_error = np.sqrt(avgsum/len(T_300))

#%% Data 500 before

T_500 = AvgDict['N500_before_avg']

T_std_500 = AvgDict['N500_before_std']


Nsum=0
for i in range(len(T_500)):
    Nsum = Nsum + T_500[i]
    
T_500_avg = Nsum/len(T_500)

stdsum=0
for i in range(len(T_500)):
    stdsum = stdsum + T_std_500[i]
    
T_500_std_avg = stdsum/len(T_500)

T_500_var_std = T_500_std_avg**2

errorsum=0
avgsum=0
for i in range(len(T_500)):
    errorsum= errorsum + (T_std_500[i]-T_500_std_avg)**2
    avgsum = avgsum + (T_500[i]-T_500_avg)**2
    
T_500_std_error = np.sqrt(errorsum/len(T_500))
T_500_avg_error = np.sqrt(avgsum/len(T_500))

#%% Data 1300 before

T_1300 = AvgDict['N1300_before_avg']

T_std_1300 = AvgDict['N1300_before_std']


Nsum=0
for i in range(len(T_1300)):
    Nsum = Nsum + T_1300[i]
    
T_1300_avg = Nsum/len(T_1300)

stdsum=0
for i in range(len(T_1300)):
    stdsum = stdsum + T_std_1300[i]
    
T_1300_std_avg = stdsum/len(T_1300)

T_1300_var_std = T_1300_std_avg**2

errorsum=0
avgsum=0
for i in range(len(T_1300)):
    errorsum= errorsum + (T_std_1300[i]-T_1300_std_avg)**2
    avgsum = avgsum + (T_1300[i]-T_1300_avg)**2
    
T_1300_std_error = np.sqrt(errorsum/len(T_1300))
T_1300_avg_error = np.sqrt(avgsum/len(T_1300))

#%% Arrays for representation of before data

SS_T = np.array([T_300_avg, T_500_avg, T_1300_avg])/Nmax

Std_T = np.array([T_300_std_avg,T_500_std_avg,T_1300_std_avg])/Nmax

Error_avg_T= np.array([T_300_avg_error,T_500_avg_error,T_1300_avg_error])/Nmax

Error_std_T = np.array([T_300_std_error,T_500_std_error,T_1300_std_error])/Nmax

#%% Representation of steady states vs variance

plt.title('Standard deviation vs steady state ($N_{max}=13$)')
plt.xlabel('$<\phi>$')
plt.ylabel('$(<\phi^2> - <\phi>^2)^{1/2}$')

plt.errorbar(SS_Res,Std_Res,yerr=Error_std_Res,xerr=Error_avg_Res,marker='^',label='Resurrection',capsize=3)
plt.errorbar(SS_Stall,Std_Stall,Error_std_Stall,Error_avg_Stall,marker='o',label='Stall',capsize=3)
plt.errorbar(SS_T,Std_T,Error_std_T,Error_avg_T,marker='^',label='Before',capsize=3)

plt.legend(loc='lower center')


#%%
files1=['/home/mariajose/Escritorio/Simulations/Without_ depletion/Glauber/Std_vs_Nss/sdt_vs_Nss_J=0_Res.dat'] 
data1=[]
files2=['/home/mariajose/Escritorio/Simulations/Without_ depletion/Glauber/Std_vs_Nss/sdt_vs_Nss_J=1_Res.dat'] 
data2=[]
files3=['/home/mariajose/Escritorio/Simulations/Without_ depletion/Glauber/Std_vs_Nss/sdt_vs_Nss_J=2_Res.dat']
data3=[]
files4=['/home/mariajose/Escritorio/Simulations/Without_ depletion/Glauber/Std_vs_Nss/sdt_vs_Nss_J=3_Res.dat'] 
data4=[]
files5=['/home/mariajose/Escritorio/Simulations/Without_ depletion/Glauber/Std_vs_Nss/sdt_vs_Nss_J=4_Res.dat'] 
data5=[]

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


J_0=data1[0]
J_1=data2[0]
J_2=data3[0]
J_3=data4[0]
J_4=data5[0]

#%%

np.savetxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/SS_Res.dat',(np.c_[SS_Res,Std_Res,Error_std_Res,Error_avg_Res]))

np.savetxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/SS_Stall.dat',(np.c_[SS_Stall,Std_Stall,Error_std_Stall,Error_avg_Stall]))

#%%
l=3

plt.xlabel('$<\phi>$',fontsize=14)
plt.ylabel('($<\phi^2> - <\phi>^2)^{1/2}$',fontsize=15)


plt.errorbar(SS_Res,Std_Res,Error_std_Res,Error_avg_Res,marker='^',label='Exp Res',capsize=3,linestyle='')
plt.errorbar(SS_Stall,Std_Stall,Error_std_Stall,Error_avg_Stall,marker='^',label='Exp Stall',capsize=3,linestyle='')
plt.errorbar(J_0[:,1],J_0[:,2],J_0[:,3],label='J=0 $k_BT$',marker='s',capsize=3,linestyle='')
plt.errorbar(J_1[:,1],J_1[:,2],J_1[:,3],label='J=1 $k_BT$',marker='s',capsize=3,linestyle='')
plt.errorbar(J_2[:,1],J_2[:,2],J_2[:,3],label='J=2 $k_BT$',marker='s',capsize=3,linestyle='')
#plt.errorbar(J_3[:,1],J_3[:,2],J_3[:,3],label='J=3 $k_BT$',marker='s',capsize=3,linestyle='')
#plt.errorbar(J_4[:,1],J_4[:,2],J_4[:,3],label='J=4 $k_BT$',marker='x',capsize=3)
plt.rcParams["figure.figsize"] = [8.0,6.0]

plt.legend(loc='upper left',prop={'size': 12})

#%%