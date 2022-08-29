#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 19:05:01 2022

@author: luke
"""

import numpy as np
import matplotlib.pyplot as plt
FPS=1000.

plt.rcParams["figure.figsize"] = [10.0,8.0]

#%% Ashley program

def LoadExpData():
    '''Loads the Experimental Data.  Change the pathnames below to match your comp'''
    
    D300 = np.load('/media/luke/Elements/Simulations/Data/NilsD300_V2.p',allow_pickle=True)
    D500 = np.load('/media/luke/Elements/Simulations/Data/NilsD500.p',allow_pickle=True)
    D1300 = np.load('/media/luke/Elements/Simulations/Data/NilsD1300.p',allow_pickle=True)
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

#%% DISTRIBUTION FOR STATES FOR 300nm BEAD

# General data from 300nm
t300 = Dict2DataMatrix(D300)
# Stators before stall
n300_before = t300[6]
# Stators after stall
n300_after = t300[7]
# Stators resurrection
n300_res = t300[8]

# Eliminate half of data to only have the equilibrium part
n300_after_eq = n300_after[:,int(len(n300_after[0,:])/2):len(n300_after[0,:])]
n300_res_eq = n300_res[:,int(len(n300_res[0,:])/2):len(n300_res[0,:])]

# Counts of each state
unique300_b, counts300_b = np.unique(n300_before,return_counts=True)
unique300_a, counts300_a = np.unique(n300_after_eq,return_counts=True)
unique300_r, counts300_r = np.unique(n300_res_eq,return_counts=True)

# Elimination of the nans
unique300_b_clean = unique300_b[0:len(unique300_b)-1]
counts300_b_clean = counts300_b[0:len(unique300_b_clean)]

unique300_a_clean = unique300_a[0:len(unique300_a)-1]
counts300_a_clean = counts300_a[0:len(unique300_a_clean)]

unique300_r_clean = unique300_r[0:len(unique300_r)-1]
counts300_r_clean = counts300_r[0:len(unique300_r_clean)]

# Matrix of frecuencies    
frecuencies300_b = np.asarray((unique300_b_clean,counts300_b_clean)).T    
total300_b = sum(frecuencies300_b[:,1])

frecuencies300_a = np.asarray((unique300_a_clean,counts300_a_clean)).T    
total300_a = sum(frecuencies300_a[:,1])

frecuencies300_r = np.asarray((unique300_r_clean,counts300_r_clean)).T    
total300_r = sum(frecuencies300_r[:,1])

# Histogram 
plt.title('Distribution states 300nm bead',fontsize=16)
plt.xlabel('N',fontsize=16)
plt.ylabel('Probability',fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.grid(alpha=0.5)
plt.bar(frecuencies300_b[:,0],frecuencies300_b[:,1]/total300_b,width=1,ls='solid',lw=2, edgecolor='royalblue',hatch='/',fc='None',label='Before stall')
plt.bar(frecuencies300_a[:,0],frecuencies300_a[:,1]/total300_a,width=1,ls='solid',lw=2, edgecolor='orange',hatch='.',fc='None',label='After stall')
plt.bar(frecuencies300_r[:,0],frecuencies300_r[:,1]/total300_r,width=1,ls='solid',lw=2, edgecolor='green',hatch='-',fc='None',label='Resurrection')

plt.legend(fontsize=16)

#%% DISTRIBUTION FOR STATES FOR 500nm BEAD

# General data from 500nm
t500 = Dict2DataMatrix(D500)
# Stators before stall
n500_before = t500[6]
# Stators after stall
n500_after = t500[7]
# Stators resurrection
n500_res = t500[8]

# Eliminate half of data to only have the equilibrium part
n500_after_eq = n500_after[:,int(len(n500_after[0,:])/2):len(n500_after[0,:])]
n500_res_eq = n500_res[:,int(len(n500_res[0,:])/2):len(n500_res[0,:])]

# Counts of each state
unique500, counts500 = np.unique(n500_before,return_counts=True)
unique500_a, counts500_a = np.unique(n500_after_eq,return_counts=True)
unique500_r, counts500_r = np.unique(n500_res_eq,return_counts=True)

# Elimination of the nans
unique500_clean = unique500[0:len(unique500)-1]
counts500_clean = counts500[0:len(unique500_clean)]

unique500_a_clean = unique500_a[0:len(unique500_a)-1]
counts500_a_clean = counts500_a[0:len(unique500_a_clean)]

unique500_r_clean = unique500_r[0:len(unique500_r)-1]
counts500_r_clean = counts500_r[0:len(unique500_r_clean)]


# Matrix of frecuencies     
frecuencies500 = np.asarray((unique500_clean,counts500_clean)).T    
total500 = sum(frecuencies500[:,1])

frecuencies500_a = np.asarray((unique500_a_clean,counts500_a_clean)).T    
total500_a = sum(frecuencies500_a[:,1])

frecuencies500_r = np.asarray((unique500_r_clean,counts500_r_clean)).T    
total500_r = sum(frecuencies500_r[:,1])

# Histogram
plt.title('Distribution states 500nm bead',fontsize=16)
plt.xlabel('N',fontsize=16)
plt.ylabel('Probability',fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.grid(alpha=0.5) 
plt.bar(frecuencies500[:,0],frecuencies500[:,1]/total500,width=1,ls='solid',lw=2, edgecolor='royalblue',hatch='/',fc='None',label='Before stall')
plt.bar(frecuencies500_a[:,0],frecuencies500_a[:,1]/total500_a,width=1,ls='solid',lw=2, edgecolor='orange',hatch='.',fc='None',label='After stall')
plt.bar(frecuencies500_r[:,0],frecuencies500_r[:,1]/total500_r,width=1,ls='solid',lw=2, edgecolor='green',hatch='-',fc='None',label='Resurrection')

plt.legend(fontsize=16)

#%% #%% DISTRIBUTION FOR STATES FOR 1300nm BEAD

# General data from 1300nm
t1300 = Dict2DataMatrix(D1300)
# Stators before stall
n1300_before = t1300[6]
# Stators after stall
n1300_after = t1300[7]
# Stators resurrection
n1300_res = t1300[8]

# Eliminate half of data to only have the equilibrium part
n1300_after_eq = n1300_after[:,int(len(n1300_after[0,:])/2):len(n1300_after[0,:])]
n1300_res_eq = n1300_res[:,int(len(n1300_res[0,:])/2):len(n1300_res[0,:])]

# Counts of each state
unique1300, counts1300 = np.unique(n1300_before,return_counts=True)
unique1300_a, counts1300_a = np.unique(n1300_after_eq,return_counts=True)
unique1300_r, counts1300_r = np.unique(n1300_res_eq,return_counts=True)

# Elimination of the nans
unique1300_clean = unique1300[0:len(unique1300)-1]
counts1300_clean = counts1300[0:len(unique1300_clean)]

unique1300_a_clean = unique1300_a[0:len(unique1300_a)-1]
counts1300_a_clean = counts1300_a[0:len(unique1300_a_clean)]

unique1300_r_clean = unique1300_r[0:len(unique1300_r)-1]
counts1300_r_clean = counts1300_r[0:len(unique1300_r_clean)]

# Matrix of frecuencies       
frecuencies1300 = np.asarray((unique1300_clean,counts1300_clean)).T    
total1300 = sum(frecuencies1300[:,1])

frecuencies1300_a = np.asarray((unique1300_a_clean,counts1300_a_clean)).T    
total1300_a = sum(frecuencies1300_a[:,1])

frecuencies1300_r = np.asarray((unique1300_r_clean,counts1300_r_clean)).T    
total1300_r = sum(frecuencies1300_r[:,1])

# Histogram
plt.title('Distribution states 1300nm bead',fontsize=16)
plt.xlabel('N',fontsize=16)
plt.ylabel('Probability',fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.grid(alpha=0.5) 
plt.bar(frecuencies1300[:,0],frecuencies1300[:,1]/total1300,width=1,ls='solid',lw=2, edgecolor='royalblue',hatch='/',fc='None',label='Before stall')
plt.bar(frecuencies1300_a[:,0],frecuencies1300_a[:,1]/total1300_a,width=1,ls='solid',lw=2, edgecolor='orange',hatch='.',fc='None',label='After stall')
plt.bar(frecuencies1300_r[:,0],frecuencies1300_r[:,1]/total1300_r,width=1,ls='solid',lw=2, edgecolor='green',hatch='-',fc='None',label='Resurrection')

plt.legend(fontsize=16)
