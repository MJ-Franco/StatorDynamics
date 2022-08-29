#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 17:13:17 2021

@author: mariajose
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

#%% Useful variables and functions

nmax=13
l=5000
#Function used for fits
def N(x,tau,nss):
    return nss + (n0-nss)*np.exp(-x/tau) 

#%% Files for resurrection

files1=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Res_J=0.00  _4.dat'] #All resurrection simulations
data1=[]
files11=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Res_J=0.25  _4.dat'] #All resurrection simulations
data11=[]
files12=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Res_J=0.50  _4.dat'] #All resurrection simulations
data12=[]
files13=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Res_J=0.75  _4.dat'] #All resurrection simulations
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

files2=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Res_J=1.00  _4.dat'] #Average of resurrection simulations
data2=[]
files21=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Res_J=1.25  _4.dat'] #Average of resurrection simulations
data21=[]
files22=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Res_J=1.50  _4.dat'] #Average of resurrection simulations
data22=[]
files23=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Res_J=1.75  _4.dat'] #Average of resurrection simulations
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

files3=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Res_J=2.00  _4.dat'] #Reservoir in each simulation
data3=[]
files31=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Res_J=2.25  _4.dat'] #Reservoir in each simulation
data31=[]
files32=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Res_J=2.50  _4.dat'] #Reservoir in each simulation
data32=[]
files33=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Res_J=2.75  _4.dat'] #Reservoir in each simulation
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

files4=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Res_J=3.00  _4.dat'] #Average of reservoir in resurrection simulations
data4=[]
files41=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Res_J=3.25  _4.dat'] #Average of reservoir in resurrection simulations
data41=[]
files42=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Res_J=3.50  _4.dat'] #Average of reservoir in resurrection simulations
data42=[]
files43=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Res_J=3.75  _4.dat'] #Average of reservoir in resurrection simulations
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

files5=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Res_J=4.00  _4.dat'] #Average of reservoir in resurrection simulations
data5=[]
files51=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Res_J=4.25  _4.dat'] #Average of reservoir in resurrection simulations
data51=[]
files52=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Res_J=4.50  _4.dat'] #Average of reservoir in resurrection simulations
data52=[]
files53=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Res_J=4.75  _4.dat'] #Average of reservoir in resurrection simulations
data53=[]
files54=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Res_J=5.00  _4.dat'] #Average of reservoir in resurrection simulations
data54=[]

for data_file in files5:
    data5.append(np.loadtxt(data_file))
for data_file in files51:
    data51.append(np.loadtxt(data_file))
for data_file in files52:
    data52.append(np.loadtxt(data_file))
for data_file in files53:
    data53.append(np.loadtxt(data_file))
for data_file in files54:
    data54.append(np.loadtxt(data_file))

Res_J_4=data5[0]
Res_J_42=data51[0]
Res_J_45=data52[0]
Res_J_47=data53[0]
Res_J_5=data54[0]


#%% Files for stall

files6=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Stall_J=0.00  _4_dist.dat'] #All resurrection simulations
data6=[]
files61=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Stall_J=0.25  _4_dist.dat'] #All resurrection simulations
data61=[]
files62=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Stall_J=0.50  _4_dist.dat'] #All resurrection simulations
data62=[]
files63=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Stall_J=0.75  _4_dist.dat'] #All resurrection simulations
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

files7=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Stall_J=1.00  _4_dist.dat'] #Average of resurrection simulations
data7=[]
files71=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Stall_J=1.25  _4_dist.dat'] #Average of resurrection simulations
data71=[]
files72=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Stall_J=1.50  _4_dist.dat'] #Average of resurrection simulations
data72=[]
files73=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Stall_J=1.75  _4_dist.dat'] #Average of resurrection simulations
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

files8=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Stall_J=2.00  _4_dist.dat'] #Reservoir in each simulation
data8=[]
files81=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Stall_J=2.25  _4_dist.dat'] #Reservoir in each simulation
data81=[]
files82=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Stall_J=2.50  _4_dist.dat'] #Reservoir in each simulation
data82=[]
files83=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Stall_J=2.75  _4_dist.dat'] #Reservoir in each simulation
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

files9=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Stall_J=3.00  _4_dist.dat'] #Average of reservoir in resurrection simulations
data9=[]
files91=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Stall_J=3.25  _4_dist.dat'] #Average of reservoir in resurrection simulations
data91=[]
files92=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Stall_J=3.50  _4_dist.dat'] #Average of reservoir in resurrection simulations
data92=[]
files93=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Stall_J=3.75  _4_dist.dat'] #Average of reservoir in resurrection simulations
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

files10=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Stall_J=4.00  _4_dist.dat'] #Average of reservoir in resurrection simulations
data10=[]
files101=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Stall_J=4.25  _4_dist.dat'] #Average of reservoir in resurrection simulations
data101=[]
files102=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Stall_J=4.50  _4_dist.dat'] #Average of reservoir in resurrection simulations
data102=[]
files103=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Stall_J=4.75  _4_dist.dat'] #Average of reservoir in resurrection simulations
data103=[]
files104=['/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/Traces_Stall_J=5.00  _4_dist.dat'] #Average of reservoir in resurrection simulations
data104=[]
    
for data_file in files10:
    data10.append(np.loadtxt(data_file))
for data_file in files101:
    data101.append(np.loadtxt(data_file))
for data_file in files102:
    data102.append(np.loadtxt(data_file))
for data_file in files103:
    data103.append(np.loadtxt(data_file))
for data_file in files104:
    data104.append(np.loadtxt(data_file))


Stall_J_4=data10[0]
Stall_J_42=data101[0]
Stall_J_45=data102[0]
Stall_J_47=data103[0]
Stall_J_5=data104[0]

#%% Files for resurrection With_depletion

files100=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Res_J=0.00  .dat'] #All resurrection simulations
data100=[]
files110=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Res_J=0.25  .dat'] #All resurrection simulations
data110=[]
files120=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Res_J=0.50  .dat'] #All resurrection simulations
data120=[]
files130=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Res_J=0.75  .dat'] #All resurrection simulations
data130=[]

for data_file in files100:
    data100.append(np.loadtxt(data_file))
    
for data_file in files110:
    data110.append(np.loadtxt(data_file))
    
for data_file in files120:
    data120.append(np.loadtxt(data_file))

for data_file in files130:
    data130.append(np.loadtxt(data_file))
    
Res_J_0_DEP=data100[0]
Res_J_02_DEP=data110[0]
Res_J_05_DEP=data120[0]
Res_J_07_DEP=data120[0]



files200=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Res_J=1.00  .dat'] #Average of resurrection simulations
data200=[]
files210=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Res_J=1.25  .dat'] #Average of resurrection simulations
data210=[]
files220=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Res_J=1.50  .dat'] #Average of resurrection simulations
data220=[]
files230=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Res_J=1.75  .dat'] #Average of resurrection simulations
data230=[]

for data_file in files200:
    data200.append(np.loadtxt(data_file))
    
for data_file in files210:
    data210.append(np.loadtxt(data_file))
    
for data_file in files220:
    data220.append(np.loadtxt(data_file))

for data_file in files230:
    data230.append(np.loadtxt(data_file))
    
Res_J_1_DEP=data200[0]
Res_J_12_DEP=data210[0]
Res_J_15_DEP=data220[0]
Res_J_17_DEP=data230[0]



files300=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Res_J=2.00  .dat'] #Reservoir in each simulation
data300=[]
files310=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Res_J=2.25  .dat'] #Reservoir in each simulation
data310=[]
files320=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Res_J=2.50  .dat'] #Reservoir in each simulation
data320=[]
files330=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Res_J=2.75  .dat'] #Reservoir in each simulation
data330=[]

for data_file in files300:
    data300.append(np.loadtxt(data_file)) 
    
for data_file in files310:
    data310.append(np.loadtxt(data_file)) 
    
for data_file in files320:
    data320.append(np.loadtxt(data_file))
    
for data_file in files330:
    data330.append(np.loadtxt(data_file)) 
    
Res_J_2_DEP=data300[0]
Res_J_22_DEP=data310[0]
Res_J_25_DEP=data320[0]
Res_J_27_DEP=data330[0]



files400=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Res_J=3.00  .dat'] #Average of reservoir in resurrection simulations
data400=[]
files410=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Res_J=3.25  .dat'] #Average of reservoir in resurrection simulations
data410=[]
files420=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Res_J=3.50  .dat'] #Average of reservoir in resurrection simulations
data420=[]
files430=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Res_J=3.75  .dat'] #Average of reservoir in resurrection simulations
data430=[]

for data_file in files400:
    data400.append(np.loadtxt(data_file))
    
for data_file in files410:
    data410.append(np.loadtxt(data_file))
    
for data_file in files420:
    data420.append(np.loadtxt(data_file))
    
for data_file in files430:
    data430.append(np.loadtxt(data_file))
    
Res_J_3_DEP=data400[0]
Res_J_32_DEP=data410[0]
Res_J_35_DEP=data420[0]
Res_J_37_DEP=data430[0]

files500=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Res_J=4.00  .dat'] #Average of reservoir in resurrection simulations
data500=[]
files510=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Res_J=4.25  .dat'] #Average of reservoir in resurrection simulations
data510=[]
files520=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Res_J=4.50  .dat'] #Average of reservoir in resurrection simulations
data520=[]
files530=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Res_J=4.75  .dat'] #Average of reservoir in resurrection simulations
data530=[]
files600=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Res_J=5.00  .dat'] #Average of reservoir in resurrection simulations
data600=[]

for data_file in files500:
    data500.append(np.loadtxt(data_file))
    
for data_file in files510:
    data510.append(np.loadtxt(data_file))

for data_file in files520:
    data520.append(np.loadtxt(data_file))

for data_file in files530:
    data530.append(np.loadtxt(data_file))

for data_file in files600:
    data600.append(np.loadtxt(data_file))

Res_J_4_DEP=data500[0]
Res_J_42_DEP=data510[0]
Res_J_45_DEP=data520[0]
Res_J_47_DEP=data530[0]
Res_J_5_DEP=data600[0]


#%% Files for stall With_depletion

files600=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Stall_J=0.00  .dat'] #All resurrection simulations
data600=[]
files610=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Stall_J=0.25  .dat'] #All resurrection simulations
data610=[]
files620=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Stall_J=0.50  .dat'] #All resurrection simulations
data620=[]
files630=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Stall_J=0.75  .dat'] #All resurrection simulations
data630=[]

for data_file in files600:
    data600.append(np.loadtxt(data_file))
    
for data_file in files610:
    data610.append(np.loadtxt(data_file))
    
for data_file in files620:
    data620.append(np.loadtxt(data_file))
    
for data_file in files630:
    data630.append(np.loadtxt(data_file))
    
    
Stall_J_0_DEP=data600[0]
Stall_J_02_DEP=data610[0]
Stall_J_05_DEP=data620[0]
Stall_J_07_DEP=data630[0]



files700=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Stall_J=1.00  .dat'] #Average of resurrection simulations
data700=[]
files710=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Stall_J=1.25  .dat'] #Average of resurrection simulations
data710=[]
files720=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Stall_J=1.50  .dat'] #Average of resurrection simulations
data720=[]
files730=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Stall_J=1.75  .dat'] #Average of resurrection simulations
data730=[]

for data_file in files700:
    data700.append(np.loadtxt(data_file))
    
for data_file in files710:
    data710.append(np.loadtxt(data_file))
    
for data_file in files720:
    data720.append(np.loadtxt(data_file))

for data_file in files730:
    data730.append(np.loadtxt(data_file))
    
    
Stall_J_1_DEP=data700[0]
Stall_J_12_DEP=data710[0]
Stall_J_15_DEP=data720[0]
Stall_J_17_DEP=data730[0]



files800=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Stall_J=2.00  .dat'] #Reservoir in each simulation
data800=[]
files810=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Stall_J=2.25  .dat'] #Reservoir in each simulation
data810=[]
files820=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Stall_J=2.50  .dat'] #Reservoir in each simulation
data820=[]
files830=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Stall_J=2.75  .dat'] #Reservoir in each simulation
data830=[]

for data_file in files800:
    data800.append(np.loadtxt(data_file))
    
for data_file in files810:
    data810.append(np.loadtxt(data_file))
    
for data_file in files820:
    data820.append(np.loadtxt(data_file))
    
for data_file in files830:
    data830.append(np.loadtxt(data_file))
    
    
Stall_J_2_DEP=data800[0]
Stall_J_22_DEP=data810[0]
Stall_J_25_DEP=data820[0]
Stall_J_27_DEP=data830[0]



files900=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Stall_J=3.00  .dat'] #Average of reservoir in resurrection simulations
data900=[]
files910=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Stall_J=3.25  .dat'] #Average of reservoir in resurrection simulations
data910=[]
files920=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Stall_J=3.50  .dat'] #Average of reservoir in resurrection simulations
data920=[]
files930=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Stall_J=3.75  .dat'] #Average of reservoir in resurrection simulations
data930=[]

for data_file in files900:
    data900.append(np.loadtxt(data_file))
    
for data_file in files910:
    data910.append(np.loadtxt(data_file))
    
for data_file in files920:
    data920.append(np.loadtxt(data_file))
    
for data_file in files930:
    data930.append(np.loadtxt(data_file))
    
Stall_J_3_DEP=data900[0]
Stall_J_32_DEP=data910[0]
Stall_J_35_DEP=data920[0]
Stall_J_37_DEP=data930[0]

files1000=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Stall_J=4.00  .dat'] #Average of reservoir in resurrection simulations
data1000=[]
files1100=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Stall_J=4.25  .dat'] #Average of reservoir in resurrection simulations
data1100=[]
files1200=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Stall_J=4.50  .dat'] #Average of reservoir in resurrection simulations
data1200=[]
files1300=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Stall_J=4.75  .dat'] #Average of reservoir in resurrection simulations
data1300=[]
files2000=['/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/Traces_Stall_J=5.00  .dat'] #Average of reservoir in resurrection simulations
data2000=[]
    
for data_file in files1000:
    data1000.append(np.loadtxt(data_file))

for data_file in files1100:
    data1100.append(np.loadtxt(data_file))

for data_file in files1200:
    data1200.append(np.loadtxt(data_file))

for data_file in files1300:
    data1300.append(np.loadtxt(data_file))

for data_file in files2000:
    data2000.append(np.loadtxt(data_file))

Stall_J_4_DEP=data1000[0]
Stall_J_42_DEP=data1100[0]
Stall_J_45_DEP=data1200[0]
Stall_J_47_DEP=data1300[0]
Stall_J_5_DEP=data2000[0]

#%% Fit for resurrection Without _depletion
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

poRes42, pcRes42 = curve_fit(N,Res_J_42[0:l,0],Res_J_42[0:l,1]*nmax)
tau_Res42=poRes42[0]
err_Res42=np.sqrt(pcRes42[0,0])

poRes45, pcRes45 = curve_fit(N,Res_J_45[0:l,0],Res_J_45[0:l,1]*nmax)
tau_Res45=poRes45[0]
err_Res45=np.sqrt(pcRes45[0,0])

poRes47, pcRes47 = curve_fit(N,Res_J_47[0:l,0],Res_J_47[0:l,1]*nmax)
tau_Res47=poRes47[0]
err_Res47=np.sqrt(pcRes47[0,0])

poRes5, pcRes5 = curve_fit(N,Res_J_5[0:l,0],Res_J_5[0:l,1]*nmax)
tau_Res5=poRes5[0]
err_Res5=np.sqrt(pcRes5[0,0])

#%%

plt.title('Resurrection Without _depletion ($N_{ss}=8$)',fontsize=17)
plt.xlabel('Simulation time',fontsize=14)
plt.ylabel('Stator number',fontsize=14)
plt.rcParams["figure.figsize"] = [8.0,6.0]

plt.plot(Res_J_0[0:l,0],Res_J_0[0:l,1]*nmax,label=r'J=0 $k_BT$, $\tau$={:1.2f}'.format(tau_Res0))
plt.plot(Res_J_1[0:l,0],Res_J_1[0:l,1]*nmax,label=r'J=1 $k_BT$, $\tau$={:1.2f}'.format(tau_Res1))
plt.plot(Res_J_2[0:l,0],Res_J_2[0:l,1]*nmax,label=r'J=2 $k_BT$, $\tau$={:1.2f}'.format(tau_Res2))
plt.plot(Res_J_3[0:l,0],Res_J_3[0:l,1]*nmax,label=r'J=3 $k_BT$, $\tau$={:1.2f}'.format(tau_Res3))
plt.plot(Res_J_4[0:l,0],Res_J_4[0:l,1]*nmax,label=r'J=4 $k_BT$, $\tau$={:1.2f}'.format(tau_Res4))
plt.plot(Res_J_5[0:l,0],Res_J_5[0:l,1]*nmax,label=r'J=5 $k_BT$, $\tau$={:1.2f}'.format(tau_Res5))      

plt.legend(loc='lower right', prop={'size': 15})

#%% Fit for resurrection With_depletion

n0=0

poRes0Dep, pcRes0Dep = curve_fit(N,Res_J_0_DEP[0:l,0],Res_J_0_DEP[0:l,1]*nmax)
tau_Res0Dep=poRes0Dep[0]
err_Res0Dep=np.sqrt(pcRes0Dep[0,0])

poRes02Dep, pcRes02Dep = curve_fit(N,Res_J_02_DEP[0:l,0],Res_J_02_DEP[0:l,1]*nmax)
tau_Res02Dep=poRes02Dep[0]
err_Res02Dep=np.sqrt(pcRes02Dep[0,0])

poRes05Dep, pcRes05Dep = curve_fit(N,Res_J_05_DEP[0:l,0],Res_J_05_DEP[0:l,1]*nmax)
tau_Res05Dep=poRes05Dep[0]
err_Res05Dep=np.sqrt(pcRes05Dep[0,0])

poRes07Dep, pcRes07Dep = curve_fit(N,Res_J_07_DEP[0:l,0],Res_J_07_DEP[0:l,1]*nmax)
tau_Res07Dep=poRes07Dep[0]
err_Res07Dep=np.sqrt(pcRes07Dep[0,0])



poRes1Dep, pcRes1Dep = curve_fit(N,Res_J_1_DEP[0:l,0],Res_J_1_DEP[0:l,1]*nmax)
tau_Res1Dep=poRes1Dep[0]
err_Res1Dep=np.sqrt(pcRes1Dep[0,0])

poRes12Dep, pcRes12Dep = curve_fit(N,Res_J_12_DEP[0:l,0],Res_J_12_DEP[0:l,1]*nmax)
tau_Res12Dep=poRes12Dep[0]
err_Res12Dep=np.sqrt(pcRes12Dep[0,0])

poRes15Dep, pcRes15Dep = curve_fit(N,Res_J_15_DEP[0:l,0],Res_J_15_DEP[0:l,1]*nmax)
tau_Res15Dep=poRes15Dep[0]
err_Res15Dep=np.sqrt(pcRes15Dep[0,0])

poRes17Dep, pcRes17Dep = curve_fit(N,Res_J_17_DEP[0:l,0],Res_J_17_DEP[0:l,1]*nmax)
tau_Res17Dep=poRes17Dep[0]
err_Res17Dep=np.sqrt(pcRes17Dep[0,0])



poRes2Dep, pcRes2Dep = curve_fit(N,Res_J_2_DEP[0:l,0],Res_J_2_DEP[0:l,1]*nmax)
tau_Res2Dep=poRes2Dep[0]
err_Res2Dep=np.sqrt(pcRes2Dep[0,0])

poRes22Dep, pcRes22Dep = curve_fit(N,Res_J_22_DEP[0:l,0],Res_J_22_DEP[0:l,1]*nmax)
tau_Res22Dep=poRes22Dep[0]
err_Res22Dep=np.sqrt(pcRes22Dep[0,0])

poRes25Dep, pcRes25Dep = curve_fit(N,Res_J_25_DEP[0:l,0],Res_J_25_DEP[0:l,1]*nmax)
tau_Res25Dep=poRes25Dep[0]
err_Res25Dep=np.sqrt(pcRes25Dep[0,0])

poRes27Dep, pcRes27Dep = curve_fit(N,Res_J_27_DEP[0:l,0],Res_J_27_DEP[0:l,1]*nmax)
tau_Res27Dep=poRes27Dep[0]
err_Res27Dep=np.sqrt(pcRes27Dep[0,0])



poRes3Dep, pcRes3Dep = curve_fit(N,Res_J_3_DEP[0:l,0],Res_J_3_DEP[0:l,1]*nmax)
tau_Res3Dep=poRes3Dep[0]
err_Res3Dep=np.sqrt(pcRes3Dep[0,0])

poRes32Dep, pcRes32Dep = curve_fit(N,Res_J_32_DEP[0:l,0],Res_J_32_DEP[0:l,1]*nmax)
tau_Res32Dep=poRes32Dep[0]
err_Res32Dep=np.sqrt(pcRes32Dep[0,0])

poRes35Dep, pcRes35Dep = curve_fit(N,Res_J_35_DEP[0:l,0],Res_J_35_DEP[0:l,1]*nmax)
tau_Res35Dep=poRes35Dep[0]
err_Res35Dep=np.sqrt(pcRes35Dep[0,0])

poRes37Dep, pcRes37Dep = curve_fit(N,Res_J_37_DEP[0:l,0],Res_J_37_DEP[0:l,1]*nmax)
tau_Res37Dep=poRes37Dep[0]
err_Res37Dep=np.sqrt(pcRes37Dep[0,0])



poRes4Dep, pcRes4Dep = curve_fit(N,Res_J_4_DEP[0:l,0],Res_J_4_DEP[0:l,1]*nmax)
tau_Res4Dep=poRes4Dep[0]
err_Res4Dep=np.sqrt(pcRes4Dep[0,0])

poRes42Dep, pcRes42Dep = curve_fit(N,Res_J_42_DEP[0:l,0],Res_J_42_DEP[0:l,1]*nmax)
tau_Res42Dep=poRes42Dep[0]
err_Res42Dep=np.sqrt(pcRes42Dep[0,0])

poRes45Dep, pcRes45Dep = curve_fit(N,Res_J_45_DEP[0:l,0],Res_J_45_DEP[0:l,1]*nmax)
tau_Res45Dep=poRes45Dep[0]
err_Res45Dep=np.sqrt(pcRes45Dep[0,0])

poRes47Dep, pcRes47Dep = curve_fit(N,Res_J_47_DEP[0:l,0],Res_J_47_DEP[0:l,1]*nmax)
tau_Res47Dep=poRes47Dep[0]
err_Res47Dep=np.sqrt(pcRes47Dep[0,0])

poRes5Dep, pcRes5Dep = curve_fit(N,Res_J_5_DEP[0:l,0],Res_J_5_DEP[0:l,1]*nmax)
tau_Res5Dep=poRes5Dep[0]
err_Res5Dep=np.sqrt(pcRes5Dep[0,0])

#%%

plt.title('Resurrection With_depletion ($N_{ss}=8$)',fontsize=17)
plt.xlabel('Simulation time',fontsize=14)
plt.ylabel('Stator number',fontsize=14)
plt.rcParams["figure.figsize"] = [8.0,6.0]

plt.plot(Res_J_0_DEP[0:l,0],Res_J_0_DEP[0:l,1]*nmax,label=r'J=0 $k_BT$, $\tau$={:1.2f}'.format(tau_Res0Dep))
plt.plot(Res_J_1_DEP[0:l,0],Res_J_1_DEP[0:l,1]*nmax,label=r'J=1 $k_BT$, $\tau$={:1.2f}'.format(tau_Res1Dep))
plt.plot(Res_J_2_DEP[0:l,0],Res_J_2_DEP[0:l,1]*nmax,label=r'J=2 $k_BT$, $\tau$={:1.2f}'.format(tau_Res2Dep))
plt.plot(Res_J_3_DEP[0:l,0],Res_J_3_DEP[0:l,1]*nmax,label=r'J=3 $k_BT$, $\tau$={:1.2f}'.format(tau_Res3Dep))
plt.plot(Res_J_4_DEP[0:l,0],Res_J_4_DEP[0:l,1]*nmax,label=r'J=4 $k_BT$, $\tau$={:1.2f}'.format(tau_Res4Dep))
plt.plot(Res_J_5_DEP[0:l,0],Res_J_5_DEP[0:l,1]*nmax,label=r'J=5 $k_BT$, $\tau$={:1.2f}'.format(tau_Res5Dep))       

plt.legend(loc='lower right', prop={'size': 15})


#%% Fit for stall Without _depletion

n0=9

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

poStall42, pcStall42 = curve_fit(N,Stall_J_42[0:l,0],Stall_J_42[0:l,1]*nmax)
tau_Stall42=poStall42[0]
err_Stall42=np.sqrt(pcStall42[0,0])

poStall45, pcStall45 = curve_fit(N,Stall_J_45[0:l,0],Stall_J_45[0:l,1]*nmax)
tau_Stall45=poStall45[0]
err_Stall45=np.sqrt(pcStall45[0,0])

poStall47, pcStall47 = curve_fit(N,Stall_J_47[0:l,0],Stall_J_47[0:l,1]*nmax)
tau_Stall47=poStall47[0]
err_Stall47=np.sqrt(pcStall47[0,0])

poStall5, pcStall5 = curve_fit(N,Stall_J_5[0:l,0],Stall_J_5[0:l,1]*nmax)
tau_Stall5=poStall5[0]
err_Stall5=np.sqrt(pcStall5[0,0])

#%%

plt.title('Stall Without _depletion ($N_{ss}=8$)',fontsize=17)
plt.xlabel('Simulation time',fontsize=14)
plt.ylabel('Stator number',fontsize=14)
plt.rcParams["figure.figsize"] = [8.0,6.0]

plt.plot(Stall_J_0[0:l,0],Stall_J_0[0:l,1]*nmax,label=r'J=0 $k_BT$, $\tau$={:1.2f}'.format(tau_Stall0))
plt.plot(Stall_J_1[0:l,0],Stall_J_1[0:l,1]*nmax,label=r'J=1 $k_BT$, $\tau$={:1.2f}'.format(tau_Stall1))
plt.plot(Stall_J_2[0:l,0],Stall_J_2[0:l,1]*nmax,label=r'J=2 $k_BT$, $\tau$={:1.2f}'.format(tau_Stall2))
plt.plot(Stall_J_3[0:l,0],Stall_J_3[0:l,1]*nmax,label=r'J=3 $k_BT$, $\tau$={:1.2f}'.format(tau_Stall3))
plt.plot(Stall_J_4[0:l,0],Stall_J_4[0:l,1]*nmax,label=r'J=4 $k_BT$, $\tau$={:1.2f}'.format(tau_Stall4))
plt.plot(Stall_J_5[0:l,0],Stall_J_5[0:l,1]*nmax,label=r'J=5 $k_BT$, $\tau$={:1.2f}'.format(tau_Stall5))

#plt.hlines(8.05,0,200, linestyles='dashed',alpha=0.4)    

plt.legend(loc='right', prop={'size': 15})

#%% Fit for stall With_depletion

n0=9

poStall0Dep, pcStall0Dep = curve_fit(N,Stall_J_0_DEP[0:l,0],Stall_J_0_DEP[0:l,1]*nmax)
tau_Stall0Dep=poStall0Dep[0]
err_Stall0Dep=np.sqrt(pcStall0Dep[0,0])

poStall02Dep, pcStall02Dep= curve_fit(N,Stall_J_02_DEP[0:l,0],Stall_J_02_DEP[0:l,1]*nmax)
tau_Stall02Dep=poStall02Dep[0]
err_Stall02Dep=np.sqrt(pcStall02Dep[0,0])

poStall05Dep, pcStall05Dep = curve_fit(N,Stall_J_05_DEP[0:l,0],Stall_J_05_DEP[0:l,1]*nmax)
tau_Stall05Dep=poStall05Dep[0]
err_Stall05Dep=np.sqrt(pcStall05Dep[0,0])

poStall07Dep, pcStall07Dep = curve_fit(N,Stall_J_07_DEP[0:l,0],Stall_J_07_DEP[0:l,1]*nmax)
tau_Stall07Dep=poStall07Dep[0]
err_Stall07Dep=np.sqrt(pcStall07Dep[0,0])



poStall1Dep, pcStall1Dep = curve_fit(N,Stall_J_1_DEP[0:l,0],Stall_J_1_DEP[0:l,1]*nmax)
tau_Stall1Dep=poStall1Dep[0]
err_Stall1Dep=np.sqrt(pcStall1Dep[0,0])

poStall12Dep, pcStall12Dep = curve_fit(N,Stall_J_12_DEP[0:l,0],Stall_J_12_DEP[0:l,1]*nmax)
tau_Stall12Dep=poStall12Dep[0]
err_Stall12Dep=np.sqrt(pcStall12Dep[0,0])

poStall15Dep, pcStall15Dep = curve_fit(N,Stall_J_15_DEP[0:l,0],Stall_J_15_DEP[0:l,1]*nmax)
tau_Stall15Dep=poStall15Dep[0]
err_Stall15Dep=np.sqrt(pcStall15Dep[0,0])

poStall17Dep, pcStall17Dep = curve_fit(N,Stall_J_17_DEP[0:l,0],Stall_J_17_DEP[0:l,1]*nmax)
tau_Stall17Dep=poStall17Dep[0]
err_Stall17Dep=np.sqrt(pcStall17Dep[0,0])



poStall2Dep, pcStall2Dep = curve_fit(N,Stall_J_2_DEP[0:l,0],Stall_J_2_DEP[0:l,1]*nmax)
tau_Stall2Dep=poStall2Dep[0]
err_Stall2Dep=np.sqrt(pcStall2Dep[0,0])

poStall22Dep, pcStall22Dep = curve_fit(N,Stall_J_22_DEP[0:l,0],Stall_J_22_DEP[0:l,1]*nmax)
tau_Stall22Dep=poStall22Dep[0]
err_Stall22Dep=np.sqrt(pcStall22Dep[0,0])

poStall25Dep, pcStall25Dep = curve_fit(N,Stall_J_25_DEP[0:l,0],Stall_J_25_DEP[0:l,1]*nmax)
tau_Stall25Dep=poStall25Dep[0]
err_Stall25Dep=np.sqrt(pcStall25Dep[0,0])

poStall27Dep, pcStall27Dep = curve_fit(N,Stall_J_27_DEP[0:l,0],Stall_J_27_DEP[0:l,1]*nmax)
tau_Stall27Dep=poStall27Dep[0]
err_Stall27Dep=np.sqrt(pcStall27Dep[0,0])



poStall3Dep, pcStall3Dep = curve_fit(N,Stall_J_3_DEP[0:l,0],Stall_J_3_DEP[0:l,1]*nmax)
tau_Stall3Dep=poStall3Dep[0]
err_Stall3Dep=np.sqrt(pcStall3Dep[0,0])

poStall32Dep, pcStall32Dep = curve_fit(N,Stall_J_32_DEP[0:l,0],Stall_J_32_DEP[0:l,1]*nmax)
tau_Stall32Dep=poStall32Dep[0]
err_Stall32Dep=np.sqrt(pcStall32Dep[0,0])

poStall35Dep, pcStall35Dep = curve_fit(N,Stall_J_35_DEP[0:l,0],Stall_J_35_DEP[0:l,1]*nmax)
tau_Stall35Dep=poStall35Dep[0]
err_Stall35Dep=np.sqrt(pcStall35Dep[0,0])

poStall37Dep, pcStall37Dep = curve_fit(N,Stall_J_37_DEP[0:l,0],Stall_J_37_DEP[0:l,1]*nmax)
tau_Stall37Dep=poStall37Dep[0]
err_Stall37Dep=np.sqrt(pcStall37Dep[0,0])



poStall4Dep, pcStall4Dep = curve_fit(N,Stall_J_4_DEP[0:l,0],Stall_J_4_DEP[0:l,1]*nmax)
tau_Stall4Dep=poStall4Dep[0]
err_Stall4Dep=np.sqrt(pcStall4Dep[0,0])

poStall42Dep, pcStall42Dep = curve_fit(N,Stall_J_42_DEP[0:l,0],Stall_J_42_DEP[0:l,1]*nmax)
tau_Stall42Dep=poStall42Dep[0]
err_Stall42Dep=np.sqrt(pcStall42Dep[0,0])

poStall45Dep, pcStall45Dep = curve_fit(N,Stall_J_45_DEP[0:l,0],Stall_J_45_DEP[0:l,1]*nmax)
tau_Stall45Dep=poStall45Dep[0]
err_Stall45Dep=np.sqrt(pcStall45Dep[0,0])

poStall47Dep, pcStall47Dep = curve_fit(N,Stall_J_47_DEP[0:l,0],Stall_J_47_DEP[0:l,1]*nmax)
tau_Stall47Dep=poStall47Dep[0]
err_Stall47Dep=np.sqrt(pcStall47Dep[0,0])

poStall5Dep, pcStall5Dep = curve_fit(N,Stall_J_5_DEP[0:l,0],Stall_J_5_DEP[0:l,1]*nmax)
tau_Stall5Dep=poStall5Dep[0]
err_Stall5Dep=np.sqrt(pcStall5Dep[0,0])

#%%

plt.title('Stall With_depletion ($N_{ss}=8$)',fontsize=17)
plt.xlabel('Simulation time',fontsize=14)
plt.ylabel('Stator number',fontsize=14)
plt.rcParams["figure.figsize"] = [8.0,6.0]

plt.plot(Stall_J_0_DEP[0:l,0],Stall_J_0_DEP[0:l,1]*nmax,label=r'J=0 $k_BT$, $\tau$={:1.2f}'.format(tau_Stall0Dep))
plt.plot(Stall_J_1_DEP[0:l,0],Stall_J_1_DEP[0:l,1]*nmax,label=r'J=1 $k_BT$, $\tau$={:1.2f}'.format(tau_Stall1Dep))
plt.plot(Stall_J_2_DEP[0:l,0],Stall_J_2_DEP[0:l,1]*nmax,label=r'J=2 $k_BT$, $\tau$={:1.2f}'.format(tau_Stall2Dep))
plt.plot(Stall_J_3_DEP[0:l,0],Stall_J_3_DEP[0:l,1]*nmax,label=r'J=3 $k_BT$, $\tau$={:1.2f}'.format(tau_Stall3Dep))
plt.plot(Stall_J_4_DEP[0:l,0],Stall_J_4_DEP[0:l,1]*nmax,label=r'J=4 $k_BT$, $\tau$={:1.2f}'.format(tau_Stall4Dep))
plt.plot(Stall_J_5_DEP[0:l,0],Stall_J_5_DEP[0:l,1]*nmax,label=r'J=5 $k_BT$, $\tau$={:1.2f}'.format(tau_Stall5Dep))

#plt.hlines(8.05,0,200, linestyles='dashed',alpha=0.4)    

plt.legend(loc='right', prop={'size': 15})


#%%No depletion

J=np.arange(0,5.25,0.25)


Relax_t_Res = np.array([tau_Res0,tau_Res02,tau_Res05,tau_Res07,tau_Res1,tau_Res12,tau_Res15,tau_Res17,tau_Res2,tau_Res22,
                        tau_Res25,tau_Res27,tau_Res3,tau_Res32,tau_Res35,tau_Res37,tau_Res4,tau_Res42,tau_Res45,tau_Res47,tau_Res5])

avg_Relax_t_Res = sum(Relax_t_Res)/len(Relax_t_Res)

#std_Relax_t_Res = np.sqrt(sum((Relax_t_Res-avg_Relax_t_Res)**2)/len(Relax_t_Res))

std_Relax_t_Res=np.array([err_Res0,err_Res02,err_Res05,err_Res07,err_Res1,err_Res12,err_Res15,err_Res17,err_Res2,err_Res22,
                          err_Res25,err_Res27,err_Res3,err_Res32,err_Res35,err_Res37,err_Res4,err_Res42,err_Res45,err_Res47,err_Res5])

Relax_t_Stall = np.array([tau_Stall0,tau_Stall02,tau_Stall05,tau_Stall07,tau_Stall1,tau_Stall12,tau_Stall15,tau_Stall17,
                          tau_Stall2,tau_Stall22,tau_Stall25,tau_Stall27,tau_Stall3,tau_Stall32,tau_Stall35,tau_Stall37,
                          tau_Stall4,tau_Stall42,tau_Stall45,tau_Stall47,tau_Stall5])

avg_Relax_t_Stall = sum(Relax_t_Stall)/len(Relax_t_Stall)

#std_Relax_t_Stall = np.sqrt(sum((Relax_t_Stall - avg_Relax_t_Stall)**2)/len(Relax_t_Stall))

std_Relax_t_Stall = np.array([err_Stall0,err_Stall02,err_Stall05,err_Stall07,err_Stall1,err_Stall12,err_Stall15,err_Stall17,
                              err_Stall2,err_Stall22,err_Stall25,err_Stall27,err_Stall3,err_Stall32,err_Stall35,err_Stall37,
                              err_Stall4,err_Stall42,err_Stall45,err_Stall47,err_Stall5])

#%% Depletion

J=np.arange(0,5.25,0.25)


Relax_t_Res_DEP = np.array([tau_Res0Dep,tau_Res02Dep,tau_Res05Dep,tau_Res07Dep,tau_Res1Dep,tau_Res12Dep,tau_Res15Dep,
                        tau_Res17Dep,tau_Res2Dep,tau_Res22Dep,tau_Res25Dep,tau_Res27Dep,tau_Res3Dep,tau_Res32Dep,
                        tau_Res35Dep,tau_Res37Dep,tau_Res4Dep,tau_Res42Dep,tau_Res45Dep,tau_Res47Dep,tau_Res5Dep])

avg_Relax_t_Res_DEP = sum(Relax_t_Res_DEP)/len(Relax_t_Res_DEP)

#std_Relax_t_Res = np.sqrt(sum((Relax_t_Res-avg_Relax_t_Res)**2)/len(Relax_t_Res))

std_Relax_t_Res_DEP=np.array([err_Res0Dep,err_Res02Dep,err_Res05Dep,err_Res07Dep,err_Res1Dep,err_Res12Dep,err_Res15Dep,
                              err_Res17Dep,err_Res2Dep,err_Res22Dep,err_Res25Dep,err_Res27Dep,err_Res3Dep,err_Res32Dep,
                              err_Res35Dep,err_Res37Dep,err_Res4Dep,err_Res42Dep,err_Res45Dep,err_Res47Dep,err_Res5Dep])

Relax_t_Stall_DEP = np.array([tau_Stall0Dep,tau_Stall02Dep,tau_Stall05Dep,tau_Stall07Dep,tau_Stall1Dep,tau_Stall12Dep,
                              tau_Stall15Dep,tau_Stall17Dep,tau_Stall2Dep,tau_Stall22Dep,tau_Stall25Dep,tau_Stall27Dep,
                              tau_Stall3Dep,tau_Stall32Dep,tau_Stall35Dep,tau_Stall37Dep,tau_Stall4Dep,tau_Stall42Dep,
                              tau_Stall45Dep,tau_Stall47Dep,tau_Stall5Dep])

avg_Relax_t_Stall_DEP = sum(Relax_t_Stall_DEP)/len(Relax_t_Stall_DEP)

#std_Relax_t_Stall = np.sqrt(sum((Relax_t_Stall - avg_Relax_t_Stall)**2)/len(Relax_t_Stall))

std_Relax_t_Stall_DEP = np.array([err_Stall0Dep,err_Stall02Dep,err_Stall05Dep,err_Stall07Dep,err_Stall1Dep,err_Stall12Dep,
                                  err_Stall15Dep,err_Stall17Dep,err_Stall2Dep,err_Stall22Dep,err_Stall25Dep,
                                  err_Stall27Dep,err_Stall3Dep,err_Stall32Dep,err_Stall35Dep,err_Stall37Dep,
                                  err_Stall4Dep,err_Stall42Dep,err_Stall45Dep,err_Stall47Dep,err_Stall5Dep])


#%%

#np.savetxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/RelaxTimes/ResTime_8_n09_dist.dat',(np.c_[J,Relax_t_Res,std_Relax_t_Res]))
#np.savetxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/RelaxTimes/StallTime_8_n09_dist.dat',(np.c_[J,Relax_t_Stall,std_Relax_t_Stall]))


np.savetxt('/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/RelaxTimes/ResTime_8.dat',(np.c_[J,Relax_t_Res_DEP,std_Relax_t_Res_DEP]))
np.savetxt('/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/RelaxTimes/StallTime_8_n09.dat',(np.c_[J,Relax_t_Stall_DEP,std_Relax_t_Stall_DEP]))

#%%

#Dtau05Res = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Relaxation times/ResTime_J=-mu.dat')
#Dtau05Stall = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Relaxation times/StallTime_J=-mu.dat')


Dtau03Res_DEP = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/RelaxTimes/ResTime_4.dat')
Dtau03Stall_DEP = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/RelaxTimes/StallTime_4_n06.dat')

#Dtau07Res = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/RelaxTimes/ResTime_8.dat')
#Dtau07Stall = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/RelaxTimes/StallTime_8_n09.dat')

#%%

Dtau04Res_dist = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/RelaxTimes/ResTime_4_n06_dist.dat')
Dtau04Stall_dist = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/RelaxTimes/StallTime_4_n06_dist.dat')

Dtau04Res = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/RelaxTimes/ResTime_4_n06.dat')
Dtau04Stall = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/RelaxTimes/StallTime_4_n06.dat')


Dtau08Res_dist = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/RelaxTimes/ResTime_8_n09_dist.dat')
Dtau08Stall_dist = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/RelaxTimes/StallTime_8_n09_dist.dat')

Dtau08Res = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/RelaxTimes/ResTime_8_n09.dat')
Dtau08Stall = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/RelaxTimes/StallTime_8_n09.dat')


Dtau10Res_dist = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/RelaxTimes/ResTime_10_n010_dist.dat')
Dtau10Stall_dist = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/RelaxTimes/StallTime_10_n010_dist.dat')

Dtau10Res = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/RelaxTimes/ResTime_10_n010.dat')
Dtau10Stall = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Glauber/RelaxTimes/StallTime_10_n010.dat')

#%%

Dtau04Res_DEP = np.loadtxt('/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/RelaxTimes/ResTime_4.dat')
Dtau04Stall_dist_DEP = np.loadtxt('/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/RelaxTimes/StallTime_4_n06.dat')

Dtau08Res_DEP = np.loadtxt('/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/RelaxTimes/ResTime_8.dat')
Dtau08Stall_dist_DEP = np.loadtxt('/home/mariajose/Escritorio/Simulations/With_depletion/Glauber/RelaxTimes/StallTime_8_n09.dat')


#%%

ax = plt.subplot(111)
#ax.set_xscale("log", nonposx='clip')
ax.set_yscale("log", nonposy='clip')

ax.set_title('$<N_{\infty}>=4$, $N_{max}=13$, $n_{tot}=13$',fontsize=17)
ax.set_xlabel(r'J  ($k_BT$)',fontsize=14)
ax.set_ylabel('$t_c$ (a.u.)',fontsize=14)
#plt.rcParams["figure.figsize"] = [8.0,6.0]

def r(p,a,b):
    return a*np.exp(b*p)

cJ=np.linspace(2, 20,1000)

poRes4, pcRes4 = curve_fit(r,J[8:],Dtau04Res_dist[8:,1])
poStall4, pcStall4 = curve_fit(r,J[8:],Dtau04Stall_dist[8:,1])
modelRes4=r(cJ,*poRes4)
modelStall4=r(cJ,*poStall4)


ax.plot(J,Dtau04Res_dist[:,1],marker='o',color='blue',linestyle='')
ax.plot(J,Dtau04Stall_dist[:,1],marker='^',color='orange',linestyle='')
ax.plot(cJ,modelRes4,color='blue',alpha=0.4)
ax.plot(cJ,modelStall4,color='orange',alpha=0.4)


black_line1 = mlines.Line2D([],[],linestyle='', color='blue', marker='o', label='Resurrection $n_0=0$')
black_line2 = mlines.Line2D([],[],linestyle='', color='orange', marker='^', label='Release $<n_0>=6$')
black_line3 = mlines.Line2D([],[],linestyle='-', color='blue', label=r'Slope={:1.2f}'.format(poRes4[0]))
black_line4 = mlines.Line2D([],[],linestyle='-', color='orange', label=r'Slope={:1.2f}'.format(poStall4[0]))
blue_patch = mpatches.Patch(color='blue', label='Resurrection')
orange_patch = mpatches.Patch(color='orange', label='Stall')
circles = mpatches.Patch( label='Stall')
plt.legend(loc='upper left', prop={'size': 15},handles=[black_line1,black_line3,black_line2,black_line4])

#%%

poRes4, pcRes4 = curve_fit(r,J[8:],Dtau04Res_dist[8:,1])
poStall4, pcStall4 = curve_fit(r,J[8:],Dtau04Stall_dist[8:,1])
modelRes4=r(cJ,*poRes4)
modelStall4=r(cJ,*poStall4)


poRes8, pcRes8 = curve_fit(r,J[8:],Dtau08Res_dist[8:,1])
poStall8, pcStall8 = curve_fit(r,J[8:],Dtau08Stall_dist[8:,1])
modelRes8=r(cJ,*poRes8)
modelStall8=r(cJ,*poStall8)

Dtau4_model = (modelStall4-modelRes4)/(modelStall4+modelRes4)
Dtau8_model = (modelStall8-modelRes8)/(modelStall8+modelRes8)

np.max(np.where(Dtau8_model<0.63))

cJ[982]


plt.plot(cJ,Dtau4_model)
plt.plot(cJ,Dtau8_model)


#%%

#Dtau05Res = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Relaxation times/ResTime_J=-mu.dat')
#Dtau05Stall = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Relaxation times/StallTime_J=-mu.dat')


#Dtau03Res = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Relaxation times/ResTime_4.dat')
#Dtau03Stall = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Relaxation times/StallTime_4.dat')

#Dtau07Res = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Relaxation times/ResTime_8.dat')
#Dtau07Stall = np.loadtxt('/home/mariajose/Escritorio/Simulations/Without _depletion/Relaxation times/StallTime_8.dat')

Dtau04 = (Dtau04Stall_dist[:,1]-Dtau04Res_dist[:,1])/((Dtau04Res_dist[:,1]+Dtau04Stall_dist[:,1]))
Dtau08 = (Dtau08Stall_dist[:,1]-Dtau08Res_dist[:,1])/((Dtau08Stall_dist[:,1]+Dtau08Res_dist[:,1]))


plt.title('$N_{max}=13$, no depletion',fontsize=17)
plt.xlabel(r'J ($k_B T$)', fontsize=14)
plt.ylabel(r'$\frac{t_c^{stall}-t_c^{res}}{<t_c>}$',fontsize=16)

plt.scatter(J, Dtau04, label=r'$<N_\infty>=4$',marker='^')
plt.scatter(J,Dtau08, label=r'$<N_\infty>=8$',marker='o')
#plt.scatter(J,Dtau07,label='$\Phi=2/3$',marker='s' )


plt.legend(loc='upper right', prop={'size': 15})
