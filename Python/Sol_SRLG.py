#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 11 12:51:51 2022

@author: luke
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
plt.rcParams["figure.figsize"] = [10.0,8.0]

#%% Parameters and rates

J = 0.0
mu = 0.0

#a = rate from 0 -> 1
#b = rate from 1 -> 0

#0 = no neighbours
#1 = one neighbour
#2 = two neigbours

a0 = 1./2.*(1 - np.tanh(-mu/2.))
b0 =  1./2.*(1 - np.tanh(mu/2.))

a1 = 1./2.*(1 - np.tanh(-(J + mu)/2.))
b1 = 1./2.*(1 - np.tanh((J + mu)/2.))

a2 = 1./2.*(1 - np.tanh(-J-mu/2.))
b2 = 1./2.*(1 - np.tanh(J+mu/2.))


#%% Create a matrix

Nmax = 5 #Number of sites

TNmicro = 2**Nmax #Total number of micro-states

macro_s = np.arange(0, Nmax+1,1) #Possible macro-states

M = np.zeros((TNmicro,TNmicro)) #Initialization matrix

Nmacro = np.zeros(TNmicro) #Array of macro-states associated to micro-states

state = np.zeros((TNmicro,Nmax))

for i in range(TNmicro):
    
    state[i] = [int(d) for d in bin((1<<Nmax)+i)[-Nmax:]]


for j in range(TNmicro):
    
    Nmacro[j] = sum(state[j])


Nmicro = np.zeros(len(macro_s)) #Array of micro-states associate to macro-states

for i in range(len(macro_s)):
    
    Nmicro[i] = np.count_nonzero(Nmacro == float(i) )
    

print(Nmacro)
print(Nmicro)
print(state)


#%%

a = np.zeros((TNmicro,TNmicro))
b = np.zeros((TNmicro,TNmicro))
c = np.zeros((TNmicro,TNmicro))
d = np.zeros((TNmicro,TNmicro))
e = np.zeros((TNmicro,TNmicro))
f = np.zeros((TNmicro,TNmicro))



for i in range(TNmicro):
    
    for j in range(TNmicro):
        
        for k in range(Nmax):
                
            if(i==j):
            
                if((k-1)<0):
                
                    if(state[i,k] == state[i,Nmax-1] == state[i,k+1] == 0):
                    
                        a[i,j] = a[i,j] - 1
                    
                    elif(state[i,k] == state[i,Nmax-1] == state[i,k+1] == 1):
                       # if(Nmax==2):
                       #     d[i,j] = d[i,j] - 1
                       # else:
                            f[i,j] = f[i,j] - 1
                            
                    elif((state[i,k]==0)and(state[i,k+1]==state[i,Nmax-1]==1)):
                        #if(Nmax==2):
                        #    c[i,j] = c[i,j] - 1
                        #else:
                            e[i,j] = e[i,j] - 1
                            
                    elif((state[i,k]==0)and((state[i,k+1]==1!=state[i,Nmax-1])or(state[i,Nmax-1]==1!=state[i,k+1]))):

                        c[i,j] = c[i,j] - 1
                            
                    elif((state[i,k]==1)and(state[i,k+1]==state[i,Nmax-1]==0)):
             
                        b[i,j] = b[i,j] - 1
                            
                    elif((state[i,k]==1)and((state[i,k+1]==1!=state[i,Nmax-1])or(state[i,Nmax-1]==1!=state[i,k+1]))):

                        d[i,j] = d[i,j] - 1
                        
                elif((k+1)==Nmax):
                            
                    if(state[i,k] == state[i,0] == state[i,k-1] == 0):
                    
                        a[i,j] = a[i,j] - 1
                    
                    elif(state[i,k] == state[i,0] == state[i,k-1] == 1):
                      #  if(Nmax==2):
                      #      d[i,j] = d[i,j] - 1
                       # else:
                            f[i,j] = f[i,j] - 1
                            
                    elif((state[i,k]==0)and(state[i,k-1]==state[i,0]==1)):
                       # if(Nmax==2):
                       #     c[i,j] = c[i,j] - 1
                       # else:
                            e[i,j] = e[i,j] - 1
                            
                    elif((state[i,k]==0)and((state[i,0]==1!=state[i,k-1])or(state[i,k-1]==1!=state[i,0]))):
                     
                        c[i,j] = c[i,j] - 1
                            
                    elif((state[i,k]==1)and(state[i,k-1]==state[i,0]==0)):
             
                            
                        b[i,j] = b[i,j] - 1
                            
                    elif((state[i,k]==1)and((state[i,0]==1!=state[i,k-1])or(state[i,k-1]==1!=state[i,0]))):

                        d[i,j] = d[i,j] - 1
            
                else:

                    if(state[i,k] == state[i,k+1] == state[i,k-1] == 0):
                    
                        a[i,j] = a[i,j] - 1
                    
                    elif(state[i,k] == state[i,k+1] == state[i,k-1] == 1):

                        f[i,j] = f[i,j] - 1
                        
                    elif((state[i,k]==0)and(state[i,k+1]==state[i,k-1]==1)):
                       # if(Nmax==2):
                       #     c[i,j] = c[i,j] - 1
                       # else:
                            e[i,j] = e[i,j] - 1
                            
                    elif((state[i,k]==0)and((state[i,k+1]==1!=state[i,k-1])or(state[i,k-1]==1!=state[i,k+1]))):

                        c[i,j] = c[i,j] - 1
                            
                    elif((state[i,k]==1)and(state[i,k+1]==state[i,k-1]==0)):
             
                        b[i,j] = b[i,j] - 1
                            
                    elif((state[i,k]==1)and(state[i,k+1]==1!=state[i,k-1])or(state[i,k-1]==1!=state[i,k+1])):

                        d[i,j] = d[i,j] - 1
                        
            elif(i!=j):
                
                if((Nmacro[i]==0)and(Nmacro[j]==1)):
                    
                    b[i,j] = 1
                    
                elif((Nmacro[j]==0)and(Nmacro[i]==1)):
                    
                    a[i,j] = 1
                    
                else:
                    
                    if((state[i,k]==state[j,k]==1)and(Nmacro[i]==Nmacro[j]-1)):
                        
                        for l in range(Nmax):
                            
                            if((state[i,l]==0)and(state[j,l]==1)):
                            
                                if((l-1)<0):
                                
                                
                                    if(state[j,l]==state[j,Nmax-1]==state[j,l+1]==1):
                            
                                        f[i,j] = 1
                            
                                    elif((state[j,l]==1)and((state[j,Nmax-1]==1!=state[j,l+1])or(state[j,l+1]==1!=state[j,Nmax-1]))):
                            
                                        d[i,j] = 1
                                        
                                    elif((state[j,l]==1)and(state[j,Nmax-1]==state[j,l+1]==0)):
                            
                                        b[i,j] = 1
                                        
                                if((l+1)==Nmax):
                                    
                                    if(state[j,l]==state[j,0]==state[j,l-1]==1):
                             
                                        f[i,j] = 1
                             
                                    elif((state[j,l]==1)and((state[j,l-1]==1!=state[j,0])or(state[j,0]==1!=state[j,l-1]))):
                             
                                        d[i,j] = 1
                                    
                                    elif((state[j,l]==1)and(state[j,0]==state[j,l-1]==0)):
                             
                                        b[i,j] = 1
                                        
                                else:
                                    
                                    if(state[j,l]==state[j,l+1]==state[j,l-1]==1):
                             
                                        f[i,j] = 1
                             
                                    elif((state[j,l]==1)and((state[j,l-1]==1!=state[j,l+1])or(state[j,l+1]==1!=state[j,l-1]))):
                             
                                        d[i,j] = 1
                     
                                    elif((state[j,l]==1)and(state[j,l-1]==state[j,l-1]==0)):
                             
                                        b[i,j] = 1
                                        
                    if((state[i,k]==state[j,k]==1)and(Nmacro[i]==Nmacro[j]+1)):
                        
                        for l in range(Nmax):
                            
                            if((state[i,l]==1)and(state[j,l]==0)):
                                
                                if((l-1)<0):
                                
                                    if(state[j,l]==state[j,Nmax-1]==state[j,l+1]==0):
                            
                                        a[i,j] = 1
                            
                                    elif((state[j,l]==0)and(((state[j,Nmax-1]==1!=state[j,l+1])or(state[j,l+1]==1!=state[j,Nmax-1])))):
                            
                                        c[i,j] = 1
                                        
                                    elif((state[j,l]==0)and(state[j,Nmax-1]==state[j,l+1]==1)):
                            
                                        e[i,j] = 1
                                        
                                if((l+1)==Nmax):
                                    
                                    if(state[j,l]==state[j,0]==state[j,l-1]==0):
                             
                                        a[i,j] = 1
                             
                                    elif((state[j,l]==0)and((state[j,0]==1!=state[j,l-1])or(state[j,l-1]==1!=state[j,0]))):
                             
                                        c[i,j] = 1
                                    
                                    elif((state[j,l]==0)and(state[j,0]==state[j,l-1]==1)):
                             
                                        e[i,j] = 1
                                        
                                else:
                                    
                                    if(state[j,l]==state[j,l+1]==state[j,l-1]==0):
                             
                                        a[i,j] = 1
                             
                                    elif((state[j,l]==0)and((state[j,l+1]==1!=state[j,l-1])or(state[j,l-1]==1!=state[j,l+1]))):
                             
                                        c[i,j] = 1
                     
                                    elif((state[j,l]==0)and(state[j,l-1]==state[j,l-1]==1)):
                             
                                        e[i,j] = 1
                        
                        
                                    
                                
                                ##### continua con el caso de los alphas. No olvides la condicion de a=b=c...=0
                
                

        M[i,j] = a[i,j]*a0 + b[i,j]*b0 + c[i,j]*a1 + d[i,j]*b1 + e[i,j]*a2 + f[i,j]*b2
        

print(M)
                
#%%
                
t = sum(M[:,15])
    
print(t)

#print(a0,b0,a1,b1,a2,b2)


#%%
                
                if((state[i,k]==0)and(state[j,k]==1)and(Nmacro[i]==Nmacro[j]-1)and(a[i,j]==b[i,j]==c[i,j]==d[i,j]==e[i,j]==f[i,j]==0)):
                    
                    if((k-1)<0):
                        
                     
                  #      if(state[j,k]==state[j,Nmax-1]==state[j,k+1]==0):
                            
                  #          a[i,j] = 1
                            
                  #      elif((state[j,k]==0)and((state[j,Nmax-1]==1)or(state[j,k+1]==1))):
                            
                  #          c[i,j] = 1 
                    
                  #      elif((state[j,k]==0)and(state[j,Nmax-1]==state[j,k+1]==1)):
                            
                  #          e[i,j] = 1
                      
                        if(state[j,k]==state[j,Nmax-1]==state[j,k+1]==1):
                            
                            f[i,j] = 1
                            
                        elif((state[j,k]==1)and((state[j,Nmax-1]==1)or(state[j,k+1]==1))):
                            
                            d[i,j] = 1
                    
                        elif((state[j,k]==1)and(state[j,Nmax-1]==state[j,k+1]==0)):
                            
                            b[i,j] = 1
                            
                    if((k+1)==Nmax):
                         
                    #    if(state[j,k]==state[j,0]==state[j,k-1]==0):
                             
                    #        a[i,j] = 1
                             
                    #    elif((state[j,k]==0)and((state[j,0]==1)or(state[j,k-1]==1))):
                             
                    #        c[i,j] =  1
                     
                    #    elif((state[j,k]==0)and(state[j,0]==state[j,k-1]==1)):
                             
                    #        e[i,j] = 1
                           
                        if(state[j,k]==state[j,0]==state[j,k-1]==1):
                             
                            f[i,j] = 1
                             
                        elif((state[j,k]==1)and((state[j,0]==1)or(state[j,k-1]==1))):
                             
                            d[i,j] = 1
                     
                        elif((state[j,k]==1)and(state[j,0]==state[j,k-1]==0)):
                             
                            b[i,j] = 1
                            
                    else:
                        '''
                        if(state[j,k]==state[j,k+1]==state[j,k-1]==0):
                             
                            a[i,j] = 1
                             
                        elif((state[j,k]==0)and((state[j,k+1]==1)or(state[j,k-1]==1))):
                             
                            c[i,j] = 1
                     
                        elif((state[j,k]==0)and(state[j,k-1]==state[j,k-1]==1)):
                             
                            e[i,j] =  1
                         '''   
                        if(state[j,k]==state[j,k+1]==state[j,k-1]==1):
                             
                            f[i,j] = 1
                             
                        elif((state[j,k]==1)and((state[j,k+1]==1)or(state[j,k-1]==1))):
                             
                            d[i,j] = 1
                     
                        elif((state[j,k]==1)and(state[j,k-1]==state[j,k-1]==0)):
                             
                            b[i,j] = 1     
                            
                        
                            
                elif((state[i,k]==1)and(state[j,k]==0)and(Nmacro[i]==Nmacro[j]+1)and(a[i,j]==b[i,j]==c[i,j]==d[i,j]==e[i,j]==f[i,j]==0)):
    
                       
                       if((k-1)<0):
                           
                           if(state[j,k]==state[j,Nmax-1]==state[j,k+1]==0):
                               
                               a[i,j] = 1
                               
                           elif((state[j,k]==0)and((state[j,Nmax-1]==1)or(state[j,k+1]==1))):
                               
                               c[i,j] = 1 
                       
                           elif((state[j,k]==0)and(state[j,Nmax-1]==state[j,k+1]==1)):
                               
                               e[i,j] = 1
                          
                        #   elif(state[j,k]==state[j,Nmax-1]==state[j,k+1]==1):
                               
                        #       f[i,j] = 1
                               
                        #   elif((state[j,k]==1)and((state[j,Nmax-1]==1)or(state[j,k+1]==1))):
                               
                        #       d[i,j] = 1
                       
                        #   elif((state[j,k]==1)and(state[j,Nmax-1]==state[j,k+1]==0)):
                               
                        #       b[i,j] = 1
                              
                       if((k+1)==Nmax):
                           
                           if(state[j,k]==state[j,0]==state[j,k-1]==0):
                                
                               a[i,j] = 1
                                
                           elif((state[j,k]==0)and((state[j,0]==1)or(state[j,k-1]==1))):
                                
                               c[i,j] =  1
                        
                           elif((state[j,k]==0)and(state[j,0]==state[j,k-1]==1)):
                                
                               e[i,j] = 1
                               
                           '''      
                           elif(state[j,k]==state[j,0]==state[j,k-1]==1):
                                
                               f[i,j] = 1
                                
                           elif((state[j,k]==1)and((state[j,0]==1)or(state[j,k-1]==1))):
                                
                               d[i,j] = 1
                        
                           elif((state[j,k]==1)and(state[j,0]==state[j,k-1]==0)):
                                
                               b[i,j] = 1
                            '''  
                       else:
                           
                           if(state[j,k]==state[j,k+1]==state[j,k-1]==0):
                                
                               a[i,j] = 1
                                
                           elif((state[j,k]==0)and((state[j,k+1]==1)or(state[j,k-1]==1))):
                                
                               c[i,j] = 1
                        
                           elif((state[j,k]==0)and(state[j,k-1]==state[j,k-1]==1)):
                                
                               e[i,j] =  1
                        
                          # elif(state[j,k]==state[j,k+1]==state[j,k-1]==1):
                                
                          #     f[i,j] = 1
                                
                          # elif((state[j,k]==1)and((state[j,k+1]==1)or(state[j,k-1]==1))):
                                
                          #     d[i,j] = 1
                        
                          # elif((state[j,k]==1)and(state[j,k-1]==state[j,k-1]==0)):
                                
                          #     b[i,j] = 1     
                        
            
        M[i,j] = a[i,j]*a0 + b[i,j]*b0 + c[i,j]*a1 + d[i,j]*b1 + e[i,j]*a2 + f[i,j]*b2
        

print(M)

print(len(M))

#%%


    


print(a0,b0,a1,b1,a2,b2)                
#%%

    
def MM(Nmax):

    TNmicro = 2**Nmax #Total number of micro-states

    macro_s = np.arange(0, Nmax+1,1) #Possible macro-states

    MM = np.zeros((TNmicro,TNmicro)) #Initialization matrix

    Nmacro = np.zeros(TNmicro) #Array of macro-states associated to micro-states

    state = np.zeros((TNmicro,Nmax))

    for i in range(TNmicro):
    
        state[i] = [int(d) for d in bin((1<<8)+i)[-Nmax:]]


    for j in range(TNmicro):
    
        Nmacro[j] = sum(state[j])


    Nmicro = np.zeros(len(macro_s)) #Array of micro-states associate to macro-states

    for i in range(len(macro_s)):
    
        Nmicro[i] = np.count_nonzero(Nmacro == float(i) )
        
        
        
    a = np.zeros((TNmicro,TNmicro))
    b = np.zeros((TNmicro,TNmicro))
    c = np.zeros((TNmicro,TNmicro))
    d = np.zeros((TNmicro,TNmicro))
    e = np.zeros((TNmicro,TNmicro))
    f = np.zeros((TNmicro,TNmicro))



    for i in range(TNmicro):
        
        for j in range(TNmicro):
        
            for k in range(Nmax):
                
                if(i==j):
            
                    if((k-1)<0):
                
                        if(state[i,k] == state[i,Nmax-1] == state[i,k+1] == 0):
                    
                            a[i,j] = a[i,j] - 1
                    
                        elif(state[i,k] == state[i,Nmax-1] == state[i,k+1] == 1):
                            if(Nmax==2):
                                d[i,j] = d[i,j] - 1
                            else:
                                f[i,j] = f[i,j] - 1
                            
                        elif((state[i,k]==0)and(state[i,k+1]==state[i,Nmax-1]==1)):
                            if(Nmax==2):
                                c[i,j] = c[i,j] - 1
                            else:
                                e[i,j] = e[i,j] - 1
                            
                        elif((state[i,k]==0)and((state[i,k+1]== 1)or(state[i,Nmax-1]==1))):

                            c[i,j] = c[i,j] - 1
                            
                        elif((state[i,k]==1)and(state[i,k+1]==state[i,Nmax-1]==0)):
             
                            b[i,j] = b[i,j] - 1
                            
                        elif((state[i,k]==1)and((state[i,k+1] == 1)or(state[i,Nmax-1]== 1))):

                            d[i,j] = d[i,j] - 1
                        
                    elif((k+1)==Nmax):
                            
                        if(state[i,k] == state[i,0] == state[i,k-1] == 0):
                    
                            a[i,j] = a[i,j] - 1
                    
                        elif(state[i,k] == state[i,0] == state[i,k-1] == 1):
                            if(Nmax==2):
                                d[i,j] = d[i,j] - 1
                            else:
                                f[i,j] = f[i,j] - 1
                            
                        elif((state[i,k]==0)and(state[i,k-1]==state[i,0]==1)):
                            if(Nmax==2):
                                c[i,j] = c[i,j] - 1
                            else:
                                e[i,j] = e[i,j] - 1
                            
                        elif((state[i,k]==0)and((state[i,k-1]== 1)or(state[i,0]==1))):
                            
                            c[i,j] = c[i,j] - 1
                            
                        elif((state[i,k]==1)and(state[i,k-1]==state[i,0]==0)):
                            
                            b[i,j] = b[i,j] - 1
                            
                        elif((state[i,k]==1)and((state[i,k-1] == 1)or(state[i,0]== 1))):

                            d[i,j] = d[i,j] - 1
            
                    else:

                        if(state[i,k] == state[i,k+1] == state[i,k-1] == 0):
                    
                            a[i,j] = a[i,j] - 1
                    
                        elif(state[i,k] == state[i,k+1] == state[i,k-1] == 1):

                            f[i,j] = f[i,j] - 1
                        
                        elif((state[i,k]==0)and(state[i,k+1]==state[i,k-1]==1)):
                            if(Nmax==2):
                                c[i,j] = c[i,j] - 1
                            else:
                                e[i,j] = e[i,j] - 1
                            
                        elif((state[i,k]==0)and((state[i,k+1]==1)or(state[i,k-1]==1))):

                            c[i,j] = c[i,j] - 1
                            
                        elif((state[i,k]==1)and(state[i,k+1]==state[i,k-1]==0)):
             
                            b[i,j] = b[i,j] - 1
                            
                        elif((state[i,k]==1)and((state[i,k+1] == 1)or(state[i,k-1]== 1))):

                            d[i,j] = d[i,j] - 1
                        
                elif(i!=j):
                
                    if((state[i,k]!=state[j,k])and(Nmacro[i]==Nmacro[j]-1)and(a[i,j]==b[i,j]==c[i,j]==d[i,j]==e[i,j]==f[i,j]==0)):
                    
                        if((k-1)<0):
                        
                            if(state[j,k]==state[j,Nmax-1]==state[j,k+1]==0):
                            
                                a[i,j] = 1
                            
                            elif((state[j,k]==0)and((state[j,Nmax-1]==1)or(state[j,k+1]==1))):
                            
                                c[i,j] = 1 
                    
                            elif((state[j,k]==0)and(state[j,Nmax-1]==state[j,k+1]==1)):
                            
                                e[i,j] = 1
                            
                            elif(state[j,k]==state[j,Nmax-1]==state[j,k+1]==1):
                            
                                f[i,j] = 1
                                
                            elif((state[j,k]==1)and((state[j,Nmax-1]==1)or(state[j,k+1]==1))):
                            
                                d[i,j] = 1
                    
                            elif((state[j,k]==1)and(state[j,Nmax-1]==state[j,k+1]==0)):
                            
                                b[i,j] = 1
                            
                        if((k+1)==Nmax):
                         
                            if(state[j,k]==state[j,0]==state[j,k-1]==0):
                             
                                a[i,j] = 1
                             
                            elif((state[j,k]==0)and((state[j,0]==1)or(state[j,k-1]==1))):
                             
                                c[i,j] =  1
                     
                            elif((state[j,k]==0)and(state[j,0]==state[j,k-1]==1)):
                             
                                e[i,j] = 1
                            
                            elif(state[j,k]==state[j,0]==state[j,k-1]==1):
                             
                                f[i,j] = 1
                             
                            elif((state[j,k]==1)and((state[j,0]==1)or(state[j,k-1]==1))):
                             
                                d[i,j] = 1
                     
                            elif((state[j,k]==1)and(state[j,0]==state[j,k-1]==0)):
                             
                                b[i,j] = 1
                            
                        else:
                        
                            if(state[j,k]==state[j,k+1]==state[j,k-1]==0):
                             
                                a[i,j] = 1
                             
                            elif((state[j,k]==0)and((state[j,k+1]==1)or(state[j,k-1]==1))):
                             
                                c[i,j] = 1
                     
                            elif((state[j,k]==0)and(state[j,k-1]==state[j,k-1]==1)):
                             
                                e[i,j] =  1
                            
                            elif(state[j,k]==state[j,k+1]==state[j,k-1]==1):
                             
                                f[i,j] = 1
                             
                            elif((state[j,k]==1)and((state[j,k+1]==1)or(state[j,k-1]==1))):
                             
                                d[i,j] = 1
                     
                            elif((state[j,k]==1)and(state[j,k-1]==state[j,k-1]==0)):
                             
                                b[i,j] = 1     
                            
                        
                            
                    elif((state[i,k]!=state[j,k])and(Nmacro[i]==Nmacro[j]+1)and(a[i,j]==b[i,j]==c[i,j]==d[i,j]==e[i,j]==f[i,j]==0)):
    
                       
                        if((k-1)<0):
                           
                           if(state[j,k]==state[j,Nmax-1]==state[j,k+1]==0):
                               
                               a[i,j] = 1
                               
                           elif((state[j,k]==0)and((state[j,Nmax-1]==1)or(state[j,k+1]==1))):
                               
                               c[i,j] = 1 
                       
                           elif((state[j,k]==0)and(state[j,Nmax-1]==state[j,k+1]==1)):
                               
                               e[i,j] = 1
                           
                           elif(state[j,k]==state[j,Nmax-1]==state[j,k+1]==1):
                               
                               f[i,j] = 1
                               
                           elif((state[j,k]==1)and((state[j,Nmax-1]==1)or(state[j,k+1]==1))):
                               
                               d[i,j] = 1
                       
                           elif((state[j,k]==1)and(state[j,Nmax-1]==state[j,k+1]==0)):
                               
                               b[i,j] = 1
                               
                        if((k+1)==Nmax):
                           
                           if(state[j,k]==state[j,0]==state[j,k-1]==0):
                                
                               a[i,j] = 1
                                
                           elif((state[j,k]==0)and((state[j,0]==1)or(state[j,k-1]==1))):
                                
                               c[i,j] =  1
                        
                           elif((state[j,k]==0)and(state[j,0]==state[j,k-1]==1)):
                                
                               e[i,j] = 1
                            
                           elif(state[j,k]==state[j,0]==state[j,k-1]==1):
                                
                               f[i,j] = 1
                                
                           elif((state[j,k]==1)and((state[j,0]==1)or(state[j,k-1]==1))):
                                
                               d[i,j] = 1
                        
                           elif((state[j,k]==1)and(state[j,0]==state[j,k-1]==0)):
                                
                               b[i,j] = 1
                               
                        else:
                           
                           if(state[j,k]==state[j,k+1]==state[j,k-1]==0):
                                
                               a[i,j] = 1
                                
                           elif((state[j,k]==0)and((state[j,k+1]==1)or(state[j,k-1]==1))):
                                
                               c[i,j] = 1
                        
                           elif((state[j,k]==0)and(state[j,k-1]==state[j,k-1]==1)):
                                
                               e[i,j] =  1
                           
                           elif(state[j,k]==state[j,k+1]==state[j,k-1]==1):
                                
                               f[i,j] = 1
                                
                           elif((state[j,k]==1)and((state[j,k+1]==1)or(state[j,k-1]==1))):
                                
                               d[i,j] = 1
                        
                           elif((state[j,k]==1)and(state[j,k-1]==state[j,k-1]==0)):
                                
                               b[i,j] = 1     
                            
                    

            
        MM[i,j] = a[i,j]*a0 + b[i,j]*b0 + c[i,j]*a1 + d[i,j]*b1 + e[i,j]*a2 + f[i,j]*b2
    
   
    return MM 
    
    
MM(2)
