#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 16:11:39 2022
Last modification on Thu Jun 9
@author: Nils-Ole Walliser
"""

import numpy as np
#import matplotlib.pyplot as plt

""" Functions """
def decimalToBinary(n):
    # converting decimal to
    # and removng the prefix 0b
    return bin(n).replace("0b","")

def intVectorToBinNumMatrix(nums):
    # https://www.w3resource.com/python-exercises/numpy/python-numpy-exercise-187.php
    # converting a vector of integer numbers
    # into a matrix whose row are the representations
    # of those integer numbers in the binary system
    # Ex.: [1,2,3] => [[0,1],[1,0],[1,1]]
    max = np.max(nums)
    dimRow = int(round(np.log(max)/np.log(2)))
    return ((nums.reshape(-1,1) & (2**np.arange(dimRow))) != 0).astype(int)

def occVector(microstateMatrix):
    # extracts the occupancy vector from the microstate Matrix
    # Ex.: [[0,0],[1,0],[0,1],[1,1]] => [0,1,1,2]
    size = len(microstateMatrix)
    return np.array([np.dot(microstateMatrix[i],microstateMatrix[i]) \
                     for i in range(size)])

def position(occupationVector, n):
    # find position of "n" inside vector occupationVector
    return np.array(np.where(occupationVector == n)[0])

def detDownStates(microstate):
    # identifies the microstates with occupationNum-1 from which the given 
    # microstate can be created by the adsorption of one stator
    # Ex. 1011 => [0011,1001,1010]
    # How it works: select among all microstates with occupationNum-1 those
    # for which their scalar product with the given microstate
    # equals occupationNum-1
    
    # occupation number of MS
    occNum = microstate.sum(axis=0)
    assert occNum > 0
    assert len(microstate) < Nmax + 1
    
    # nb of microstates with occNum-1
    scan = len(position(occupationVector,occNum-1)) 
    downStates = []
    
    for i in range(scan):
        msTrans = microstateMatrix[position(occupationVector,occNum-1)][i]
        prod = np.dot(msTrans, microstate)
        if prod == occNum-1:
            #print(msTrans)
            input = np.where((microstateMatrix == msTrans).all(axis=1))[0]
            #print(input)
            downStates.append(input)
    return np.reshape(np.array(downStates),-1)                


def detUpStates(microstate):
    # identifies the microstates with occupationNum+1 from which the given 
    # microstate can be created by the desorption of one stator
    # Ex. 1011 => [1001,1010,0011]
    # How it works: select among all microstates with occupationNum+1 those
    # for which their scalar product with the given microstate
    # equals occupationNum
    
    # occupation number of MS
    occNum = microstate.sum(axis=0)
    assert occNum < Nmax
    assert len(microstate) < Nmax + 1
    
    # nb of microstates with occNum-1
    scan = len(position(occupationVector,occNum+1)) 
    upStates = []
    
    for i in range(scan):
        msTrans = microstateMatrix[position(occupationVector,occNum+1)][i]
        prod = np.dot(msTrans, microstate)
        if prod == occNum:
            #print(msTrans)
            input = np.where((microstateMatrix == msTrans).all(axis=1))[0]
            #print(input)
            upStates.append(input)
    return np.reshape(np.array(upStates),-1)     
    


def stoMatFill(microstateMatrix, stoMatrix, microstateNum):
    # CORE ROUTINE: creates the row of the stochastic matrix corresponding
    # to the microstate microstateNum
    
    lenB = len(microstateMatrix[1]); lenSM = len(stoMatrix)
    
    assert isinstance(microstateNum,int)
    assert microstateNum > -1
    assert microstateNum < lenSM

    # first row corresponding to the single zero-occupation microstate
    if microstateNum == 0:
        # loss
        stoMatrix[0,0] = -lenB*a0
        # gain
        for i in range(len(position(occupationVector,1))):
            stoMatrix[0,position(occupationVector,1)[i]] = b0
            
    # last row corresponding to the single full-occupation microstate
    elif microstateNum == lenSM-1:
        # loss
        stoMatrix[lenSM-1,lenSM-1] = -lenB*b2
        # gain (remember PBC)
        for i in range(len(position(occupationVector,lenB-1))):
            stoMatrix[lenSM-1,position(occupationVector,lenB-1)[i]] = a2
    
    # fill the bulk of the matrix ... this is the diffult part
    else:
        pass
    
    return stoMatrix

""" Control parameters """
# number of binding sites (> 2 needed to avoid ambiguity with PBC)
Nmax = 4
#assert Nmax > 2, f"number greater than 2 expected, got: Nmax = {Nmax}"

# kinetic rates
# TODO: express them i.t.o. chem. pot. mu and interaction pot. J
a0 = 0.5; a1 = 1; a2 = 2*a1
b0 = 0.5; b1 = 1; b2 = 2*b1

""" Main """
enumMicrostates =  np.arange(0,2**Nmax)
microstateMatrix = intVectorToBinNumMatrix(enumMicrostates)
occupationVector = occVector(microstateMatrix)

# initialize stochastic matrix
stoMatrix = np.zeros([len(microstateMatrix),len(microstateMatrix)])

# add entries to stochastic matrix

# 0th row
stoMatFill(microstateMatrix, stoMatrix, 0)
stoMatFill(microstateMatrix, stoMatrix, len(enumMicrostates)-1)

""" Print outs / diagnostics """
print(enumMicrostates)
print("microstates:\n",microstateMatrix)
print("occupation:\n",occupationVector)
#print(position(occupationVector,3), len(position(occupationVector,3)))
#print(stoMatrix)


#%%

Nmicro = len(enumMicrostates)

for i in range(Nmicro):
    
    # State we are looking at
    state = microstateMatrix[i] 
    
    # Down states
    dstates = detDownStates(state)
    
    # Up states
    ustates = detUpStates(state)
    
    # Study of transitions from and to down states
    for k in range(len(dstates)):
        
        # Row of the microstate matrix of the transition state
        j = dstates[k]
        
        # Transition state we are looking at
        transState = microstateMatrix[j]
        
        # Which site has gain or loose a stator 
        # final state - old state
        UpTrans = state - transState
        
        DownTrans = transState - state
        
        for l in range(Nmax):
            
            # There is a jump from the transition state to the state: 0 -> 1
            # We study then the neighbours of the transition state
            if(UpTrans[l] == 1):
                
                if(l==0):
                    
                    # No neigbours
                    if(state[j,Nmax-1] == state[j,l+1] == 0):
                        
                        a0=1
                    
                    # One neighbour
                    elif(state[j,Nmax-1] != state[j,l+1]):
                        
                        a1=1
                    
                    # Two neighbours
                    elif(state[j,Nmax-1] == state[j,l+1] == 1):
                    
                        a2=1
                        
                elif(l==Nmax):
                    
                    # No neigbours
                    if(state[j,0] == state[j,l-1] == 0):
                        
                        a0=1
                    
                    # One neighbour
                    elif(state[j,0] != state[j,l-1]):
                        
                        a1=1
                    
                    # Two neighbours
                    elif(state[j,0] == state[j,l-1] == 1):
                    
                        a2=1
                    
                else:
                     
                     # No neigbours
                    if(state[j,l+1] == state[j,l-1] == 0):
                        
                        a0=1
                    
                    # One neighbour
                    elif(state[j,l+1] != state[j,l-1]):
                        
                        a1=1
                    
                    # Two neighbours
                    elif(state[j,l+1] == state[j,l-1] == 1):
                    
                        a2=1
                             
            # There is a jump from the state to the transition state: 1 -> 0
            # We study then the neighbours of the state            
            if(DownTrans[l]== -1):
                
                if(l==0):
                    
                    # No neigbours
                    if(state[i,Nmax-1] == state[i,l+1] == 0):
                        
                        b0=-1
                    
                    # One neighbour
                    elif(state[i,Nmax-1] != state[i,l+1]):
                        
                        b1=-1
                    
                    # Two neighbours
                    elif(state[i,Nmax-1] == state[i,l+1] == 1):
                    
                        b2=-1
                        
                elif(l==Nmax):
                    
                    # No neigbours
                    if(state[i,0] == state[i,l-1] == 0):
                        
                        b0=-1
                    
                    # One neighbour
                    elif(state[i,0] != state[i,l-1]):
                        
                        b1=-1
                    
                    # Two neighbours
                    elif(state[i,0] == state[i,l-1] == 1):
                    
                        b2=-1
                    
                else:
                     
                     # No neigbours
                    if(state[i,l+1] == state[i,l-1] == 0):
                        
                        b0=-1
                    
                    # One neighbour
                    elif(state[i,l+1] != state[i,l-1]):
                        
                        b1=-1
                    
                    # Two neighbours
                    elif(state[i,l+1] == state[i,l-1] == 1):
                    
                        b2=-1
                        
    for k in range(len(ustates)):
        
        # Row of the microstate matrix of the transition state
        j = ustates[k]
        
        # Transition state we are looking at
        transState = microstateMatrix[j]
        
        # Which site has gain or loose a stator 
        # final state - old state
        UpTrans = transState - state
        
        DownTrans = state - transState
        
        for l in range(Nmax):
            
            # There is a jump from the state to the transition state: 0 -> 1
            # We study then the neighbours of the state
            if(UpTrans[l] == 1):
                
                if(l==0):
                    
                    # No neigbours
                    if(state[i,Nmax-1] == state[i,l+1] == 0):
                        
                        a0=-1
                    
                    # One neighbour
                    elif(state[i,Nmax-1] != state[i,l+1]):
                        
                        a1=-1
                    
                    # Two neighbours
                    elif(state[i,Nmax-1] == state[i,l+1] == 1):
                    
                        a2=-1
                        
                elif(l==Nmax):
                    
                    # No neigbours
                    if(state[i,0] == state[i,l-1] == 0):
                        
                        a0=-1
                    
                    # One neighbour
                    elif(state[i,0] != state[i,l-1]):
                        
                        a1=-1
                    
                    # Two neighbours
                    elif(state[i,0] == state[i,l-1] == 1):
                    
                        a2=-1
                    
                else:
                     
                     # No neigbours
                    if(state[i,l+1] == state[i,l-1] == 0):
                        
                        a0=-1
                    
                    # One neighbour
                    elif(state[i,l+1] != state[i,l-1]):
                        
                        a1=-1
                    
                    # Two neighbours
                    elif(state[i,l+1] == state[i,l-1] == 1):
                    
                        a2=-1
                             
            # There is a jump from the transition state to the state: 1 -> 0
            # We study then the neighbours of the transition state            
            if(DownTrans[l]== -1):
                
                if(l==0):
                    
                    # No neigbours
                    if(state[j,Nmax-1] == state[j,l+1] == 0):
                        
                        b0=1
                    
                    # One neighbour
                    elif(state[j,Nmax-1] != state[j,l+1]):
                        
                        b1=1
                    
                    # Two neighbours
                    elif(state[j,Nmax-1] == state[j,l+1] == 1):
                    
                        b2=1
                        
                elif(l==Nmax):
                    
                    # No neigbours
                    if(state[j,0] == state[j,l-1] == 0):
                        
                        b0=1
                    
                    # One neighbour
                    elif(state[j,0] != state[j,l-1]):
                        
                        b1=1
                    
                    # Two neighbours
                    elif(state[j,0] == state[j,l-1] == 1):
                    
                        b2=1
                    
                else:
                     
                     # No neigbours
                    if(state[j,l+1] == state[j,l-1] == 0):
                        
                        b0=1
                    
                    # One neighbour
                    elif(state[j,l+1] != state[j,l-1]):
                        
                        b1=1
                    
                    # Two neighbours
                    elif(state[j,l+1] == state[j,l-1] == 1):
                    
                        b2=1
                
                
                        
                    
                
        
        