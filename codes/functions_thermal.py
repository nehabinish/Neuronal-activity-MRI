#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 14:31:12 2021

@author: nehabinish
"""


import numpy as np
import init
from scipy import linalg as LA2


def control(time):
    
    """
    
      Generating Control fields to generate spin-echo sequences 
      and revive the relaxation of the magnetization vector.

        |  :param
        |-----------
        |      time      : array - length of evolution of time
        |                  time for temporal evolution of the spin ensemble
        |                  in the presence of thermal noises. 
        |
        |  :return
        |-----------
        |      controlx : float
        |                 control field of the x component of 
        |                 the pertubative magnetic field.
        |
        |      controly : float
        |                 control field of the y component of 
        |                 the pertubative magnetic field.
        

    """
    
    wx = 10*10**(5)                         # amplitude of control pulse on x axis
    wy = 100*10**(5)                         # amplitude of control pulse on y axis
    

    
    Np = 5
    
    controlx = np.zeros_like(time)
    controly = np.zeros_like(time)
    
    Fs = len(time)/Np

    controlx[::int(Fs)] = wx
    controly[::int(Fs)] = wy
            
        
    return (controlx, controly)


def integratemat(Mxyz, time, omegax, omegay):
    
    """
    
      Models the dynamics of the thermal evolution using the Bloch matrix

        |  :param
        |-----------
        |      Mxyz      : array - shape (3)
        |                  The initial components of the Magnetization vector
        |                  [Mx, My, Mz]
        |
        |      time      : array - length of evolution of time
        |                  time for temporal evolution of the spin ensemble
        |                  in the presence of thermal noises.
        |
        |      omegax    : float
        |                 control field of the x component of 
        |                 the pertubative magnetic field.
        |
        |      omegay    : float
        |                 control field of the y component of 
        |                 the pertubative magnetic field.
        |
        |  :return
        |-----------
        |      M        : array, shape( length(time), length(Mxyz) )
        |                 Temporal Evolution of the magnetization vector
        |                 in the presence of thermal noises.
        |
     
    """
    
    norm0   = 1 / np.sqrt( 1 + Mxyz[0]**2 + Mxyz[1]**2 + Mxyz[2]**2)
    
    initial = norm0 * np.array([1, Mxyz[0], Mxyz[1], Mxyz[2]])
    
    M       = np.zeros( ( len(time), len(initial) ) )
    
    N       = len(time)
    
    M[0]    = initial

    
    for i in range(N-1):
        
        
        A = np.array( [[ 1          , 0                              , 0                            , 0                          ], 
                       [ 0          , (-1/init.T2)                   , - init.omega_th              ,   (init.alpha * omegay[i]) ], 
                       [ 0          , init.omega_th                  , (-1/init.T2)                 , - (init.alpha * omegax[i]) ], 
                       [ (1/init.T1), - (init.alpha * omegay[i])     , (init.alpha * omegax[i])     , - (1/init.T1)              ] ])
         
        norm = 1 / np.sqrt( M[i,0]**2 + M[i,1]**2 + M[i,2]**2 + M[i,3]**2) 
        
        M[i+1] = norm * np.dot( LA2.expm(A * init.Deltat),  M[i])
        
        
    return M   
        



def integratemat_noises(Mxyz, time, omegax, omegay, f_noise):
    
    """
    
      Models the dynamics of the thermal evolution using the Bloch matrix

        |  :param
        |-----------
        |      Mxyz      : array - shape (3)
        |                  The initial components of the Magnetization vector
        |                  [Mx, My, Mz]
        |
        |      time      : array - length of evolution of time
        |                  time for temporal evolution of the spin ensemble
        |                  in the presence of thermal noises.
        |
        |      omegax    : float
        |                 control field of the x component of 
        |                 the pertubative magnetic field.
        |
        |      omegay    : float
        |                 control field of the y component of 
        |                 the pertubative magnetic field.
        |
        |      f_noise   : array
        |	               chaotic electromagnetic noises 
        |                  generated from the dynamical systems.
        |
        |  :return
        |-----------
        |      M        : array, shape( length(time), length(Mxyz) )
        |                 Temporal Evolution of the magnetization vector
        |                 in the presence of thermal noises.
        |
     
    """
    
    
    
    norm0    = 1 / np.sqrt( 1 + Mxyz[0]**2 + Mxyz[1]**2 + Mxyz[2]**2)
    
    initial  = norm0 * np.array([1, Mxyz[0], Mxyz[1], Mxyz[2]])
    
    M        = np.zeros( ( len(time), len(initial) ) )
    
    time_len = len(time)
    
    M[0]     = initial
    
    
    for n in range( time_len - 1):
        
        for j in range (0, init.N):

##---------------------------------------------------------------------------------------------------
#
### Uncomment of the choice is hindmarsh-rose model.
#            
#            noisey = (init.alpha_n * omegay[n]) + ( init.epsilon * f_noise[j,n,1] )
#            noisex = (init.alpha_n * omegax[n]) + ( init.epsilon * f_noise[j,n,0] )
#            noisez = init.omega_th + ( init.epsilon * f_noise[j,n,2] )
#            
#            H       = np.array([ [ 1          ,  0            , 0           , 0             ], 
#                                 [ 0          ,  (-1/init.T2) , - noisez    ,   noisey      ], 
#                                 [ 0          ,   noisez      , (-1/init.T2), - noisex      ], 
#                                 [ (1/init.T1), - noisey      , noisex      , - (1/init.T1) ] ])
#                               
#             
 ##--------------------------------------------------------------------------------------------------
 
            noisey = (init.alpha_n * omegay[n]) + ( init.epsilon * f_noise[j,n,1] )
            noisex = (init.alpha_n * omegax[n]) + ( init.epsilon * f_noise[j,n,0] )
            
            H       = np.array([ [ 1          ,  0             , 0              , 0             ], 
                                 [ 0          ,  (-1/init.T2)  , - init.omega_th,   noisey      ], 
                                 [ 0          ,   init.omega_th, (-1/init.T2)   , - noisex      ], 
                                 [ (1/init.T1), - noisey       ,   noisex       , - (1/init.T1) ] ])
                
                
            norm   = 1 / np.sqrt( M[n,0]**2 + M[n,1]**2 + M[n,2]**2 + M[n,3]**2) 
                
                
            M[n+1] = norm * np.dot( LA2.expm( H * init.Deltat_n),  M[n])
            
        
    return M  
