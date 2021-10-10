#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 11:33:51 2021

@author: nehabinish
"""

import numpy as np


hbar  = 1.                                         # in a.u.

# Pauli Matrices
Sx = np.array( [[0 ,   1], [1  ,  0]] )
Sy = np.array( [[0 , -1j], [1j ,  0]] )    
Sz = np.array( [[1 ,   0], [0  , -1]] )

N          = 10    

''' ========================================================================

     INITIAL PARAMETERS FOR FITZHUGH NAGUMO CLASSICAL DYNAMICAL SYSTEM

========================================================================== ''' 
                            

ntime_fitz = 10000

t0_fitz    = 0                                     # initial time 
tf_fitz    = 600                                   # final time 
dt_fitz    = (tf_fitz - t0_fitz) / ntime_fitz      # time step 
t_fitz     = np.arange(t0_fitz, tf_fitz, dt_fitz)  # time series for 
                                                   # chaotic evolution

''' ========================================================================

    INITIAL PARAMETERS FOR HINDMARSH-ROSE CLASSICAL DYNAMICAL SYSTEM

========================================================================== ''' 
                           

ntime_hind = 10000

t0_hind    = 0                                     # initial time 
tf_hind    = 200                                   # final time 
dt_hind    = (tf_hind - t0_hind) / ntime_hind      # time step 
t_hind     = np.arange(t0_hind, tf_hind, dt_hind)  # time series for 
                                                   # chaotic evolution


''' ========================================================================

         INITIAL PARAMETERS FOR QUANTUM SYSTEM - SCHRÖDINGER EQUATIONS

========================================================================== ''' 


ntime_sch   = ntime_fitz



t0_sch      = -30.                                   # initial time (a.u.)
tf_sch      =  200.                                  # final time (a.u.)

## ------------------------------------------------------------------------ ##
# t0_sch = 30. * 2.41888425 * 10**(-17)              # seconds
# tf_sch = 200 * 2.41888425 * 10**(-17)              # seconds
## ------------------------------------------------------------------------ ##

## ------------------------------------------------------------------------ ##
### Uncomment if the choice is hndmarsh-rose model
# t0_sch = -30.                                      # a.u.
# tf_sch = 5000.                                     # a.u.
## ------------------------------------------------------------------------ ##


Deltat_sch  = (tf_sch - t0_sch) / ntime_sch          # time step

t_sch       = np.linspace(t0_sch, tf_sch, ntime_sch) # time series for 
                                                     # integration of the 
                                                     # schrodinger equation 
                                                     
gamma = 1.414 * (10 ** (-3))                         # in a.u.
## ------------------------------------------------------------------------ ##
# gamma =  2.6752218744 * 10**(8)                    # rad⋅s−1⋅T−1.
## ------------------------------------------------------------------------ ##


epsilon = 10**10*gamma                                      # electromagnetic 
                                                     # coupling constant
                                                     
B0      = 10 ** (-5)                                 # a.u.
omega   = gamma * B0                                 # angular frequency 


''' ========================================================================

        INITIAL PARAMETERS FOR THE ADDITION OF ISOLATED THERMAL NOISES

========================================================================== ''' 


M_init    =  np.array([-1.5, 1.5, 0])               # initial position of 
                                                    # the magnetization vector

gamma_T   = 2.6752218744 * 10**(8)                  # gyromagnetic ratio
                                                    # rad⋅s−1⋅T−1.

B0_th     = 0.5                                     # a.u.
omega_th  = gamma_T * np.abs(B0_th)                 # larmor frequency

alpha     = 1                                       # coupling constant

T1        = 0.3                                     # T1 relaxation time (s)
T2        = 0.2                                     # T2 relaxation time (s)

time      = np.linspace(0,6, 7100)                  # time for the temporal
                                                    # evolution of the 
                                                    # magnetization vector

Deltat    = (time[1] - time[0])                     # time-step of evolution

''' ========================================================================

  INITIAL PARAMETERS FOR THE ADDITION OF THERMAL + ELECTROMAGNETIC NOISES

========================================================================== ''' 


alpha_n = gamma                                     # coupling constan

time_n   = np.linspace(0,6,ntime_fitz)              # time for the temporal
                                                    # evolution of the 
                                                    # magnetization vector
                                                    
##---------------------------------------------------------------------------- 
###  uncomment if the choice is the hindmarshrose model
                                               
#time_n   = np.linspace(0,6,ntime_hind)             # time for the temporal
                                                    # evolution of the 
                                                    # magnetization vector
##---------------------------------------------------------------------------- 


Deltat_n    = (time_n[1] - time_n[0])               # time-step of evolution



##---------------------------------------------------------------------------- 




