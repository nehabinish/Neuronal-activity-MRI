#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 11:14:21 2021

@author: nehabinish
"""

import numpy as np


''' ========================================================================

                            FITZHUGH NAGUMO MODEL

========================================================================== '''

def fitzhugh(x, t):
    """
    
      Time derivative of the Fitzhugh-Nagumo neural model.

        |  :param
        |-----------
        |      X      : array - size 2
        |                 [ Membrane potential, Recovery variable ]
        |      a, b   : float
        |                 Fitzhugh-Nagumo model parameters
        |      tau    : float   
        |                 Time scale. (Relaxation time)
        |      eps, I : float   
        |                 Bifurcation parameter
        |
        |  :return
        |-----------
        |        : array
        |          dv/dt, du/dt
        |          Time derivative of v and u 
        |
        |
        

    """

    alph = 10
    b = 1
    a = 0.1
    f = 0.1271
    
    omega = 2 * np.pi * f
    
    I = (a/omega) * np.cos(omega * t)
    
    vdot = ( x[0] * (x[0] - 1) * (1 - (alph * x[0])) ) - x[1] + I
    wdot = b * x[0]
    
## ------------------------------------------------------------------------ ##
#
###  uncomment to use the chaotic phase space different from the
###  chaotic attractor example used.
#   
    # a = 0.267
    # T = 2.025
    # f = 2 * np.pi / T

    # I = ( a/2 * np.pi * f ) * np.cos( 2 * np.pi * f * t )

    # # differential equations for the fitzhugh nagumo model
  
    # vdot = 10 * ( x[0] + x[1] - x[0]**3 / 3 + I )
    # wdot = x[0] - ( 0.8 * x[1] ) + 0.7

## ------------------------------------------------------------------------ ##
     
    

    return np.array([ vdot, wdot ])





''' ========================================================================

                    HINDMARSH - ROSE MODEL 

========================================================================== '''


def hindmarsh(X, t):
    
    
    """
    
      Time derivative of the Hindmarsh-Rose neural model.

        |  :param
        |-----------
        |      X           : array - size 3
        |                    [Membrane potential, Rate of fast ion channel,
        |                             Spiking variable]
        |      a, b, c, d  : float
        |                    Hindmarsh-rose model Parameters.
        |      ts, x0, b   : float   
        |                    Hindmarsh-rose model Parameters.
        |      eps, I      : float   
        |                    Bifurcation parameter
        |
        |  :return
        |-----------
        |          : array
        |          dx/dt, dy/dt, dz/dt
        |          Time derivative of x, y and z
        |
        |

    """
    
    # a  = 1
    # c  = 1.0
    # d  = 5.0
    # s  = 4.0
    # x0 = - 1.6
    
    # # Bifurcation parameters
    
    # b   = 3.09
    # I   = 3.25
    # eps = 0.15
    
    # d1  = 0.8 * np.sin(t)
    # d2  = 0.5
    # d3  = 0.5 * np.cos(t)
    
    
    x,y,z = X
    
    taux = 0.03
    tauz = 0.8
    a = 0.49
    b = d = s = 1.0
    c = 0.0322
    
    # dxdt  = y - (a * x**3) + (b * x**2) + I - z 
    # dydt  = c - (d * x**2) - y 
    # dzdt  = eps * ( (s * (x - x0)) - z) 
    
    dxdt  = 1/taux * ( y - (a * x**3) + (b * x**2) + z - (taux * x))
    dydt  = - (a * x**3) - ((d - b) * x**2) + z
    dzdt  = 1/tauz * ( -(s*x) - z + c )
    
    
    return np.array([dxdt, dydt, dzdt])

