#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 10:50:21 2021

@author: nehabinish
"""

import numpy as np



''' ========================================================================

                    RUNGE KUTTA FOURTH ORDER INTEGRATOR

========================================================================== '''

def RK4(f, x0, t):

    """
        Runge Kutta fourth order integrator

        |  :param
        |-----------
        |      f  : A function callable(x, t)
        |           computes the derivative of x at t
        |      x0 : array
        |           Initial condition on x (can be a vector).
        |      t  : array
        |           A sequence of time points for which to solve for y.
        |           The initial value point should be the first element
        |           of this sequence. This sequence must be monotonically
        |           increasing or monotonically decreasing; repeated values
        |           are allowed.
        |
        |  :return
        |-----------
        |     X  : array, shape (len(t), len(x0))
        |          Array containing the value of y for each desired time in t,
        |          with the initial value x0 in the first row.

    """


    dt   = t[2] - t[1]                                     # time step
    N    = len(t)                                         # length of t
    X    = np.empty((len(t), len(x0)))
    X[0] = x0

    for i in range(1, N):
        
        # Runge-Kutta fourth order algorithm

        k1   = f(X[i-1], t[i-1])
        k2   = f(X[i-1] + dt/2*k1, t[i-1] + dt/2)
        k3   = f(X[i-1] + dt/2*k2, t[i-1] + dt/2)
        k4   = f(X[i-1] + dt*k3, t[i-1] + dt)

        X[i] = X[i-1] + dt/6*(k1 + 2*k2 + 2*k3 + k4)

    return X
