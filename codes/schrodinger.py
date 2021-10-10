#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 12:10:13 2021

@author: nehabinish
"""


import numpy as np
import init
from scipy import linalg as LA2



''' ========================================================================

                DENSITY MATRIX FOR SCHRÖDINGER EQUATION.

========================================================================== ''' 

def Densitymatrix( psi ):
    
       
    """
    
      Function to generate density matrix |psi_n> <psi_n|
        
        |
        |  :param
        |-----------
        |
        |      psi  : array - size [ 2, len(t) ]
        |             Wave Functions
        |                      
        |  :return
        |-----------
        |
        |          : array - size [ N, 2, len(t) ]
        |            Density Matrix
        |
        
    """
    
    rho = []

    for n in range( init.ntime_sch ):
        
        val = 0
        
        for k in range( 0, init.N ):
                
            a   = np.matrix( psi[k, :, n] ).T
            b   = a.getH()
            val += np.dot( a,b) 
            
        rho.append(val)
    
    rho = 1/init.N * ( np.array(rho) )     
    
    return rho


''' ========================================================================

                        SCHRÖDINGER EQUATION

========================================================================== '''


# Time dependent Hamiltonian 

def schr(f_noise):
    
          
    """
    
      Evolution of the density matrix of the spin ensemble
      in the presence of electromagnetic noises.
        
        |
        |  :param
        |-----------
        |
        |      f_noise  : array - size [ N, len(t), 2 ]
        |                 Chaotic variables
        |                      
        |  :return
        |-----------
        |
        |      psi_den  : array - size [ N, 2, len(t) ]
        |                 Density Matrix
        |
        
    """
    
    psi_den = []
    
    for j in range ( 0, init.N ):
        
        H = np.zeros( ( 2, 2, len(init.t_sch) ), dtype = complex )
        
        for n in range( init.ntime_sch ):
            
            S = np.array([ [np.cos(init.omega * n),   np.sin(init.omega * n)   ],
                           [np.sin(init.omega * n), - np.cos(init.omega * n)   ]  ])
           
            H[:,:,n] = ( ( init.gamma * init.B0 * init.hbar / 2 ) * S ) \
                         + init.epsilon * ( f_noise[j,n,0] * init.Sx )  \
                         + ( f_noise[j,n,1] * init.Sy )


## ------------------------------------------------------------------------ ##
### Uncomment if the choice is Hindmarsh-Rose model
#
#             H[:,:,n] = ( ( init.gamma * init.B0 * init.hbar / 2 ) * S ) \
#                         + init.epsilon * ( f_noise[j,n,0] * init.Sx )  \
#                         + ( f_noise[j,n,1] * init.Sy ) \ 
#                         + ( f_noise[j,n,2] * init.Sz )  
## ------------------------------------------------------------------------ ##

               
             
        
        #initialising a vector of psi values, each a column vector of 2 elements
        psi = np.zeros( ( 2,len( init.t_sch ) ), dtype = complex )
        
        # psi is a column vector [1,0]
        psi[0, 0] = 1/np.sqrt(2)
        psi[1, 0] = 1/np.sqrt(2)
        
        
        for i in range( init.ntime_sch - 1 ):
            
            psi[:, i+1] = np.dot( LA2.expm(-1j * (1/init.hbar) * H[:,:,i] \
                                           * init.Deltat_sch), psi[:,i])
    
    
        psi_den.append(psi)
        
    return psi_den