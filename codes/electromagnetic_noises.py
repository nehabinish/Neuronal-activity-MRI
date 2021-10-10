#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 11:21:23 2021

@author: nehabinish
"""

import numpy as np
import matplotlib.pyplot as plt

from integrator import RK4
from dynamical_sys import fitzhugh
import init
import schrodinger


## ------------------------------------------------------------------------ ##
#                                                  |
# plt.rc('xtick', labelsize= 16)                   | uncomment to increase 
# plt.rc('ytick', labelsize= 16)                   | the axis tick size and
# plt.rcParams.update({'font.size': 18})           | increase font size
#                                                  | for plot labels.
## ------------------------------------------------------------------------ ##




''' ========================================================================

         GENERATION OF CHAOTIC VARAIBLES FROM CELL ELECTRIC ACTIVITY

========================================================================== '''


chaotic_var = []

for k in range ( 0, init.N ):
    
## ------------------------------------------------------------------------ ##
#
### Uncomment if the choice is hindmarsh-rose model.
 
    # rand1 = np.random.uniform(-1,1) 
    # rand2 = np.random.uniform(-1,1) 
    # rand3 = np.random.uniform(-1,1) 
    # chaotic_var.append(np.array([rand1,rand2]))
    
## ------------------------------------------------------------------------ ##
   
    
    rand1 = np.random.uniform(-1,2) 
    rand2 = np.random.uniform(-1,2) 
 
    chaotic_var.append(np.array([rand1,rand2]))
    
 


f_xyz = []
   
for x0 in chaotic_var:
     
    # chaotic noise variables 
    
## ------------------------------------------------------------------------ ##
### Uncomment if the choice is hindmarsh-rose model.
    # f_vw.append(RK4(hindmarsh, x0, t_hind))  
## ------------------------------------------------------------------------ ##
   
    f_xyz.append( RK4(fitzhugh, x0, init.t_fitz))
 
f_xyz = np.array(f_xyz) 


## ------------------------------------------------------------------------ ##
#  
###  uncomment to generate the time series and phase space plots of the 
###  chosen classical dynamical system.
#
# for k in range ( 0, init.N ):
#                                               
#     plt.figure(figsize=(10,5))            
#     plt.plot(init.t_fitz, f_vw[k,:,0])    
#     plt.title("time series")              
#     plt.ylabel('v')                         
#     plt.xlabel('t')                         
#     plt.xlim(0,10)                        
#     plt.grid()                            
#                                               
#     plt.figure(figsize=(10,5))            
#     plt.plot(f_vw[k,:,0], f_vw[k,:,1])    
#     plt.title("Phase portrait")           
#     plt.ylabel('x')                       
#     plt.xlabel('y')                       
#     plt.grid()                            
#                                           
## ------------------------------------------------------------------------ ##

''' ========================================================================

                    GENERATION OF THE DENSITY MATRIXX

========================================================================== '''

densitymatrix = schrodinger.schr(f_xyz)

  
''' ========================================================================

            DYNAMICS OF THE SPIN ENSEMBLE IN THE PRESENCE
                    OF ELECTROMAGNETIC NOISES.

========================================================================== '''  
    
densitymatrix = np.array(densitymatrix)      
rho_t         = schrodinger.Densitymatrix(densitymatrix)


plt.figure(figsize=(10,5))
plt.plot(init.t_sch, rho_t[:,0,0], label =' $ \left< \u2191 | \\rho | \u2191 \\right> $' )
plt.plot(init.t_sch, rho_t[:,1,1], label = '$ \left< \u2193 | \\rho | \u2193 \\right> $')
plt.ylabel('survival probability')
plt.xlabel(' time( a.u.) ')
plt.title('Evolution of the density matrix with N = {}'.format(init.N))
plt.grid()
plt.legend()


plt.figure(figsize=(10,5))
plt.plot(init.t_sch, np.abs(rho_t[:,0,1]), label = '$ \left|\left< \u2191 | \\rho | \u2193 \\right>  \\right| $')
plt.ylabel('coherence')
plt.xlabel(' time (a.u.) ')
plt.title('Evolution of the density matrix with N = {}'.format(init.N))
plt.grid()
plt.legend()


 



    