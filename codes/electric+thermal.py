#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 27 10:02:55 2021

@author: nehabinish
"""

import numpy as np
import matplotlib.pyplot as plt

import init

from dynamical_sys import fitzhugh
from integrator import RK4
from functions_thermal import control
from functions_thermal import integratemat_noises

from scipy import linalg as LA2



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



''' ========================================================================

                                 GENERATE CONTROL FIELDS

========================================================================== ''' 


omegax, omegay = control( init.time_n)



''' ========================================================================

            EVOLUTION OF THE MAGNETIC FIELD IN THE PRESENCE OF 
                    THERMAL AND ELECTROMAGNETIC NOISES.

========================================================================== ''' 


M  = integratemat_noises( init.M_init, init.time_n, omegax, omegay, f_xyz )
 

''' ========================================================================

                     PLOTS FOR THE  EVOLUTION OF MAGNETIZATION

========================================================================== ''' 

    
plt.figure()
ax = plt.axes(projection='3d')

# Data for a three-dimensional line
zline = M[:,3]
xline = M[:,2]
yline = M[:,1]

ax.set_xlabel('Mx', labelpad = 15)
ax.set_ylabel('My', labelpad = 15)
ax.set_zlabel('Mz', labelpad = 15)
plt.title('Evolution of the Magnetization.')
ax.plot3D(xline, yline, zline, 'red')


plt.figure(figsize=(10,5))
plt.plot(init.time_n, M[:,1], 'r', label = 'Mx')
plt.plot(init.time_n, M[:,2], 'k', label = 'My')
plt.ylabel('Magnetization')
plt.xlabel(' time(a.u.) ')
plt.title('Evolution of the Magnetization.')
plt.grid()
plt.legend()
 
plt.figure(figsize=(10,5))
plt.plot(init.time_n, M[:,3], 'k',  label = 'Mz')
plt.ylabel('Magnetization')
plt.xlabel(' time (a.u.) ')
plt.title('Evolution of the Magnetization.')
plt.grid()
plt.legend()

plt.show()

