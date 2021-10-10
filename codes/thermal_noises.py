#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 20 12:29:51 2021

@author: nehabinish
"""


import numpy as np
import matplotlib.pyplot as plt
from functions_thermal import control
from functions_thermal import integratemat
import init


## ------------------------------------------------------------------------ ##
#                                                  |
# plt.rc('xtick', labelsize= 16)                   | uncomment to increase 
# plt.rc('ytick', labelsize= 16)                   | the axis tick size and
# plt.rcParams.update({'font.size': 18})           | increase font size
#                                                  | for plot labels.
## ------------------------------------------------------------------------ ##


''' ========================================================================

                                 GENERATE CONTROL FIELDS

========================================================================== ''' 


omegax, omegay = control( init.time )


''' ========================================================================

    EVOLUTION OF THE MAGNETIC FIELD IN THE PRESENCE OF THERMAL NOISES.

========================================================================== ''' 


M  = integratemat( init.M_init, init.time, omegax, omegay )
 

''' ========================================================================

                    PLOT THE EVOLUTION OF THE MAGNETIZATION VECTOR

========================================================================== ''' 

   

plt.figure()
ax = plt.axes(projection='3d')

# Data for a three-dimensional line
zline = M[100:,3]
xline = M[100:,2]
yline = M[100:,1]

ax.set_xlabel('Mx', labelpad=20)
ax.set_ylabel('My', labelpad=20 )
ax.set_zlabel('Mz', labelpad=20 )
plt.title('Evolution of the Magnetization.', fontsize = 18 )
ax.plot3D(xline, yline, zline, 'red')


plt.figure(figsize=(10,5))
plt.plot(init.time, M[:,1], 'r', label = 'Mx')
plt.plot(init.time, M[:,2], 'k', label = 'My')
plt.ylabel('Magnetization')
plt.xlabel(' time (s) ')
plt.title('Evolution of the Magnetization.')
plt.grid()
plt.legend()


plt.figure(figsize=(10,5))
plt.plot(init.time, M[:,3], 'k', label = 'Mz')
plt.ylabel('Magnetization')
plt.xlabel(' time (s) ')
plt.title('Evolution of the Magnetization.')
plt.grid()
plt.legend()

plt.show()


