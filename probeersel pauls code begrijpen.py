# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 13:30:26 2024

@author: sven
"""
import numpy as np
import matplotlib.pyplot as plt

def gridpoints(R0, Rmax, Nr, Ngamma):
    """
    Takes:
    R0: minimal radius
    Rmax: maximum radius
    Nr: amount of radius points to take
    Ngamma: amount of angles to take.

    Returns nested array, each array is an angle and contains all radii. 
    """
    gamma = np.linspace(0, 2*np.pi, Ngamma) #array for angles
    R = np.linspace(R0, Rmax, Nr) #array for radii
    rr, gg = np.meshgrid(R,gamma) #combine
    points = rr*(np.cos(gg) - 1j*np.sin(gg)) 
    return points.ravel() #flatten your output into a 1D array. 

def complex_potential(gamma, z):
    """
    Calculate complex potential for U = 1. 
    gamma = vorticity
    mu = viscocity
    z = complex coordinate. 
    """
    mu = - 2*np.pi*(np.sqrt(z.real**2 + z.imag**2))
    return z - 1j*gamma*np.log(z)/(2*np.pi) + mu/(2*np.pi*z)

g = gridpoints(1.12, 2, 500, 500)
# plt.plot(g.real, g.imag)


clist = np.zeros(len(g))
for i in range(len(g)): 
    clist[i] = complex_potential(0, g[i])
# plt.scatter(g, clist)
    
print(g[0])

# g.real, g.imag = box_df.x, box_df.y
  
# make scalar field



# Plotting scalar field with tricontour
plt.tricontourf(g.imag, g.real, clist)

plt.title('Scalar Field')
  
# Show plot with gird
plt.grid()