import numpy as np
import matplotlib.pyplot as plt
import Joukowski

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
    points = rr*(np.cos(gg) + np.sin(gg)*1j) 
    print(points)
    print(points.ravel())
    return points.ravel() #flatten your output into a 1D array. 

def complex_potential(gamma, z):
    """
    Calculate complex potential for U = 1. 
    gamma = vorticity
    mu = viscocity
    z = complex coordinate. 
    """
    mu = - 2*np.pi*(np.sqrt(z.real**2 + z.imag**2))
    return z - gamma*np.log(z)/(2*np.pi) + mu/(2*np.pi*z)


g = gridpoints(1.12, 2, 250, 250) 


clist = np.zeros(len(g))
for i in range(len(g)): 
    clist[i] = complex_potential(-3, g[i])

# print(g[0])

a = Joukowski.circle(complex(0,0), 1.12, 100)
 

# Plotting scalar field with tricontour
plt.tricontourf(g.imag -0.1, g.real +0.22, clist)
plt.plot(a.real - 0.1, a.imag + 0.22)
plt.title('Scalar Field')
  
# # Show plot with gird
# plt.grid()

