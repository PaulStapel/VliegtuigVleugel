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
q = g.real*1j + g.imag ### Switch de real en imaginary waardes omdat ze op een of andere manier verkeerd om staan?

clist = np.zeros(len(g)) ### Genereer een lijst die gevuld gaat worden met waardes van de complexe potential voor coordinaten
for i in range(len(g)): 
    clist[i] = complex_potential(-3, g[i])
    print(clist[i]) 

print(clist)

a = Joukowski.circle(complex(0,0), 1.12, 100)  ### Cirkel

# qlist = Joukowski.joukowski(q)


# Plotting scalar field with tricontour
plt.tricontourf(g.real -0.1, g.imag +0.22, clist) ### Plotting the complex potential
plt.plot(a.real - 0.1, a.imag + 0.22) ### Plotting the circle
plt.title('Scalar Field')
  
# # Show plot with gird
plt.grid()
plt.show()

