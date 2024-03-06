import numpy as np
import matplotlib.pyplot as plt
import Joukowski
from scipy.interpolate import griddata
from matplotlib import cm

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
    points = rr*(np.cos(gg) + 1j*np.sin(gg)) 
    return points.ravel() #flatten your output into a 1D array. 

# TEST. Moet een array geven met stralen voor alle 4 de hoeken. 
# print(gridpoints(0, 1, 2, 4))

def complex_potential(gamma, z):
    """
    Calculate complex potential for U = 1. 
    gamma = vorticity
    mu = viscocity
    z = complex coordinate. 
    """
    mu = - 2*np.pi*(np.sqrt(z.real**2 + z.imag**2))
    return z - 1j*gamma*np.log(z)/(2*np.pi) - mu/(2*np.pi*z)


def plot_stream_function(points, gamma):
    x = points.real
    y = points.imag

    potential = complex_potential(gamma, points)
    streamfunction = potential.imag

    figure, axes = plt.subplots() 
    cc = plt.Circle(( -0.1, 0.22 ), 1.22 ,color='black') 
 
    axes.set_aspect( 1 ) 
    axes.add_artist( cc ) 
    axes.tricontour(x, y, streamfunction, levels=20)
    # plt.scatter(-0.1,0.22) #plotting the middle of the circle
    plt.show()

grid = gridpoints(1.12,2, 50, 50) - 0.1 + 0.22j #correction for center of circle
plot_stream_function(grid, 0)

 

