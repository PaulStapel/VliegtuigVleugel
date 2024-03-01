import numpy as np
import matplotlib.pyplot as plt
from Joukowski import joukowski

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

# TEST. Moet een array geven met stralen voor alle 4 de hoeken. 
# print(gridpoints(0, 1, 2, 4))

def complex_potential(gamma, mu, z):
    """
    Calculate complex potential for U = 1. 
    gamma = vorticity
    mu = viscocity
    z = complex coordinate. 
    """
    return z - 1j*gamma*np.log(z)/(2*np.pi) + mu/(2*np.pi*z)

#comment