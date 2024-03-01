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
    return rr*(np.cos(gg) - 1j*np.sin(gg)) #calculate points and return

# TEST. Moet 4 arrays geven (voor de 4 hoeken), met 2 elementen op R = 0 (alles 0) en R = 1.
#print(gridpoints(0, 1, 2, 4)) 

def complex_potential(U, gamma, mu, z):
    