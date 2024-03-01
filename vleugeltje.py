import numpy as np
import matplotlib.pyplot as plt
from Joukowski import joukowski

def gridpoints(R0, Rmax, Nr, Ngamma):
    gamma = np.linspace(0, 2*np.pi, Ngamma)
    R = np.linspace(R0, Rmax, Nr)
    rr, gg = np.meshgrid(R,gamma)
    return rr*(np.cos(gg) - 1j*np.sin(gg))

print(gridpoints(0, 1, 2, 4)) #Moet 4 arrays geven (voor de 4 hoeken), met 2 elementen op R = 0 (alles 0) en R = 1.


