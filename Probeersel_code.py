import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from vleugeltje import *
import Joukowski
from matplotlib import cm


Radius = 1.12

def alternate_flowlines(gamma, alpha=0):
    circle = Joukowski.circle(complex(-0.1,+0.22), Radius, 100) 
    wing = Joukowski.joukowski(circle)

    X = np.arange(-3, 3, 0.05) 
    Y = np.arange(-3, 3, 0.05)
    x,y = np.meshgrid(X, Y)

    z = x+1j*y
    mask = np.absolute(z - (-0.1 + 0.22j)) >= Radius
    #z -= -0.1 + 0.22j 
    z = z[mask]
    x,y = x[mask], y[mask]

    w = Joukowski.joukowski(z)

    pot = complex_potential(gamma, z, alpha)
    stream = pot.imag

    plt.tricontour(x, y, stream, levels=20, color='blue')
    # Plot the circle and the transformed wing
    plt.plot(circle.real, circle.imag, label='Circle')
    plt.scatter(-0.1,0.22) #plotting the middle of the circle
    
    # Show the plot
    plt.legend()
    plt.axis('equal')
    plt.show()

    plt.tricontour(w.real, w.imag, stream, levels=20)
    plt.plot(wing.real, wing.imag, label='Joukowski Wing')
    plt.show()
