import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import Joukowski
from matplotlib import cm

def deg2rad(deg):
    return deg*np.pi/180

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
    return points.ravel()

def complex_centering(x0,y0,points):
    return points + x0 + 1j*y0 

def complex_potential(gamma, z, alpha):
    """
    Calculate complex potential for U = 1. 
    gamma = vorticity
    mu = viscocity
    z = complex coordinate. 
    """
    alpha1 = np.cos(alpha) - 1j*np.sin(alpha)
    alpha2 = np.cos(alpha) + 1j*np.sin(alpha)
    return (z*alpha1 + 1.12**2*alpha2/z) - 1j*gamma/(2*np.pi)*np.log(z)

def plot_stream_function(points, gamma):
    x = points.real - 0.1
    y = points.imag + 0.22

    potential = complex_potential(gamma, points, alpha)
    streamfunction = potential.imag

    figure, axes = plt.subplots() 
    cc = plt.Circle(( -0.1, 0.22 ), 1.12 ,color='black') 
 
    axes.set_aspect( 1 ) 
    axes.add_artist( cc ) 
    axes.tricontour(x, y, streamfunction, levels=20)
    plt.scatter(-0.1,0.22) #plotting the middle of the circle
    plt.show()
    return potential


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



alpha = deg2rad(-10) #angle of attack.
radius = 1.12


grid = gridpoints(1.12 ,2.5, 200, 100) #Generate gridpoints

potential = plot_stream_function(grid, -2)

grid = complex_centering(-0.1,0.22, grid) 

circle = Joukowski.circle(complex(-0.1,0.22), 1.12, 1000) 
wing = Joukowski.joukowski(circle)
grid_trans = Joukowski.joukowski(grid)

x,y = grid_trans.real, grid_trans.imag

plt.tricontour(x,y, potential.imag, levels=20) #bug
plt.plot(wing.real, wing.imag, color='black')
plt.show()




 

