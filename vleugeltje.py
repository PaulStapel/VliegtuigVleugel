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
    points = np.vstack([rr.ravel(), gg.ravel()])
    return points #flatten your output into a 1D array. 


def polar_centering(x0, y0, points):
    r = points[0]
    theta = points[1]

    # Convert polar coordinates to Cartesian coordinates
    x = r * np.cos(theta)
    y = r * np.sin(theta)

    # Centering in Cartesian coordinates
    new_x = x + x0
    new_y = y + y0

    # Convert back to polar coordinates
    new_r = np.sqrt(new_x**2 + new_y**2)
    new_theta = np.arctan2(new_y, new_x)

    return [new_r, new_theta]

def complex_potential(gamma, r, theta):
    """
    Calculate complex potential for U = 1. 
    gamma = vorticity
    mu = viscocity
    z = complex coordinate. 
    """
    return (r- (1.22/r)**2)*np.sin(theta) - gamma*np.log(r)/(2*np.pi)

def plot_stream_function(points, gamma):
    r = points[0]
    theta = points[1]

    x = r*np.cos(theta)
    y = r*np.sin(theta)

    potential = complex_potential(gamma, r, theta)
    streamfunction = potential

    figure, axes = plt.subplots() 
    cc = plt.Circle(( -0.1, 0.22 ), 1.22 ,color='black') 
 
    axes.set_aspect( 1 ) 
    axes.add_artist( cc ) 
    axes.tricontour(x, y, streamfunction, levels=20)
    # plt.scatter(-0.1,0.22) #plotting the middle of the circle
    plt.show()
    return potential


grid = gridpoints(1.12,2, 200, 200) 
grid = polar_centering(-0.1, 0.22, grid) #correction for center of circle
potential = plot_stream_function(grid, 0)

circle = Joukowski.circle(complex(-0.1,0.22), 1.12, 1000) 
wing = Joukowski.joukowski(circle)
grid_transx = Joukowski.joukowski(grid[0])
grid_transy = Joukowski.joukowski(grid[1])
pot_trans = Joukowski.joukowski(potential)

x,y = grid_transx, grid_transy
stream_trans = pot_trans
plt.tricontour(x,y, stream_trans, levels=10) 
plt.plot(wing.real, wing.imag, color='black')
plt.show()




 

