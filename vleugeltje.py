import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
# import Joukowski
from matplotlib import cm
from scipy import constants as con


def circle(center, R, points): ### Circle function
    """
    Takes
    R: Radius
    center: center (can be complex)
    points: amount of points to plot the circle
    
    Returns complex circle coordinates.
    The real part is x, the imaginary part is y.
    """ 
    theta = np.linspace(0, 2*np.pi, points) #define angles
    circle = R* (np.cos(theta) + 1j*np.sin(theta))+ center #calculate points on circle
    return circle

def joukowski(circle): ### Jouwkowski function
    """ 
    Takes
    Circle
    
    Returns transformed jouwkowski coordinates.
    The real part is x, the imaginary part is y.
    """
    jk = circle+ 1/circle # calculate joukowski transform 
    return jk


def deg2rad(deg): ### degrees to pi
    return deg*np.pi/180 

def gridpoints(R0, Rmax, Nr, Ngamma): ### Generating gridpoints
    """
    Takes:
    R0: minimal radius
    Rmax: maximum radius
    Nr: amount of radius points to take
    Ngamma: amount of angles to take.

    Returns array containing all coordinates as complex numbers. 
    """
    gamma = np.linspace(0, 2*np.pi, Ngamma) #array for angles
    R = np.linspace(R0, Rmax, Nr) #array for radii
    rr, gg = np.meshgrid(R,gamma) #combine
    points = rr*(np.cos(gg) + 1j*np.sin(gg)) #Transform to cartesian coordinates
    return points.ravel()

def complex_centering(x0,y0,points):  
    return points + x0 + 1j*y0 #center the points to some coordinate (x0,y0)

def complex_potential(gamma, z, alpha): ###
    """
    Calculate complex potential for U = 1. 
    1st term: parallel flow
    2nd term: dipole
    3rd term: vortex with vorticity gamma
    z = complex coordinate. 
    alpha = angle
    """
    #Correct for angle of attack with terms as stated in excercise 1
    alpha1 = np.cos(alpha) - 1j*np.sin(alpha) #parallel angle term
    alpha2 = np.cos(alpha) + 1j*np.sin(alpha) #dipole angle term

    return z*alpha1 + 1.12**2*alpha2/z - 1j*gamma/(2*np.pi)*np.log(z) #for formula see reader

def find_tail(joukowski_wing):
    """Returns the boundaries of the local frame where the tail can be found"""
    xmax_wing,xmin_wing = max(joukowski_wing.real), min(joukowski_wing.real)
    ymax_wing,ymin_wing = max(joukowski_wing.imag), min(joukowski_wing.imag)
    
    xlength,ylength = (xmax_wing - xmin_wing), (ymax_wing - ymin_wing)

    max_index = np.argmax(joukowski_wing.real)
    ypoint = joukowski_wing.imag[max_index]
    
    #take small differences around the wing to find a small area. 
    xmax,xmin = xmax_wing + 0.1*xlength, xmax_wing - 0.1*xlength
    ymax,ymin = ypoint + 0.1*ylength,  ypoint - 0.1*ylength

    return (xmin,xmax,ymin,ymax)

def zoomplot(joukowski_wing, x,y,z, gamma):
    """Create a zoomed in plot of the tail"""
    
    xmin, xmax, ymin, ymax = find_tail(joukowski_wing) # find a small frame for the tail  
    
    plt.plot()
    plt.title(f'Zoomplot of tail, gamma = {gamma}.')
    plt.xlabel("Length in m.")
    plt.ylabel("Length in m.")
    plt.tricontour(x,y,z, levels = 240) ### Added more levels so you can see the streamfunction
    
    plt.plot(joukowski_wing.real,joukowski_wing.imag, color='black')
    plt.axis([xmin,xmax,ymin,ymax])
    plt.show()

def extract_potential(gamma, alpha, Nr, Ngamma):
    """Calculate coordinates and respective streamfunction for
    gamma = voriticity
    alpha = angle of attack
    N = amount of points per axis of the grid
    """
    grid = gridpoints(radius, 3*radius, Nr, Ngamma) # Take gridpoints to consider
    potential = complex_potential(gamma, grid, alpha) # Calculate potential for all points

    z = complex_centering(x0,y0, grid)#Center the x and y-axis to the cylinder. 
    z_trans = joukowski(z) #transform coordinates

    return z, z_trans, potential

def plot_streamfunction(z,z_trans,potential, gamma):
    """Graphically visualise the streamfunctions for the circle and wing."""

    # circle = circle(complex(x0,y0), radius, points) #Create a circle
    # wing = joukowski(circle) #Create a wing by transforming the circle
    
    
    ### Defining the x and y coordinates for the streamfunction
    x,y = z.real, z.imag 
    x_trans,y_trans = z_trans.real, z_trans.imag
    
    ### The streamfunction is the imaginary part of the complex potential
    streamfunction = potential.imag 
    
    #plot cylinder streamfunction
    plt.title(f'The streamfunction around the cylinder with gamma = {gamma}.')
    plt.xlabel("Length in m.")
    plt.ylabel("Length in m.")
    plt.tricontour(x, y, streamfunction, levels=20, color='black')
    plt.plot(circle.real, circle.imag, color='black')
    
    plt.scatter(x0, y0) #plotting the middle of the circle
    plt.show()

    # Plot wing streamfunction
    plt.title(f'The streamfunction around the wing with gamma = {gamma}.')
    plt.xlabel("Length in m.")
    plt.ylabel("Length in m.")
    plt.tricontour(x_trans,y_trans, streamfunction, levels=20)
    plt.plot(wing.real, wing.imag, color='black')
    plt.show()

    #Plot zoomed-in tail
    zoomplot(wing, x_trans, y_trans, streamfunction,gamma)

def get_velocity(z, potential, conjugate=False):
    """Create a list of velocities for every complex coordinate from the potential function"""
    z_new = z[:-1] # Making the z list one shorter because you lose 1 point with the approximation method
    speed = np.zeros(len(z_new), dtype = np.complex_) # list to store speeds in

    for i in range(len(z_new)):

        speed[i] = (potential[i] - potential[i + 1])/(z[i] - z[i + 1])
        if conjugate: 
            speed[i] = np.conj(speed[i])

    return z_new, speed.real, speed.imag

def vel_field(gamma, alpha):
    """From the velocities around the wing and circle, visualise the velocity-field."""

    z, z_trans, potential = extract_potential(gamma,alpha, 10, 35)

    ### Plotting the velocity field for the cylinder
    z,u,v = get_velocity(z, potential, conjugate=True) 
    x,y = z.real, z.imag
    plt.figure()
    plt.title(f'Velocity around a cylinder for gamma = {gamma}.')
    plt.xlabel("Length in m.")
    plt.ylabel("Length in m.")
    plt.plot(circle.real, circle.imag, color='black')
    plt.quiver(x,y,u,v)
    plt.show() 

    ### Plotting the velocity field for the wing
    z_trans,wu,wv = get_velocity(z_trans, potential, conjugate=True)
    wx,wy = z_trans.real, z_trans.imag
    plt.figure()
    plt.title(f'Velocity around a wing for gamma = {gamma}.')
    plt.xlabel("Length in m.")
    plt.ylabel("Length in m.")
    plt.quiver(wx,wy,wu,wv)
    plt.plot(wing.real, wing.imag, color='black')
    plt.show()
    

def find_lift(z, wing, speed, N):
    """Calculate the lift around some wing-like object"""
    circle_speed = np.zeros(N, dtype =np.complex128)
    e = np.zeros(N, dtype =np.complex128)
    for i in range (N):
        a = abs(z - wing[i]) 
        b = min(a) ### Minimal value between the gridpoints and the circle point
        c = np.where(a == b)[0] ### indexing where that value is
        d = int(c[0])
        e[i] = z[d] ### Saving the gridpoints
        circle_speed[i] = speed[d] ### saving the corresponding pressure

    dz = np.roll(e,-1)/2 - np.roll(e,1)/2 ### Dz value 
    Fx_min_Fy = 1.225j/2*sum(circle_speed**2 * dz) ## multiplying the pressure with 10 and dz to make it force by multiplying it with the area
    
    Fy = -Fx_min_Fy.imag
    Fx = Fx_min_Fy.real

    return Fx, Fy


def pressure_field(gamma, alpha, plot_pressure=False):
    """For the desired shapes, calculate the pressure and lift. Visualise the lift in a contour plot."""

    z, z_trans, potential = extract_potential(gamma, alpha, 200, 200)

    N = len(circle)
    H = 1**2/2  + 101325/1.225 ### Im leaving the gravitational constant out

    #Lift for cylinder
    new_z, u ,v = get_velocity(z,potential)   
    speed = u + 1j*v 
    Fx_circle, Fy_circle = find_lift(z, circle, speed, N)

    #Lift for wing
    new_z_trans, wu ,wv = get_velocity(z_trans,potential)
    wspeed = wu + 1j*wv
    Fx_wing, Fy_wing = find_lift(z_trans, wing, wspeed, N)
    
    #plot
    if plot_pressure:
        x,y = new_z.real, new_z.imag
        pressure = (H - (u**2 + v**2)/2)*1.225 
        plt.figure()
        plt.title(f'Pressure cylinder, gamma = {gamma}.')
        plt.xlabel("Length in m.")
        plt.ylabel("Length in m.")
        plt.tricontourf(x, y, pressure, levels=500, color='black')
        plt.colorbar(label = 'pressure')
        plt.plot(circle.real, circle.imag, color='black')
        plt.show()
        
        x_trans,y_trans = new_z_trans.real, new_z_trans.imag
        wpressure = (H - (wu**2 + wv**2)/2)*1.225 
        plt.figure()
        plt.title(f'Pressure wing, gamma = {gamma}.')
        plt.xlabel("Length in m.")
        plt.ylabel("Length in m.")
        plt.tricontourf(x_trans, y_trans, wpressure, levels=500, color='black')
        plt.legend()
        plt.colorbar(label = 'pressure')
        plt.plot(wing.real, wing.imag, color='black')
        plt.show()
    
    #print values
    print(f' simulation upward lift cylinder {Fy_circle:.3f}, dragg {Fx_circle:.3f}')
    print(f' simulation upward lift wing {Fy_wing:.3f},  dragg {Fx_wing:.3f}')
    print(f'analytical wing {-1.225 * 1 * gamma:.3f}')

#define circle, global variables
x0 = -0.1
y0 = 0.22
radius = 1.12
points = 1000 


circle = (circle(complex(x0,y0), radius, points)) ### Creating the circle
wing = joukowski(circle ) ### Transforming the wing


### Plot for the jouwkowski wing                                                                                       
plt.plot(wing.real, wing.imag) 
plt.xlim(-2.5,2.5)  
plt.ylim(-2.5,2.5)
plt.title(f'Wing in {points} points.')
plt.xlabel("Length in m.")
plt.ylabel("Length in m.")
plt.show()


#Create vector field parameters
gamma1,gamma2, gamma3 = 0, -3, round(-np.e,3) ### Different gammas 
alpha = deg2rad(0) # Angle of Attack


z, z_trans, potential = extract_potential(gamma1, alpha, 100, 100) ### defining the coordinates, z.
### the wing coordinates, z_trans,  and the potential, for a certain gamma and alpha
plot_streamfunction(z,z_trans,potential, gamma1) ### Plotting the streamfunction for a certain gamma and alpha

z, z_trans, potential = extract_potential(gamma2, alpha, 100, 100)
plot_streamfunction(z,z_trans,potential, gamma2)

z, z_trans, potential = extract_potential(gamma3, alpha, 100, 100)
plot_streamfunction(z,z_trans,potential, gamma3)

vel_field(gamma3, alpha) ### Plotting the velocity field for the gamma with smooth trainling edge

pressure_field(gamma3, alpha, plot_pressure=True)### Plotting the pressure for gamma


 
