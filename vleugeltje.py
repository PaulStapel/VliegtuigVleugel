import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import Joukowski
from matplotlib import cm
from scipy import constants as con

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
    points = rr*(np.cos(gg) + 1j*np.sin(gg)) #Transform to cartesian coordinates
    return points.ravel()

def complex_centering(x0,y0,points):
    return points + x0 + 1j*y0 #center the points to analyse

def complex_potential(gamma, z, alpha):
    """
    Calculate complex potential for U = 1. 
    1st term: parallel flow
    2nd term: dipole
    3rd term: vortex with vorticity gamma
    z = complex coordinate. 
    """
    #Correct for angle of attack with terms as stated in the excercise
    alpha1 = np.cos(alpha) - 1j*np.sin(alpha) #first angle term
    alpha2 = np.cos(alpha) + 1j*np.sin(alpha) #second angle term

    return (z*alpha1 + 1.12**2*alpha2/z) - 1j*gamma/(2*np.pi)*np.log(z) #for formula see reader

def zoomplotje(joukowski_wing, x,y,z):
    xmax_wing,xmin_wing = max(joukowski_wing.real), min(joukowski_wing.real)
    ymax_wing,ymin_wing = max(joukowski_wing.imag), min(joukowski_wing.imag)
    
    xl,yl = (xmax_wing - xmin_wing), (ymax_wing - ymin_wing)

    max_index = np.argmax(joukowski_wing.real)
    ypoint = joukowski_wing.imag[max_index]
    
    xmax,xmin = xmax_wing + 0.1*xl, xmax_wing - 0.1*xl
    ymax,ymin = ypoint + 0.1*yl,  ypoint - 0.1*yl
    
    plt.plot()
    plt.tricontour(x,y,z, levels = 240) ### Added more levels so you can see the streamfunction
    ### in the zoomplot as well
    plt.plot(joukowski_wing.real,joukowski_wing.imag, color='black')
    plt.axis([xmin,xmax,ymin,ymax])
    plt.show()

def extract_streamfunction(x0, y0, radius, gamma, alpha):
    circle = Joukowski.circle(complex(x0,y0), radius, 1000) #Create a circle
    wing = Joukowski.joukowski(circle) #Create a wing by transforming the circle

    grid = gridpoints(radius, 3*radius, 200, 100) # Take gridpoints to consider
    potential = complex_potential(gamma, grid, alpha) # Calculate potential for all points
    streamfunction = potential.imag # Straemfunction is imaginary term of potential 

    z = complex_centering(x0,y0, grid)#Center the x and y-axis to the cylinder. 
    x,y = z.real, z.imag 


    z_trans = Joukowski.joukowski(z) #transform coordinates
    x_trans,y_trans = z_trans.real, z_trans.imag

    # Plot cylinder streamfunction
    
    streams = plt.tricontour(x, y, streamfunction, levels=20, color='black')
    plt.plot(circle.real, circle.imag, color='black')
    plt.scatter(x0, y0) #plotting the middle of the circle
    plt.show()

    # Plot wing streamfunction
    trans_streams = plt.tricontour(x_trans,y_trans, streamfunction, levels=20) #bug
    plt.plot(wing.real, wing.imag, color='black')
    plt.show()

    #Plot zoom-in
    zoomplotje(wing, x_trans, y_trans, streamfunction)

    return z, streams, z_trans, trans_streams

def vel_field(x0, y0, radius, gamma, alpha):
    circle = Joukowski.circle(complex(x0,y0), radius, 1000) #Create a circle
    wing = Joukowski.joukowski(circle) #Create a wing by transforming the circle

    grid = gridpoints(radius, 3*radius, 10, 35) # Take gridpoints to consider
    potential = complex_potential(gamma, grid, alpha) # Calculate potential for all points
    
    z = complex_centering(x0,y0, grid)#Center the x and y-axis to the cylinder.
    new_z = z[:-1] ### Making the z list one shorter because you lose 1 point with the approximation method
    x,y = new_z.real, new_z.imag ### x,y coordinates
    speed = np.zeros(len(new_z), dtype = np.complex_) ### Empty lists

    for i in range(len(new_z)):

        speed[i] = (potential[i] - potential[i + 1])/(z[i] - z[i + 1])
        speed[i] = np.conj( speed[i])
        
    plt.figure()
    u,v = speed.real, speed.imag
    plt.plot(circle.real, circle.imag, color='black')
    plt.quiver(x,y,u,v)
    plt.show()
    
    new_z = Joukowski.joukowski(new_z) #transform coordinates    
    wz = Joukowski.joukowski(z) #transform coordinates
    new_wz = wz[: - 1]
    wx,wy = new_wz.real, new_wz.imag     
    wspeed = np.zeros(len(new_wz), dtype = np.complex_)

    
    for i in range(len(new_wz)):
        ws = (potential[i] - potential[i + 1])/(wz[i] - wz[i + 1])
        wspeed[i] = np.conj( ws)

    plt.figure()
    wu,wv = wspeed.real, wspeed.imag
    plt.quiver(wx,wy,wu,wv)
    plt.plot(wing.real, wing.imag, color='black')
    plt.show()
    

def pressure_field(x0, y0, radius, gamma, alpha):
    circle = Joukowski.circle(complex(x0,y0), radius, 1000) #Create a circle
    wing = Joukowski.joukowski(circle) #Create a wing by transforming the circle

    grid = gridpoints(radius, 3*radius, 400, 400) # Take gridpoints to consider
    potential = complex_potential(gamma, grid, alpha) # Calculate potential for all points
    
    H = 1**2/2  + 101325/1.225 ### Im leaving the gravitational constant out

    z = complex_centering(x0,y0, grid)#Center the x and y-axis to the cylinder.
    new_z = z[:-1] ### Making the z list one shorter because you lose 1 point with the approximation method
    x,y = new_z.real, new_z.imag ### x,y coordinates
    
    speed = np.zeros(len(new_z), dtype = np.complex_) ### Empty lists
    pressure = np.zeros(len(speed), dtype = np.complex_)
    
    for i in range(len(new_z)):

        speed[i] = (potential[i] - potential[i + 1])/(z[i] - z[i + 1])
        speed[i] = np.conj( speed[i])
        pressure[i] = (H - (speed[i].imag**2 + speed[i].real**2)/2)*1.225 
    
    
    
    # circle_pressure= np.zeros(1000)
    # e = np.zeros(1000)
    # for i in range (1000):
    #     a = abs(grid - circle[i]) 
    #     b = min(a) ### Minimal value between the gridpoints and the circle point
    #     c = np.where(a == b)[0] ### indexing where that value is
    #     d = int(c[0])
    #     e[i] = grid[d] ### Saving the gridpoints
    #     circle_pressure[i] = pressure[d] ### saving the corresponding pressure
    
    # dz = np.zeros(1000)
    # dz = np.roll(e,-1)/2 - np.roll(e,1)/2 ### Dz value 
    # Fx_min_Fy_circle = 1.225*sum((circle_pressure * dz*10)**2) * 1j/2 ## multiplying the pressure with 10 and dz to make it force by multiplying it with the area
    
    # Fy_circle = Fx_min_Fy_circle.imag
    # Fx_circle = Fx_min_Fy_circle.real
    # print(Fy_circle,Fx_circle)
    
    
    circle_speed= np.zeros(1000, dtype =np.complex128)
    e = np.zeros(1000, dtype =np.complex128)
    for i in range (1000):
        a = abs(z - circle[i]) 
        b = min(a) ### Minimal value between the gridpoints and the circle point
        c = np.where(a == b)[0] ### indexing where that value is
     
        d = int(c[0])
        e[i] = z[d] ### Saving the gridpoints
        circle_speed[i] = speed[d] ### saving the corresponding pressure
    

    dz = np.roll(e,-1)/2 - np.roll(e,1)/2 ### Dz value 
    Fx_min_Fy_circle = 1.225j/2*sum(circle_speed**2 * dz) ## multiplying the pressure with 10 and dz to make it force by multiplying it with the area
    
    Fy_circle = -Fx_min_Fy_circle.imag
    Fx_circle = Fx_min_Fy_circle.real
    print(' simulation upward lift cylinder',Fy_circle, 'dragg',Fx_circle)
    

    

    # plt.figure()
    # plt.tricontour(x, y, pressure, levels=500, color='black')
    # plt.plot(circle.real, circle.imag, color='black')
    # plt.show()
    
    
    wz = Joukowski.joukowski(z) #transform coordinates
    new_wz = wz[: - 1]
    wx,wy = new_wz.real, new_wz.imag 
    
    wspeed = np.zeros(len(new_wz), dtype = np.complex128)
    wpressure = np.zeros(len(speed), dtype = np.complex128)
    
    for i in range(len(new_wz)):
        ws = (potential[i] - potential[i + 1])/(wz[i] - wz[i + 1])
        wspeed[i] = np.conj( ws)
        wpressure[i] = (H - (wspeed[i].imag**2 + wspeed[i].real**2)/2)*1225 
        
        
        
    wing_speed= np.zeros(1000, dtype =np.complex128)
    e = np.zeros(1000, dtype=np.complex128)
    for i in range (1000):
         a = abs(wz - wing[i]) 
         b = min(a) ### Minimal value between the gridpoints and the circle point
         c = np.where(a==b)[0] ### indexing where that value is
         d = int(c[0])
         e[i] = wz[d] ### Saving the gridpoints
         wing_speed[i] = wspeed[d] ### saving the corresponding pressure
     
    dz = np.roll(e,-1)/2 - np.roll(e,1)/2 ### Dz value 
    Fx_min_Fy_circle = 1.225j/2*sum(wing_speed**2 * dz) ## multiplying the pressure with 10 and dz to make it force by multiplying it with the area
     
    Fy_circle = -Fx_min_Fy_circle.imag
    Fx_circle = Fx_min_Fy_circle.real
    print(' simulation upward lift wing',Fy_circle, 'dragg',Fx_circle)
    print('analytical wing', -1.225 * 1 * gamma)
    
    # plt.figure()
    # plt.tricontour(wx, wy, wpressure, levels=5000, color='black')
    # plt.plot(wing.real, wing.imag, color='black')
    # plt.show()

x0 = -0.1
y0 = 0.22
radius = 1.12 
gamma = -np.e  ### Makes for a smooth flow at the trailing edge
alpha = deg2rad(0) #degrees

# z, streams, z_trans, trans_streams = extract_streamfunction(x0, y0, radius, gamma, alpha)


pressure_field(x0, y0, radius, gamma, alpha)
# vel_field(x0, y0, radius, gamma, alpha)
 