# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 15:24:31 2024

@author: sven
"""
import numpy as np
import matplotlib.pyplot as plt



### defining the constants
rho = 1000 ### kg/m^3
mu = 1 ### Ns/m^3
T = 0.1
T2 = T/2
w = 2*np.pi/T
con = mu/rho
Uamp = 1

zlength = 0.5
dz = 0.05
zsteps = int(zlength/dz)

tlength = 50### s
dt =  0.25*dz**2/(2*mu/rho)
tsteps = int(tlength/dt)

### Velocity, time, and space arrays
initial_vel = np.zeros(zsteps)
time = np.arange(0,tlength, dt)
space = np.arange(0,zlength, dz)


### derivatie function
def derivatives(temp):
    der = np.zeros(len(temp)) ### An empty list with the right length
    der[1:-1] = con*(temp[2:] - 2 * temp[1:-1] + temp[:-2]) / dz**2  ### The derivative array.
    return der


def forward_euler(initial, function ):
    matrix = np.zeros((tsteps,zsteps))  ### An empty array with dimensions of space, and time such that it can be filled with all Temperatures in the space array through all timesteps
    temp1 = np.array(initial.copy()) ### The initial temperatures, copied, the reason that there is an np,array is because otherwise
    ### -the code wont work because it sees temp1 as a function or some such thing, this might be a deepnote error because it doens't seem to happen in spyder.

    matrix[0] = np.array(temp1.copy()) ### The first array in the matrix is the initial temperature array
    
    for i in range(1, tsteps): ### Looping from t = 1 till t = tsteps
        t = time[i] ### making the time accurate, not integer but dt*timestep = time
        temp2 = temp1 + dt*derivatives(temp1) ### Calculate new array with the derviative in accordance with foward euler
        
        temp2[0] = function(t) ###  Uamp*np.sin(w*t) ### Enforcing the boundary conditions
        temp2[-1] = 0
        
        temp1 = np.array(temp2.copy()) ### Updating the array such that the new temperaturearray becomes the old one, and with this array the next temperaturearray can be calculated 
        matrix[i] = np.array(temp2.copy())  ### Saving the temperature array in the matrix for each timestep
    
    return matrix

def RK4(initial, function):
    
    def runge_kutta(old_state, t, dt, derivatives):
        k1 =  derivatives(old_state)
        k2 =  derivatives(old_state + (0.5 * k1))
        k3 =  derivatives(old_state + (0.5 * k2))
        k4 =  derivatives(old_state + k3)
        new_state = old_state + (k1 + 2*k2 + 2*k3 + k4) * dt/6. #multiplying by kappa has to be done because of it being in the heat diffusion equation
        return new_state

    matrix_r = np.zeros((tsteps, zsteps))
    temp1 = np.array(initial.copy())
        
    for i in range(tsteps):
        t = time[i]
        temp2 = runge_kutta(temp1, t, dt, derivatives)
        temp2[0] =  function(t)
        temp2[-1] = 0

        temp1 = np.array(temp2.copy())
        matrix_r[i] = np.array(temp1.copy())

    return matrix_r


def osscilation(t):
    oss = Uamp*np.sin(w*t)
    return oss


def schuif(t):
    if t%T > 0 and t%T < 1*T2 :
        v =1
    else:
        v=-1
    return v


def plotje(matrix, naam, vel): ### Plotting function for a matrix
    plt.figure()
    plt.imshow(matrix, aspect='auto', extent=[0, zsteps, 0, tsteps])
    plt.colorbar(label='x-velocity')
    plt.xlabel('z-axis')
    plt.ylabel('Time')
    plt.title(f'{vel} Velocity plot {naam}.')
    plt.show()


plotje(forward_euler(initial_vel, osscilation), "Forward euler", "Sin sscilating")
plotje(RK4(initial_vel,osscilation), "RK4", "Sin sscilating")


plotje(forward_euler(initial_vel, schuif), "Forward euler", "Con osscilating")
plotje(RK4(initial_vel, schuif), "RK4", "Con osscilating")


