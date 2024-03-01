import numpy as np
import matplotlib.pyplot as plt




x0, y0 = complex(0.1,0.4)
def joukowski(center, R ,points):
    theta = np.linspace(0, 2*np.pi, points)
    x , y = np.zeros(points),  np.zeros(points)
    
    for i in range(points):
        circle = R* complex(np.cos(theta[i]), np.sin(theta[i])) + center
        jk =  circle #+1/circle 
        x[i] = jk.real
        y[i] = jk.imag

    plt.xlim(-2.1,2.1)
    plt.ylim(-2.1,2.1)
    plt.plot(x,y)
    return np.array([x,y])
    
b = joukowski(complex(0.1, 0.40), 1, 1000)

j = 500
k = 500 
matrix = np.zeros((j,k,2))

def psi(x,y):
    
    R = np.sqrt((x - x0)**2 + (y-y0)**2)
    Theta = np.arccos((x - x0)/R)
    
    psi = (R-(Theta/R)**2)*np.sin(Theta)
    return np.array(psi)
    
    