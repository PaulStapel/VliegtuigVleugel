import numpy as np
import matplotlib.pyplot as plt


def joukowski(center, R ,points):
    theta = np.linspace(0, 2*np.pi, points)
    
    x = np.zeros(points)
    y = np.zeros(points)
    
    for i in range(points):
        circle = R* complex(np.cos(theta[i]), np.sin(theta[i])) + center
        jk =  1/circle + circle
        x[i] = jk.real
        y[i] = jk.imag


    plt.plot(x,y)
    
joukowski(complex(0.5,0.5), 1, 1000)
    
    
    
    
    