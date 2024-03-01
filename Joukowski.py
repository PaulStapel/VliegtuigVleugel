import numpy as np
import matplotlib.pyplot as plt


def joukowski(center, R ,points):
    theta = np.linspace(0, 2*np.pi, points)
    x , y = np.zeros(points),  np.zeros(points)
    
    for i in range(points):
        circle = R* complex(np.cos(theta[i]), np.sin(theta[i])) + center
        jk =  1/circle + circle
        x[i] = jk.real
        y[i] = jk.imag

    plt.xlim(-2.1,2.1)
    plt.ylim(-2.1,2.1)
    plt.plot(x,y)
    
joukowski(complex(-0.1, 0.22), 1.12, 1000)




    
    
