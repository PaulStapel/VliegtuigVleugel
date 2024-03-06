import numpy as np
import matplotlib.pyplot as plt


def circle(center, R, points):
    theta = np.linspace(0, 2*np.pi, points)
    circle = np.zeros(points, dtype = np.complex_)
   
    for i in range(points):
        circle[i] = R* complex(np.cos(theta[i]), np.sin(theta[i])) + center
    return circle



def joukowski(circle):
    jk = np.zeros(len(circle), dtype = np.complex_)
    
    for i in range(len(circle)):

        jk[i] =  circle[i] + 1/circle[i] 


    #plt.xlim(-2.1,2.1)
    #plt.ylim(-2.1,2.1)

    #plt.plot(jk.real,jk.imag)
    
    return jk

# a = (circle(complex(-0.1,0.22), 1.12, 1000))

# # plt.plot(a.real, a.imag)

# joukowski(a)
