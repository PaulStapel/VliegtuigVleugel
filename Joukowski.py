import numpy as np
import matplotlib.pyplot as plt


def circle(center, R, points):
    theta = np.linspace(0, 2*np.pi, points) #define angles
    circle = R* (np.cos(theta) + 1j*np.sin(theta))+ center #calculate points on circle
    return circle

def joukowski(circle):
    jk = circle+ 1/circle # calculate joukowski transform 

    #plt.xlim(-2.1,2.1)
    #plt.ylim(-2.1,2.1)
    #plt.plot(jk.real,jk.imag)
    
    return jk

# circle = (circle(center, R, points))

### Plot for the circle

# plt.plot(circle.real, circle.imag)
# plt.xlim(-1.5,1.5)  
# plt.ylim(-1.5,1.5)
# plt.title(f'Circle in {points} points.')
# plt.xlabel("Length in m.")
# plt.ylabel("Length in m.")
# plt.show()

#Test functies
# a = (circle(complex(-0.1,0.22), 1.12, 1000))
# # plt.plot(a.real, a.imag)
# joukowski(a)
