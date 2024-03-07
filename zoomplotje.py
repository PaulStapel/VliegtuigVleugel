import matplotlib.pyplot as plt


def zoomplotje(joukowski_wing, x,y,z):
    xmax_wing,xmin_wing = max(joukowski_wing.real), min(joukowski_wing.real)
    ymax_wing,ymin_wing = max(joukowski_wing.imag), min(joukowski_wing.imag)
    
    xl,yl = (xmax_wing - xmin_wing), (ymax_wing - ymin_wing)

    
    ypoint = joukowski_wing.imag[joukowski_wing.real.index(xmax_wing)]
    
    xmax,xmin = xmax_wing + 0.1*xl, xmax_wing - 0.1*xl
    ymax,ymin = ypoint + 0.1*yl,  ypoint - 0.1*yl
    
    plt.plot()
    plt.tricontour(x,y,z)
    plt.axis([xmin,xmax,ymin,ymax])
    plt.show()
    
    