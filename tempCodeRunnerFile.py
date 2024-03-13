    plt.figure()
    plt.tricontour(x, y, pressure, levels=5000, color='black')
    plt.plot(circle.real, circle.imag, color='black')
    plt.show()
    
    plt.figure()
    plt.tricontour(x_trans, y_trans, wpressure, levels=5000, color='black')
    plt.plot(wing.real, wing.imag, color='black')
    plt.show()