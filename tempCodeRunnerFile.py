    
    plt.figure()
    plt.tricontour(wx, wy, wpressure, levels=5000, color='black')
    plt.plot(wing.real, wing.imag, color='black')
    plt.show()