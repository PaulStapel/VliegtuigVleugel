import numpy as np
import matplotlib.pyplot as plt

class VectorField2D():   
    def __init__(self, xfunc, yfunc, numx=None, numy=None, xmin=None, ymin=None, xmax=None, ymax=None):
        self.xfunc = xfunc
        self.yfunc = yfunc
        self.numx = numx if numx is not None else 5 #5 gridpoints to start on
        self.numy = numy if numy is not None else 5 #5 gridpoints to start on
        self.xmin = xmin if xmin is not None else -5
        self.ymin = ymin if ymin is not None else -5
        self.xmax = xmax if xmax is not None else 11
        self.ymax = ymax if ymax is not None else 11

        ## Initialise gridspace
        self.xlin = np.linspace(self.xmin, self.xmax, self.numx)
        self.ylin = np.linspace(self.ymin, self.ymax, self.numy)
        self.xx, self.yy = np.meshgrid(self.xlin, self.ylin)

    def __repr__(self) -> str:
        return f"{type(self).__name__}"
    
    def get_velocity(self, t, x, y) -> tuple:
        return self.xfunc(t,x,y), self.yfunc(t,x,y)
    
    def trajectory(self, x0, y0, tmin, tmax, nt) -> tuple:
        # this function computes the trajectory of a particle for starting point (x0,y0) at time tmin,
        # up to time tmax, and with nt time steps.  
        
        # define some vectors to store results
        xtraj=np.zeros(nt+1)  #indices of this vector are 0,1,... nt-1,nt
        ytraj=np.zeros(nt+1)
        ttraj = np.linspace(tmin,tmax,nt+1)
        dt=(tmax-tmin)/nt
        
        # fix the initial conditions 
        xtraj[0]=x0
        ytraj[0]=y0
            
        # now loop over all time steps 
        for it in range (0,nt):    #this in fact means: 0,1,2,... nt-1
            (u,v) = self.get_velocity(ttraj[it],xtraj[it],ytraj[it]) 
        
            xtraj[it+1]=xtraj[it]+dt*u
            ytraj[it+1]=ytraj[it]+dt*v
            
        return ttraj, xtraj,ytraj

    def streakline(self, x0, y0, tmin, tmax, nt) -> tuple:
        # this function computes the streakline for starting point (x0,y0) at time tmin,
        # with nt time steps. Velocity field will be generated using the functions from above. 
        
        # define some vectors to store results
        xstreak =np.zeros(nt+1)  #indices of this vector are 0,1,... nt-1,nt
        ystreak =np.zeros(nt+1)
        time = np.linspace(tmin,tmax, nt+1)
        dt=(tmax-tmin)/nt

        for index in range(0,nt+1):
            tstart = tmin + index*dt
            t, xes, yes = self.trajectory(x0, y0, tstart, tmax, nt)
            xstreak[index] = xes[-1]
            ystreak[index] = yes[-1]
            
        return time, xstreak, ystreak 

    def streamline(self, x0, y0, t, L, nl) -> tuple:

        # this function computes the streamline for starting point (x0,y0) at time tmin,
        # ,with nt time steps. Velocity field will be generated using the functions from above. 
        
        # define some vectors to store results
        xstream=np.zeros((nl+1,nl+1))  #indices of this vector are 0,1,... nl-1,nl
        ystream=np.zeros((nl+1,nl+1))
        dl=L/nl

        xstream[0] = x0
        ystream[0] = y0
        
        # now loop over all time steps 
        for it in range (0,nl+1):    #this in fact means: 0,1,2,... nl-1
            
            (u,v) = self.get_velocity(t, xstream[it], ystream[it])

            dx = u*dl/np.sqrt(u**2 + v**2)
            dy = v*dl/np.sqrt(u**2 + v**2)

            xstream[it+1] = xstream[it]+dx
            ystream[it+1] = ystream[it]+dy

        return xstream, ystream

    def gradient(self, t, x, y, h) -> tuple:
        """Take the gradient of a coordinate for one dimension using the central derivative"""
        gradx = (self.xfunc(t, x+h,y) - self.xfunc(t, x-h,y))/(2*h)
        grady = (self.yfunc(t, x, y+h) - self.yfunc(t, x, y-h))/(2*h)
        return gradx, grady
    
    def curl(self, t, x, y, h=0.01) -> float:
        dydx = (self.yfunc(t, x+h, y) - self.yfunc(t, x-h,y))/(2*h)
        dxdy = (self.xfunc(t, x, y+h) - self.xfunc(t, x, y-h))/(2*h)
        return (dydx - dxdy)
    
    def divergence(self, t, x, y, h=0.01) -> float:
        gradx, grady = self.gradient(t, x,y,h)
        return gradx+grady

    def show_streamlines(self, t, L, nl) -> None:
        xstreams, ystreams = self.streamline(self.xx, self.yy, t, L, nl)
        for x,y in zip(xstreams, ystreams):
            plt.plot(x,y)
            plt.xlim(self.xmin, self.xmax)
            plt.ylim(self.ymin, self.ymax)
        plt.show()

    def show_streakline(self, x0, y0, tmin, tmax, nt) -> None:
        time, xstreak, ystreak = self.streakline(x0,y0,tmin,tmax,nt)
        plt.plot(xstreak,ystreak, ctt = time)
        plt.xlim(min(xstreak)*0.9, max(xstreak)*1.1)
        plt.ylim(min(ystreak)*0.9, max(ystreak)*1.1)
        plt.show()
    
    def show_trajectory(self, x0, y0, tmin, tmax, nt) -> None: 
        time, xtraj, ytraj = self.trajectory(x0,y0,tmin,tmax,nt)
        plt.plot(xtraj,ytraj, ctt = time)
        plt.xlim(min(xtraj)*0.9, max(xtraj)*1.1)
        plt.ylim(min(ytraj)*0.9, max(ytraj)*1.1)
        plt.show()

    def fun():
        print('Yay')
        
