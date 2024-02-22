import numpy as np
import math
from matplotlib import pyplot as plt
from matplotlib import animation



#Constants
hbar = 1
m = 1
k0 = 20
w = math.sqrt(k0/m)
xs=250
sigma=80
L=20
Nx=1000


dx=L/(Nx-1)
dt=(dx**2)/5



xvalues=np.linspace(0,20,Nx)
potential=np.zeros(Nx)
psiarray=np.array(np.zeros(Nx),dtype=complex)


def psi0(x,t):
    A = (2*np.pi*sigma**2)**(-0.25)
    verdi = A*(np.exp(-(x-xs)**2/(2*sigma**2)))*np.exp(1j*(k0*x-w*t))
    return verdi
    
    
    
    
    
    
def psi_real(array):
    newarray=np.zeros(Nx, dtype=complex)
    for i in range(0,Nx):
        newarray[i]= 1j
    
    return newarray







    
def main():
    for i in range(0,Nx):
        psiarray[i]= psi0(i,0).real + 1j*psi0(i,-dt/2).imag
        
    plt.plot(xvalues,psiarray.real)
    plt.plot(xvalues,psiarray.imag)
    plt.show()
    
    
    
    

    # First set up the figure, the axis, and the plot element we want to animate
    fig = plt.figure()
    ax = plt.axes(xlim=(0, 20), ylim=(-0.8, 0.8))
    line, = ax.plot([], [], lw=2)
    
    # initialization function: plot the background of each frame
    def init():
        line.set_data([], [])
        return line,
    
    # animation function.  This is called sequentially
    def animate(i):
        
        
        
        line.set_data(xvalues, psiarray.real)
        return line,
    
    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                frames=500, interval=20, blit=True)
    

    
    
main()
    
    
    
    
    
    
    
    