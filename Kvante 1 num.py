import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.animation import FuncAnimation
from scipy.linalg.lapack import ssbevd
from scipy import stats
import scipy as sp


#Constants
hbar = 1
m = 1
k0 = 20 
xs=5
sigma=1
L=20
Nx=1000
E0=k0**2 /2

dx=L/(Nx-1)
dt=30*(dx**2)/5

T=L/2*k0



xvalues=np.linspace(0,L,Nx)



#Funksjon som lager potensialet
def lagV(v0,bredde):
    
    l=bredde/50
    r = np.zeros(Nx)
    
    width = int(l * Nx/10 )
    r[500-width/2:500+width/2]=v0
    return r


#Funksjon som genererer grunntilstanden til elektronet
def lagPsi0(z):
    A = (np.pi*sigma**2)**(-0.25)
    Psi_0 = A * np.exp(-(z-xs)**2/(2*sigma**2)) * np.exp(1j*k0*z)
    #print("Normalization:", dx*np.sum(np.abs(Psi_0)**2))
    return Psi_0

#Hamilton operatoren
def lagH(v):
    diag = hbar**2/(m*dx**2) + v # Diagonal
    sup_diag = np.ones(Nx)*(-hbar**2/(2*m*dx**2)) # Superdiagonal
    E, psi_n, _ = ssbevd([sup_diag, diag]) # Løser matrisen
    return E, psi_n
    
    
#Utviklingskoeffisientene
def lagc(psi_n,psi0):
    psi_n = psi_n.astype(complex)
    c = np.zeros(Nx, dtype=complex)
    for n in range(Nx):
        c[n] = np.vdot(psi_n[:,n], psi0)
    return c

#Tidsavhengige bølgefunksjonen
def Psi(t, c, psi_n, E):
    list = np.zeros(len(c), dtype= np.complex128)
    for i in range(len(c)):
        list += c[i]*psi_n[:,i]*np.exp(-1j*E[i]*t/hbar) 
    return list
    return np.dot(psi_n,c*np.exp(-1j*E*t/hbar))
    
    
    

    


#Denne funksjonen animerer bølgen, plotter den og regner sannysnlighet for transmisjon og refleksjon

def main():
    V=lagV(E0/2,20)
    #V=lagV(0,0)
    Psi_0=lagPsi0(xvalues)
    psi_n=lagH(V)[1]
    E=lagH(V)[0]
    c=lagc(psi_n,Psi_0)
    
    
    #Animering av bølgepakken
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
    
        y = np.abs(Psi(i*dt, c, psi_n, E))**2
        line.set_data(xvalues, y)
        return line,
    
    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = FuncAnimation(fig, animate, init_func=init,
                                frames=209, interval=20,repeat=False, blit=True)
    
    plt.show()
    """
    #Plotter bølgefunksjonen etter den er sendt gjennom potensialet
    plt.figure(1)
    ax=plt.subplot(111)
    ax.plot(xvalues,Psi(209*dt, c, psi_n, E).real,label='$\Psi_R$')
    ax.plot(xvalues,Psi(209*dt, c, psi_n, E).imag,label ='$\Psi_I$')
    ax.legend()
    
    

    plt.figure(2)
    ax = plt.axes(xlim=(0, 20), ylim=(-1, 1))
    line, = ax.plot(xvalues,np.abs(Psi(209*dt, c, psi_n, E))**2,label='$|\Psi(x,t)|^2$')
    ax.legend()
    #ax = plt.axes(xlim=(0, 20), ylim=(-1, 1.2))
    #line, = ax.plot(xvalues,np.abs(Psi(0*dt, c, psi_n, E))**2,label='$|\Psi(x,t)|^2$')
    ax.legend()
    
    plt.show()
    """
    trans = np.trapz(np.abs(Psi(209*dt, c, psi_n, E)[500:1000])**2,xvalues[500:1000],dx)
    refle = np.trapz(np.abs(Psi(209*dt, c, psi_n, E)[0:500])**2,xvalues[0:500],dx) 
    
    
    
    print(np.trapz(np.abs(Psi(0*dt, c, psi_n, E))**2,xvalues,dx))
    print(trans,refle)
    
    

#Regner ut sannsynligheten for transmisjon og refleksjon for 50 forskjellige styrker og plotter det 
def main2():
    
    transmisjon=[]
    refleksjon=[]
    energies=[]
    
    for j in range(0,50):
        V=lagV(0 + (3*E0/2)*j/50,20)
        Psi_0=lagPsi0(xvalues)
        psi_n=lagH(V)[1]
        E=lagH(V)[0]
        c=lagc(psi_n,Psi_0)
        
        energies.append(  E0/ (3*E0/2)*j/50)
        
        transmisjon.append( np.trapz(np.abs(Psi(209*dt, c, psi_n, E)[500:1000])**2,xvalues[500:1000],dx))
        refleksjon.append( np.trapz(np.abs(Psi(209*dt, c, psi_n, E)[0:500])**2,xvalues[0:500],dx))
        
    
    
    transmisjon = np.array(transmisjon)
    refleksjon = np.array(refleksjon)
    energies=np.array(energies)
    
    plt.figure(1)
    ax=plt.subplot(111)
    ax.plot(energies,transmisjon,label='P(Transmission)')
    ax.plot(energies,refleksjon, label='P(Reflcetion)')
    ax.legend()
    plt.show()
    
    
#Regner ut sannsynligheten for transmisjon og refleksjon for 50 forskjellige bredder og plotter det 
def main3():
    transmisjon=[]
    refleksjon=[]
    bwidth=[]
    
    for j in range(0,50):
        V=lagV(9*E0/10,j)
        Psi_0=lagPsi0(xvalues)
        psi_n=lagH(V)[1]
        E=lagH(V)[0]
        c=lagc(psi_n,Psi_0)
        
        bwidth.append(j/50 )
        
        transmisjon.append( np.trapz(np.abs(Psi(209*dt, c, psi_n, E)[500:1000])**2,xvalues[500:1000],dx))
        refleksjon.append( np.trapz(np.abs(Psi(209*dt, c, psi_n, E)[0:500])**2,xvalues[0:500],dx))
        
    
    
    transmisjon = np.array(transmisjon)
    refleksjon = np.array(refleksjon)
    bwidth=np.array(bwidth)
    
    ax=plt.subplot(111)
    ax.plot(bwidth,transmisjon,label='P(Transmission)')
    ax.plot(bwidth,refleksjon,label='P(Reflection)')
    ax.legend()
    plt.show()
    
    
main()
    









