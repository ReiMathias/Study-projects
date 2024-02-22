import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.animation import FuncAnimation
from scipy.linalg.lapack import ssbevd
from scipy import stats
import scipy as sp


#WOOO Denne funker

hbar = 1.05E-34  # J⋅s. Redusert planckkonstant
m = 9.11E-31     # kg.  Electronmasse
V0 = 1.5*1.6E-21 # J. Potensialhøyde
dx = 1e-9        # m. 
N = 10            # #. Varierende potensial
N_konstant = 100*N  # #. Område med konstant potensial
x0=0.5*dx*N_konstant #Startposisjon
E0=1.390*V0
k0 = np.sqrt(2.0*m*E0)/hbar  
SIGMA = dx*N_konstant/20.


Ntot = N + 2*N_konstant # Totallengde
print(Ntot)
x = np.linspace(0, dx*Ntot, Ntot) # x-akse
print(x)

#Lager potensialet

def lagV(v0):
    tingtang= [v0*((n-N)/(N*1.0))**2 for n in range(N)]
    V = [v0]*N_konstant
    return np.array(V + tingtang + V)
    



def lagPsi0(sigma,z):
    A = (2*np.pi*sigma**2)**(-0.25)
    Psi_0 = A * np.exp(-(z-x0)**2/(4*SIGMA**2)) * np.exp(1j*k0*z)
    print("Normalization:", dx*np.sum(np.abs(Psi_0)**2))
    return Psi_0


def lagH(Nmax,v):
    diag = hbar**2/(m*dx**2) + v # Diagonal
    sup_diag = np.ones(Nmax)*(-hbar**2/(2*m*dx**2)) # Superdiagonal
    E, psi_n, _ = ssbevd([sup_diag, diag]) # Løser matrisen
    return E, psi_n
    
def lagc(psi_n,psi0):
    psi_n = psi_n.astype(complex)
    c = np.zeros(Ntot, dtype=complex)
    for n in range(Ntot):
        c[n] = np.vdot(psi_n[:,n], psi0)
    return c


def Psi(t, c, psi_n, E):
    list = np.zeros(len(c), dtype= np.complex128)
    for i in range(len(c)):
        list += c[i]*psi_n[:,i]*np.exp(-1j*E[i]*t/hbar) 
    return list
    return np.dot(psi_n,c*np.exp(-1j*E*t/hbar))

def usikkerhet(t):
    funk=SIGMA**2 + ((hbar*t)/(2*m*SIGMA))**2
    return np.sqrt(funk)


def main():
    #plt.rcParams.update({'animation.html':'html5', 'savefig.dpi': 50})
                
        
    dt=1*dx*(m/(2*E0))**.5
    V=lagV(0)
    Psi_0=lagPsi0(SIGMA,x)
    psi_n=lagH(Ntot,V)[1]
    E=lagH(Ntot,V)[0]
    c=lagc(psi_n,Psi_0)
    
    """
    tiden=200
    
    x2mean=np.zeros(tiden, dtype=float)
    xmean2=np.zeros(tiden, dtype=float)
    t=np.zeros(tiden,dtype=float)
    deltaxanal=np.zeros(tiden, dtype=float)
    for i in range(tiden):
        t[i]=i*dt
        prob = np.abs(Psi(i*dt, c, psi_n, E))**2
        x2mean[i]=sp.integrate.simps(prob*(x**2),x)
        xmean2[i]=(sp.integrate.simps(prob*x,x))**2
        deltaxanal[i]=usikkerhet(i*dt)
    deltax = np.sqrt(x2mean-xmean2)
    print(deltax[0])
    
    
    plt.figure('Oppgave 1')
    plt.title('Analytisk mot numerisk ∆x')
    plt.plot(t,deltaxanal, label='Analytisk')
    plt.plot(t,deltax, label='Numerisk')
    plt.xlabel('$t$ (s)',fontsize=16)
    plt.ylabel('$∆x$',fontsize=16)
    plt.legend(loc='best')
    plt.grid()
    plt.show()
    
    """
    fig = plt.figure('Wave packet animation' , figsize=(16,8))
    ymax = 1E7
    ax = plt.axes(xlim=(0, Ntot*dx), ylim=(0, ymax))
    line, = ax.plot([], [], lw=1)
        
    def init():
        line.set_data([], [])
        return line,
        
    
        
    def animate(i):
        t = i*dt
        #Beregner Psi(x,t)
        prob = np.abs(Psi(i*dt, c, psi_n, E))**2
        line.set_data(x, prob)
        return line,
    plt.plot(x,V*ymax)
    plt.xlabel('$x$ (m)',fontsize=20)
        
        
    anim = FuncAnimation(fig, animate, init_func=init, repeat=False,
    frames=200, interval=1, blit=True)
    plt.show()
    
    
    
main()

