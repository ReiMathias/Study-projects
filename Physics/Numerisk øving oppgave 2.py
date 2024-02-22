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
V0 = 1.6E-19 # J. Potensialhøyde
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


#Lager potensialet

def lagV(v0):
    tingtang= [v0*((n-N)/(N*1.0))**2 for n in range(N)]
    V = [v0]*N_konstant
    return np.array(V + tingtang + V)
    


#Funksjon som lager Psi0
def lagPsi0(sigma,z):
    A = (2*np.pi*sigma**2)**(-0.25)
    Psi_0 = A * np.exp(-(z-x0)**2/(4*sigma**2)) * np.exp(1j*k0*z)
    print("Normalization:", dx*np.sum(np.abs(Psi_0)**2))
    return Psi_0

#Numerisk hamilton operator
def lagH(Nmax,v):
    diag = hbar**2/(m*dx**2) + v # Diagonal
    sup_diag = np.ones(Nmax)*(-hbar**2/(2*m*dx**2)) # Superdiagonal
    E, psi_n, _ = ssbevd([sup_diag, diag]) # Løser matrisen
    return E, psi_n
    
    
#Utvilkingskoeffisientene for å lage tidsavhengig bølgefunksjon
def lagc(psi_n,psi0):
    psi_n = psi_n.astype(complex)
    c = np.zeros(Ntot, dtype=complex)
    for n in range(Ntot):
        c[n] = np.vdot(psi_n[:,n], psi0)
    return c

#Funksjoon for å lage tidsavhengigbølgefunkjson
def Psi(t, c, psi_n, E):
    list = np.zeros(len(c), dtype= np.complex128)
    for i in range(len(c)):
        list += c[i]*psi_n[:,i]*np.exp(-1j*E[i]*t/hbar)
    return list
    return np.dot(psi_n,c*np.exp(-1j*E*t/hbar))

#Funksjon som deriverer numerisk
def derivera(y):
    f=np.zeros(len(y)+2,dtype=complex)
    for i in range(1,len(y)-1):
        f[i]=(y[i+1]-y[i-1])/2*dx
    return f[0:len(y)]
    


def main():
                        
    dt=3*dx*(m/(2*E0))**.5
    V=lagV(V0)
    Psi_0=lagPsi0(SIGMA,x)
    psi_n=lagH(Ntot,V)[1]
    E=lagH(Ntot,V)[0]
    c=lagc(psi_n,Psi_0)
    
    print(x0)
    tiden=200
    
    t=np.zeros(tiden,dtype=float)
    x2mean=np.zeros(tiden, dtype=float)
    xmean2=np.zeros(tiden, dtype=float)
    p2mean=np.zeros(tiden, dtype=float)
    pmean2=np.zeros(tiden, dtype=float)
    
    
    
    for i in range(tiden):
        t[i]=i*dt
        prob = np.abs(Psi(i*dt, c, psi_n, E))**2
        cpsi=np.conjugate(Psi(i*dt, c, psi_n, E))
        dpsi=derivera(Psi(i*dt, c, psi_n, E))
        d2psi=derivera(dpsi)
        
        x2mean[i]=sp.integrate.simps(prob*(x**2),x)
        xmean2[i]=(sp.integrate.simps(prob*x,x))**2
        p2mean[i]=(sp.integrate.simps(np.abs(cpsi*d2psi),x,dx))
        pmean2[i]=(sp.integrate.simps(np.abs(cpsi*dpsi),x,dx))**2
    deltax=np.sqrt(x2mean-xmean2)
    deltap=np.sqrt(p2mean-pmean2)
    
    prod=np.zeros(tiden, dtype=float)
    #hlist=np.zeros(tiden, dtype=float)
    for q in range(tiden):
        prod[q]=deltax[q]*deltap[q]
        #hlist[q]=hbar/2
    
    plt.figure('Oppgave 2.1')
    plt.title('Usikkerhetsprodukt')
    plt.plot(t,prod)
    plt.xlabel('$t$ (s)',fontsize=16)
    plt.ylabel('$Δx(t)Δp(t)$',fontsize=16)
    plt.grid()
    plt.figure('Oppgave 2.2')
    plt.title('Usikkerhet i posisjon')
    plt.plot(t, deltax)
    plt.xlabel('$t$ (s)',fontsize=16)
    plt.ylabel('$Δx(t)$',fontsize=16)
    plt.grid()
    plt.figure('Oppgave 2.3')
    plt.title('Usikkerhet i impuls')
    plt.plot(t,deltap)
    plt.xlabel('$t$ (s)',fontsize=16)
    plt.ylabel('$Δp(t)$',fontsize=16)
    plt.grid()
    plt.show()
    
    
    
main()