import numpy as np
import math
from matplotlib import pyplot as plt
from matplotlib import animation



Nx=1000
L=20
psiarray= np.zeros(Nx, dtype=complex)
xs=5
hbar = 1
m = 1
k0 = 20
sigma=1
L=20
vg=k0
w=k0*vg

dx=L/(Nx-1)
dt=(dx**2)/5

def potensial():
    V = np.zeros(Nx)
    
    return np.array(V)
    

def psi0(x,t):
    A = (2*np.pi*sigma**2)**(-0.25)
    #A = (1/(np.exp(((L-xs)**2)/sigma**2))*( sigma**2 /-2*L+xs  ) -sigma**2/xs)
    verdi = A*(np.exp(-(x-xs)**2/(2*sigma**2)))*np.exp(1j*(k0*x-w*t))
    return verdi
    
    

xvalues=np.linspace(0,L,Nx)
psiarray=np.zeros(Nx, dtype=complex)
for i in range(Nx):
    if i>=1 and i<=Nx-2:
        psiarray[i]= psi0(xvalues[i],0).real + 1j*psi0(xvalues[i],dt/2).imag

V=potensial()
#print(V)
#print(psiarray[30].real)
#print(xvalues)



def main():
    plt.figure(1)
    plt.plot(xvalues,psiarray.real)
    plt.show()
    
    
    for j in range(100):
        for i in range(1,Nx-1):
            psiarray[i] = psiarray[i].real - (dt/2*dx**2) * (  psiarray[i+1].imag + psiarray[i-1].imag - 2*psiarray[i].imag) + dt*V[i]*psiarray[i].imag + 1j*psiarray[i].imag
    
    plt.figure(2)
    plt.plot(xvalues,psiarray.real)
    plt.show()
    
    
    
main()
