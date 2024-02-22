import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt




#Constants
d = 3 
if d==1:
    Omega = 2
elif d==2:
    Omega= 2*np.pi
elif d==3:
    Omega= 4*np.pi



def integrand(r,betaU,R):
    return r**(d-1)  *  (1 - np.exp(betaU*(R/r)**alpha))
    
    
def B_2(betaU,J):
    return Omega/2 *( (R**d /d + J[0]))
    



for j in range(5,8):
    alpha=j
    betaU = np.linspace(2,0.1,100)
    R = 1
    B2=np.zeros(len(betaU))
    
    
    for i in range(len(betaU)):
        I = integrate.quad(integrand, R, np.inf, args=(betaU[i], R))
        B2[i]=B_2(betaU[i],I)
        
    
    plt.plot(1/betaU, 2*d*B2/(Omega*R**d), label=r'$\alpha$ ='+ str(j) )
    plt.grid(True)
    plt.xlabel(r'$\frac{1}{\beta U_0}$')
    plt.ylabel(r'$\frac{2dB_2}{\Omega R^{d}}$')
    plt.legend()
    plt.title('d=3')
    plt.show()

    
    
    
