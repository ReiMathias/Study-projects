import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.animation import FuncAnimation

hbar = 1.05E-34
m = 9.11E-31  
V0 = 1.6E-19       
dx = 1E-10
k0 = ((2*m*V0)**.5)/hbar
KB = 400
BB = 3
B =  5
                                #Kjør animasjon med V0 = 1.2*1.6E-19 og E_mid = [2.25*V0] for å se hvordan animajsonen blir når T går mot 1.
                                
#----------------------------VIKTIG--------------------------------------------------    
"""    
For å kjøre animasjon er det 2 ting du først må gjøre:
    -1: Gjør E_mid om til en liste med et element (energien du ønsker å animere for) med å kommentere ut np.linspace
    -2: Endre slik at animasjonskoden som er i bunnen av funksjonen beregnT ikke er kommentert ut
    -3: Kan være greit å kommentere ut plottingen av T som skjer i main ettersom dette plottet ikke vil vise noe når E_mid bare er en enkelt verdi (Dette er ikke noe du må gjøre).
"""
#------------------------------------------------------------------------------------
                                
                                
                                
                                
    
def main():
    tider = np.linspace(0,8E-13,100)
    T_num = []
    T_ana = []
    for i in range(len(tider)):
        T_num.append(beregnT(tider[i]))
    
    plt.figure('Oppgave 2')
    plt.title('Sannsynlighet for å finne elektronet mellom barrierene')
    plt.plot(tider,T_num,'-b',label='$T$')
    plt.xlabel('$t$ (s)')
    plt.ylabel('$P_b(t)$')
    plt.grid()
    plt.show()
    

    
def beregnT(tid):
    E = 3.506952568E-17
    p0 = math.sqrt(2.0*m*E)
    
    V = lagV(V0,KB,BB,B)
    H = lagH(dx,V)
    energy,psi_matrix = np.linalg.eigh(H)
    E_verdier = energy
    psi_matrix_complex = psi_matrix*(1.0 + 0.0j)
    N = len(V)
    x = np.asarray([dx*n for n in range(N)])
    Psi0 = lagPsi0(dx,x,p0)
    #print("Normalization:",dx*np.sum(np.abs(Psi0)**2))
    
    c = np.zeros(N,dtype = np.complex128)
    for n in range(N):
        c[n] = np.vdot(psi_matrix_complex[:,n],Psi0)
    
    Psit = lagPsit(c,psi_matrix_complex,E_verdier,tid)    # tidssteg = 5.9617E-16
    
    #print("v:",p0/m, " \nE: ",E)
    SSF = np.abs(Psit)**2
    
    """
    navn = str(E)
    plt.figure('test' + navn)
    plt.plot(x,.75*V*np.max(np.abs(Psi0)**2)/max(1E-30,np.max(V)),'--')
    plt.plot(x,SSF)
    plt.grid()
    plt.show()
    """
    
    """
    plt.figure('eksempel')
    for i in [0,1,3]:
        plt.plot(x,psi_matrix[:,i],label=r"$\psi_{%.of}(x)$"%(i))
    plt.plot(x,.75*V*np.max(np.abs(Psi0)**2)/max(1E-30,np.max(V)),'--',label='Potensial')
    plt.title('Eigenmodes for the given potential')
    plt.xlabel('$x$ [m]')
    plt.ylabel("$\psi_n(x)$")
    plt.xlim([0,dx*N])
    plt.legend()
    plt.show()
    """
    t = midt(dx,KB+BB-1,KB+BB+B-1,SSF)
    tot = midt(dx,0,N-1,SSF)
    T = t/tot
    
    
    
    """
    fig = plt.figure('Wave packet animation' , figsize=(16,8))
    ymax = 1E8 
    ax = plt.axes(xlim=(0, N*dx), ylim=(0, ymax))
    line, = ax.plot([], [], lw=1)
    
    def init():
        line.set_data([], [])
        return line,
    
    tidssteg = 5E-16
    
    def animate(i):
        t = i*tidssteg
        #Beregner Psi(x,t)
        Psi_t = np.zeros(N,dtype=np.complex128)
        for n in range(N):
            Psi_t = Psi_t + c[n]*psi_matrix_complex[:,n]*np.exp(-1j*E_verdier[n]*t/hbar)
            
        rho_t = np.abs(Psi_t)**2
        line.set_data(x, rho_t)
        return line,
    plt.plot(x,V*ymax/V0)
    plt.xlabel('$x$ (m)',fontsize=20)
    
    
    anim = FuncAnimation(fig, animate, init_func=init, repeat=False,
    frames=int(tid/tidssteg), interval=1, blit=True)
    plt.show()
    """
    
    return T
    

    
    
def lagV(v0,kb,bb,b):
    left = [0]*kb
    bar = [v0]*bb
    midt = [0]*b
    right = [0]*kb
    return np.asarray(left + bar +midt + bar + right)
    
def lagH(dx,V):
    d = [v + hbar**2/(m*dx**2) for v in V]
    e = -hbar**2/(2*m*dx**2)
    N = len(V)
    H = [[0]*N for n in range(N)]
    for i in range(N):
        for j in range(N):
            if i==j:
                H[i][j] = d[i]
            if abs(i-j)==1:
                H[i][j] = e
    return np.asarray(H)
    
def lagPsi0(dx,x,p0):
    sigma = 50*dx
    x0 = 0.5*KB*dx
    normfactor = (2*np.pi*sigma**2)**(-0.25)
    gaussinit = np.exp(-(x-x0)**2/(4*sigma**2))
    planewavefactor = np.exp(1j*x*p0/hbar)
    return normfactor*gaussinit*planewavefactor
    
def lagPsit(c,psi_matrix,E,t):
    list = np.zeros(len(c), dtype= np.complex128)
    for i in range(len(c)):
        list += c[i]*psi_matrix[:,i]*np.exp(-1j*E[i]*t/hbar)
    return list
    
def midt(dx,xa,xb,y):
    return dx*np.sum(np.abs(y[xa:xb])**2)
    
    
        
main()