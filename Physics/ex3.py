import itertools
from matplotlib import pyplot as plt
import numpy as np



beta_list = np.linspace(0.005,1,100)
J_par = 1
J_perp = 1
B = 1
B_values = np.linspace(0,1,100)



def spin_conf(n):
     return( list(itertools.product([1,-1], repeat = n)))


 
def P_l_m(l,m,spin_config,beta,betaB):
    A,Bt,C = 0,0,0
    sigma_i = spin_config[l]
    sigma_k = spin_config[m]
    A = beta*J_par*np.matmul(sigma_i,sigma_k)
    Bt = beta*J_perp*np.matmul(sigma_i,np.roll(sigma_i,1))
    C = betaB/2*np.sum(sigma_i+sigma_k)
    
    return(np.exp(A+Bt+C))


def P(n,beta,betaB):
    matrix = np.zeros([2**n,2**n])
    spin_config = spin_conf(n)
    for i in range (2**n):
        for j in range (2**n):
            matrix[i][j] = P_l_m(i,j,spin_config,beta,betaB)
    return matrix

    
def get_eig(mat):
    e = np.linalg.eig(mat)[0]
    return(np.real(e))
  
    
def plot_eig(n):
    for beta in beta_list:
        P = P(n,beta,B)
        eigen_val = get_eig(P)
        b = np.ones([eigen_val.shape[0]])*beta
        #print(b)
            
        plt.semilogy(b,eigen_val,'bx')
    plt.title("eig for n = "+str(n) + r" for $J_{perp} = J_{par} = 1, B = 1$")
    plt.ylabel(r"All eig for given $\beta$")
    plt.xlabel(r"$\beta$")
    plt.show()
    


   
def get_max_eig(n,betaB):
    max_eig = np.zeros([100])
    for i in range(100):
        Matrix = P(n,beta_list[i],betaB[i])
        eigen_val = get_eig(Matrix)
        max_eig[i] = np.amax(eigen_val)
    return max_eig



def diff(n,Bval):
    delta_betaB = (beta_list[3]-beta_list[2])*(Bval[3]-Bval[2])
    print(delta_betaB)
    eigvalpluss = np.log( get_max_eig(n,Bval*beta_list + delta_betaB))
    eigval = np.log( get_max_eig(n,Bval*beta_list))
    return((eigvalpluss - eigval)/(delta_betaB))
    

def double_diff(n):
    delta_beta = beta_list[3]-beta_list[2]
    eigvalpp = np.log(get_max_eig(n,beta_list + 2*delta_beta))
    eigvalp = np.log(get_max_eig(n,beta_list + delta_beta))
    eigval = np.log(get_max_eig(n,beta_list))
    return( (eigvalpp - 2*eigvalp + eigval)/(delta_beta**2))



def plot_mag(n,Bval):
    deriv = diff(n,Bval)
    plt.plot(Bval*beta_list, deriv/n)
    plt.title(r"Magnetisation  for $J_{perp} = J_{par} = 1$.")
    plt.xlabel(r"$\beta\cdot B$")
    plt.ylabel(r"$m$")
    plt.show()
    



for i in range(1,6):
    plt.plot(beta_list,beta_list**2*double_diff(i)/i)
plt.xlabel(r"$\frac{1}{kT}$")
plt.ylabel(r"$\frac{C_B(T)}{Nnk}$")
plt.title(r"Specific heat for $J_{perp} = J_{par} = 1$, B = 1.")
plt.show()

for conf in range(0,2**5):
    print(spin_conf(5)[conf])

"""

for j in range(1,6):
    plot_derivatives(j,B_values)



"""



