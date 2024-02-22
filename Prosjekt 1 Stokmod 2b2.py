import numpy as np
import matplotlib.pyplot as plt

R = 12 #This is how manny times a day whe update the graph
N = 59+1/R #Notice how int(R*N) is used, N*R = 12*59+1 
graphs = 10
lamda = [1.5 for i in range(graphs)] #Simple list with 10 elements all equal to 1.5
X_T = [np.random.poisson(i*(1/R), size=int(N*R)) for i in lamda] #This is a 2D list (C and S are also 2D lists). it consists of 10 coppies of numbers drawn from a poisson distrobution such that there are N*R entries. 
#The larger you set "R" the less likely it is that the entries is grater than 1 in X_T. 
C = [[np.sum(np.random.exponential(scale=1/10, size = X_T[j][i])) for i in range(len(X_T[j]))] for j in range(graphs)] #This list coppies X_T and replaces its 0,1,2,... (that came from the poisson) with a number drawn from an exponential function.
#The reason np.sum is used is because it is possible that the poisson distrobution gives us a nuber larger than 1, such that you reseave multiple values form the exponential distrobution.
S = [[np.sum(X[0:i]) for i in range(int(N*R))] for X in C] #S is the list that accumulates all the values given by the list C such that our graph Z(t) will represent the sum of exponential functions where they accumulate as the poisson distrobution jumps. 

X = np.linspace(0, N, int(N*R)) #We plot the graph here
graphs = [plt.step(X, S[i])[0] for i in range(graphs)]
plt.legend(handles=graphs, loc=2)
plt.title("Z(t)", fontdict={'fontname': 'Times New Roman', 'fontsize': 21}, y=1.03)
plt.xlabel("Days")
plt.ylabel("Dollar in million(s)") 
plt.ylim(0)
plt.xlim(0)
plt.show()