import numpy as np
import matplotlib.pyplot as plt

def kron(x): #A Function that returns 1 if the input is larger than 8
    if float(x) > 8:
        return 1
    else:
        return 0

sant = 0
for i in range(1000): #1000 Simulations
    sum = 0
    for i in range(int(np.random.poisson(lam=1.5*59, size=1))): #Draw a number "L" from a poisson distrobution where t = 59 and lamda =1.5
        sum += float(np.random.exponential(scale=1/10, size=1)) #Sum "L" times a exponential function whit parameter lamda = 10
    sant += kron(sum) #Run the sum in the function kron() to see if it is over 8
print(sant/1000)