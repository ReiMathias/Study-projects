# -*- coding: utf-8 -*-
"""
Created on Mon Nov 15 21:36:32 2021

@author: mathi
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

def jumper(posisjon, tid): #Function simulating X(t). Takes in list of stae and time, then adds the new jump and the time stayed
    state = posisjon[len(posisjon)-1]
    if state == 0: #special case, no death rate
        return (posisjon + [1], tid + [np.random.exponential(scale=2/10,size = 1)+tid[len(tid)-1]])
    else: #general, when in state i/{0}
        sojurn = np.random.exponential(scale=1/11, size=1) #lamda + mu, this is the sojurn time
        if np.random.uniform(low=0, high=1, size=1) > 5/11: #prob for lambda wining 5/11 = lambda/(lambda+mu). Desiding if a death or a birth occured
            pos_1 = state - 1
            return (posisjon + [pos_1], tid + [sojurn+tid[len(tid)-1]])
        else:
            pos_2 = state + 1
            return (posisjon + [pos_2], tid + [sojurn+tid[len(tid)-1]])

#50 day simulation
(posisjon_2, tid_2) = ([0],[0])
while tid_2[len(tid_2)-1] < 24*50: #running the simulations for 50 days
    (posisjon_2, tid_2) = jumper(posisjon_2, tid_2)
#print(posisjon_2,tid_2)

#12 hour graph
(posisjon_0, tid_0) = ([0,0],[0]) #Two zerows so that the simulation starts in (0,0)
while tid_0[len(tid_0)-1] < 12:
    (posisjon_0, tid_0) = jumper(posisjon_0, tid_0)


fig = plt.figure()
ax = fig.add_subplot(111)
thegraph = ax.step(tid_0+[tid_0[len(tid_0)-1]], posisjon_0)
ax.legend(handles=thegraph, loc=2)
ax.set_title("Urgent care center", fontdict={'fontname': 'Times New Roman', 'fontsize': 21}, y=1.03)
ax.set_xlabel("Time in hours")
ax.set_ylabel("Patients in the UCC") 
ax.set_ylim(0)
ax.set_xlim(0)

#50 days simulation to fin the long-run mean
def longmean(posisjon, tid, sum_ts): #Inputs arre current state, current time and the sum of (sujurn*current state). This simulation is a verry close copy to the function jumper
    if posisjon == 0: #Special case
        sojurn_0 = np.random.exponential(scale=2/10,size = 1) 
        return (1, tid + sojurn_0, sum_ts)
    else:
        sojurn = np.random.exponential(scale=1/11, size=1) 
        if np.random.uniform(low=0, high=1, size=1) > 5/11: 
            pos_1 = posisjon - 1
            return (pos_1, tid + sojurn, sum_ts + (posisjon*sojurn))
        else:
            pos_2 = posisjon + 1
            return (pos_2, tid + sojurn, sum_ts + (posisjon*sojurn))

lamba = 5
for i in range(30): #Here we simulate 30 times to crate a confidence interval
    (posisjon_1, tid_1, mean) = (0,0,0)
    while tid_1 < 50*24:
        (posisjon_1, tid_1, mean) = longmean(posisjon_1, tid_1, mean)
    #print((float(mean)/(float(tid_1))/(lamba)))

#plotting two graphs in (1f)
x = np.arange(0, 1.05, 0.05)
y = 1/(6-5*x)
z = 6/(6-5*x)

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(x, y, color='b', label='Urgent patient')
ax1.plot(x, z, color='g', label='Normal patient')
ax1.set_xlabel("P-value")
ax1.set_ylabel("Hours")
ax1.set_title("Waiting times")
ax1.legend()
ax1.set_ylim(0)
ax1.set_xlim(0)


#Creating the simulation of U(t) and N(t), This can be graphed
def UNjumper(normal, urgent, tid): #Because the exponential function has no memory a normal patients treatment beeing abropted is like p*lambda winning
    state_1 = normal[len(normal)-1]
    state_2 = urgent[len(urgent)-1]
    if state_1 == 0 and state_2 == 0: #Special case where there is no death
        sojurn = np.random.exponential(scale=2/10,size = 1)
        if np.random.uniform(low=0, high=1, size=1) < 0.8: #Probability for there beeing an urgent patient
            return (normal + [0], urgent + [1], tid + [sojurn+tid[len(tid)-1]])
        else:
            return (normal + [1], urgent + [0], tid + [sojurn+tid[len(tid)-1]])
    elif state_2 != 0: #Case when there are urgent patients that has first priority
        sojurn = np.random.exponential(scale=1/11, size=1)
        if np.random.uniform(low=0, high=1, size=1) < 6/11: #there is a death and birth rate mu+p*lambda+(p-1)*lambda = mu+lambda
            return (normal + [state_1], urgent + [state_2-1], tid + [sojurn+tid[len(tid)-1]])
        else:
            if np.random.uniform(low=0, high=1, size=1) < 0.8:
                return (normal + [state_1], urgent + [state_2+1], tid + [sojurn+tid[len(tid)-1]])
            else:
                return (normal + [state_1+1], urgent + [state_2], tid + [sojurn+tid[len(tid)-1]])
    else: #Case when there are no urgent patients only normal patients
        sojurn = np.random.exponential(scale=1/11, size=1) #Code the same as in the elif statment abouve, just here normal patints get served insted
        if np.random.uniform(low=0, high=1, size=1) < 6/11:
            return (normal + [state_1-1], urgent + [state_2], tid + [sojurn+tid[len(tid)-1]])
        else:
            if np.random.uniform(low=0, high=1, size=1) < 0.8:
                return (normal + [state_1], urgent + [state_2+1], tid + [sojurn+tid[len(tid)-1]])
            else:
                return (normal + [state_1+1], urgent + [state_2], tid + [sojurn+tid[len(tid)-1]])

(normal_3, urgent_3, tid_3) = ([0,0],[0,0],[0])
while tid_3[len(tid_3)-1] < 12: #Plotting a joint realisation of N(t) and U(t) in the same cordinat system
    (normal_3, urgent_3, tid_3) = UNjumper(normal_3, urgent_3, tid_3)


fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

graph_n, = ax2.step(tid_3+[tid_3[len(tid_3)-1]], normal_3)
graph_u, = ax2.step(tid_3+[tid_3[len(tid_3)-1]], urgent_3)
graph_n.set_label('Normal patients')
graph_u.set_label('Urgent patients')
ax2.legend()
ax2.set_title("Simulating N(t) and U(t)", fontdict={'fontname': 'Times New Roman', 'fontsize': 21}, y=1.03)
ax2.set_xlabel("Time in hours")
ax2.set_ylabel("Number of patients") 
ax2.set_ylim(0)
ax2.set_xlim(0)


#Approximating expected time in ucc for normal and urgent patient
def UNlongmean(normal, urgent, tid, sum_n, sum_u): #Function only remember current state and the cumulative time. However it stores cumulative (sojurn*state)
    if normal == 0 and urgent == 0: #Function much like the one above. This is the special case where mu = 0
        sojurn = np.random.exponential(scale=2/10,size = 1)
        if np.random.uniform(low=0, high=1, size=1) < 0.8:
            return (0, 1, tid + sojurn, sum_n, sum_u)
        else:
            return (1, 0, tid + sojurn, sum_n, sum_u)
    elif urgent != 0: #If there are urgent patients
        sojurn = np.random.exponential(scale=1/11, size=1)
        ret_sum_n = sum_n + sojurn*normal
        ret_sum_u = sum_u + sojurn*urgent
        if np.random.uniform(low=0, high=1, size=1) < 6/11:
            return (normal, urgent-1, tid + sojurn, ret_sum_n, ret_sum_u)
        else:
            if np.random.uniform(low=0, high=1, size=1) < 0.8:
                return (normal, urgent+1, tid + sojurn, ret_sum_n, ret_sum_u)
            else:
                return (normal+1, urgent, tid + sojurn, ret_sum_n, ret_sum_u)
    else: #If there are none urgent pations the normal will be treated. This code is much like the one in the elif above.
        sojurn = np.random.exponential(scale=1/11, size=1)
        ret_sum_n = sum_n + sojurn*normal
        ret_sum_u = sum_u + sojurn*urgent
        if np.random.uniform(low=0, high=1, size=1) < 6/11:
            return (normal-1, urgent, tid + sojurn, ret_sum_n, ret_sum_u)
        else:
            if np.random.uniform(low=0, high=1, size=1) < 0.8:
                return (normal, urgent+1, tid + sojurn, ret_sum_n, ret_sum_u)
            else:
                return (normal+1, urgent, tid + sojurn, ret_sum_n, ret_sum_u)

for i in range(30):
    (normal_4, urgent_4, tid_4, mean_n_4, mean_u_4) = (0,0,0,0,0) #Simulating 30 simulations to make a confidence intervall
    while tid_4 < 24*50:
        (normal_4, urgent_4, tid_4, mean_n_4, mean_u_4) = UNlongmean(normal_4, urgent_4, tid_4, mean_n_4, mean_u_4)
    #print((float(mean_u_4)/(float(tid_4)))/((5)*(0.8)), (float(mean_n_4)/(float(tid_4)))/((5)*(0.2)))
    #print((float(mean_u_4)/(float(tid_4)))/((5)*(0.8)))
    #print((float(mean_n_4)/(float(tid_4)))/((5)*(0.2)))
    
    
def make_H(t,s):
    n_t = len(t)
    n_s = len(s)
    H = np.zeros((n_t,n_s))
    for i in range(n_t):
        for j in range(n_s):
            H[i][j] = np.abs(t[i]-s[j])
            
    return H


def relation(grid, theta_B, x_B, mu, sigma, phi):
    
    theta_A = np.setdiff1d(grid,theta_B) # Seperating out blocks
    print(theta_A, len(theta_A))
    n_A, n_B = len(theta_A), len(theta_B)
    mu_A = np.array([mu for i in range(n_A)]) # Mean vector.
    mu_B = np.array([mu for i in range(n_B)])
    
    H_A, H_B, H_AB= make_H(theta_A,theta_A), make_H(theta_B,theta_B), make_H(theta_A,theta_B) # Constructing distance matrices
    
    # Making covariance matrices with correlation function as in task \phi = 15
    Sigma_A, Sigma_B, Sigma_AB = sigma*((1+phi*H_A)*np.exp(-phi*H_A)), sigma*((1+phi*H_B)*np.exp(-phi*H_B)), sigma*((1+phi*H_AB)*np.exp(-phi*H_AB))
    
    
    Sigma_Binv= np.linalg.inv(Sigma_B)
    SabSb = np.dot(Sigma_AB,Sigma_Binv)
    
    # Calculating conditional mean and variance
    E_ab = mu_A + np.dot(SabSb,(x_B-mu_B))
    Sigma = Sigma_A - np.dot(SabSb,np.transpose(Sigma_AB))
    
    
    # Making confindence intervals and prediction array
    CIs = []
    pred = []
    for k in range(n_A):
        CIs.append([E_ab[k]-1.64*np.sqrt(Sigma[k][k]),E_ab[k]+1.64*np.sqrt(Sigma[k][k])])
#        print(Sigma[k][k], k)
        pred.append(norm.cdf((0.3-E_ab[k])/(np.sqrt(Sigma[k][k]))))
    
    print('The highest probability is ', max(pred),' at value ', theta_A[np.array(np.where(pred == max(pred)))[0][0] ])
    CIs_ = np.transpose(CIs)
    
    # Plotting
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ax.plot(theta_A,E_ab, label= 'Conditional mean')
    ax.plot(theta_B,x_B,'o', color='red', label= 'Known values')
    ax.fill_between(theta_A, CIs_[0],CIs_[1], facecolor='red', alpha=0.2, label= '90% Confidence interval')
    ax.set_xlabel('Parameter: θ')
    ax.set_ylabel('Fit: y(θ)')
    ax.legend()
    
    fig1 = plt.figure()
    
    ax1 = fig1.add_subplot(111)
    ax1.plot(theta_A,pred)
    ax1.set_ylabel('Probability')
    ax1.set_xlabel('Paramter value: θ')
    return 'Done'

gridspace = np.linspace(0.25,0.5,51)

theta_Ba, x_Ba = np.array([0.3,0.35,0.39, 0.41, 0.45]), np.array([0.5,0.32,0.4,0.35,0.6])

theta_Bc, x_Bc =  np.array([0.3, 0.33, 0.35 ,0.39, 0.41, 0.45]), np.array([0.5, 0.4, 0.32, 0.4, 0.35, 0.6])

##theta_fake, x_fake =  np.array([0.25, 0.3, 0.33, 0.35 ,0.39, 0.41, 0.45]), np.array([0.6, 0.5, 0.4, 0.32, 0.4, 0.35, 0.6]) # Fake 0.25
#
#theta_fake, x_fake =  np.array([0.36, 0.3, 0.33, 0.35 ,0.39, 0.41, 0.45]), np.array([0.38, 0.5, 0.4, 0.32, 0.4, 0.35, 0.6]) # Fake 0.36


relation(gridspace, theta_Ba, x_Ba, mu = 0.5, sigma=0.5**2, phi= 15)

relation(gridspace, theta_Bc, x_Bc, mu = 0.5, sigma = 0.5**2, phi = 15)

