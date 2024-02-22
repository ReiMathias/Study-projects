# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 10:58:44 2021

@author: mathi
"""

import numpy as np
import matplotlib.pyplot as plt


L = 10

pE = np.sqrt(8*0.5)/3
sig = 1

xx = np.linspace(0,L, 1000)

def make_H(t,s):
    n_t = len(t)
    n_s = len(s)
    H = np.zeros((n_t,n_s))
    for i in range(n_t):
        for j in range(n_s):
            H[i][j] = np.abs(t[i]-s[j])
            
    return H

# Exponetial type correlation
mu = np.zeros(len(xx))
Sig = sig**2*np.exp(-np.abs(make_H(xx,xx))*pE)

# Simulation
nR = 3
yMat = np.zeros((nR,len(xx)))
Lt = np.transpose(np.linalg.cholesky(Sig))
sample = np.random.normal(0,sig,len(xx))

for i in range(nR):
    yMat[i] = np.dot(Lt,sample)
yAll = yMat[1]

for i in range(nR):
    plt.plot(xx,yMat[i])
