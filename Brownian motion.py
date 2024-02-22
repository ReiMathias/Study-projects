# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 12:09:26 2021

@author: mathi
"""

import numpy as np
import matplotlib.pyplot as plt



xx = np.linspace(0,100,10001)
mu, sigma  = 0, 1


def brown():
    z = np.random.normal(mu, sigma, len(xx))
    
    y = np.zeros(len(xx))
    for i in range(len(xx)-1):
        y[i+1] = y[i] + np.sqrt(xx[i+1]-xx[i])*sigma*z[i]
        
    plt.plot(xx,y)



#brown()
j=0
while j<10:
    brown()
    j+=1