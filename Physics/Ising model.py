import numpy as np
import scipy as sp
import itertools as it






def configs(n):

    lst = list(itertools.product([-1, 1], repeat=n))
    
    print(lst)
    
    

configs(5)