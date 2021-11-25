# Used to generate figures x, y, and z

import scipy.special #Used for Bessel I 
import numpy as np #Used for arrays and matrix like structures
import matplotlib.pyplot as plt #Used for plotting 
from datetime import datetime #Used in runtime calculation

def v_analytic(r,z,gamma, terms):

    
    a= r*(1-z/gamma) 
    sums = 0
    n = 0
    for n in range(1,int(terms)):
        sums += ((scipy.special.iv(1,(n*np.pi*r)/gamma))/(n*scipy.special.iv(1,(n*np.pi)/gamma)))*np.sin(n*np.pi*z/gamma)    
    
    return a - (2/np.pi)*sums


def l_2(v, nr , nz):
    
    v_norm_l = np.sqrt(1/((nr+1)*(nz+1)) * np.sum(np.power(v,2)))
    return v_norm_l

terms =list(np.linspace(2,50,49))

v0_5 = np.zeros((100,50))

v1 = np.zeros((100,100))

v2_5 = np.zeros((100,250))

l2_0_5 = []

l2_1 = []

l2_2_5 = []

loop_start = datetime.now()
for term in terms:
    for r in range(100): 
        for z in range(50):
            
            v0_5[r][z] = v_analytic(1/(101)*r, 1/(50+1)*z*0.5, 0.5, term)
            
        for j in range(100): 
            
            v1[r][z] = v_analytic(1/(101)*r, 1/(100+1)*z*1, 1, term)
            
        for k in range(250):
            
            v2_5[r][z] = v_analytic(1/(101)*r, 1/(250+1)*z*2.5, 2.5, term)
        
    l2_0_5.append(l_2(v0_5, 100, 50))
    l2_1.append(l_2(v1, 100, 100))
    l2_2_5.append(l_2(v2_5, 100, 250))
    
print(datetime.now() - loop_start)  

plot1 = plt.figure(1)
plt.title(r'l2 norm for $\Gamma$ = 0.5')
plt.plot(terms, l2_0_5, label = r"$\Gamma$ = 0.5")

plot1 = plt.figure(2)
plt.title(r'l2 norm for $\Gamma$ = 1')
plt.plot(terms, l2_1, label = r"$\Gamma$ = 1")

plot1 = plt.figure(3)
plt.title(r'l2 norm for $\Gamma$ = 2.5')
plt.plot(terms, l2_2_5, label = r"$\Gamma$ = 2.5")
plt.legend()
plt.show()
