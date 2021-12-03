# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 16:42:06 2021

@author: Nick Navas
@class: MAT 462

PDE WE ARE SOLVING

v_t = (1/Re)[v_rr +(1/r)v_r - v/r^2 +v_zz]
TO-DO
------
Things to fix:
    
Things to add:
    Non-Linear System
"""
# IMPORTS #
import scipy.special  
import numpy as np #Used for arrays and matrix like structures
import matplotlib.pyplot as plt #Used for plotting 
from mpl_toolkits.mplot3d import Axes3D #Used for some 3D-plotting 
from datetime import datetime #Used in runtime calculation
import matplotlib.animation as animation

import IMPTRI #Part II improved

#%matplotlib qt #RUN THIS IN CONSOLE BEFORE RUNNING CODE#

# CONSTANTS #
start = datetime.now()
Re = 1 #Reynolds number 
dt = 1e-5 #Time step
Size = 50 #Size of matrix
dr = 1/Size #delta-r
dz = 1/Size #delta-z
w = IMPTRI.v_exact
Gamma = 1.5 #Aspect ratio
tolerance= 1e-5 #Min difference to test for steady state
r = np.linspace(0,1,Size)
max_interation = 7500 #sort of a manual killswitch to stop code from 
inter_num = 0 #interation counter (used in checking loop)

# FUNCTIONS #
def plot_contour(Phi, gamma, title, filename=None, cmap=plt.cm.turbo):
    """
    Credit to Dr. Oliver Beckstein here. I took this from the resources of my  
    PHY 202 class on solving Laplace's equation numerically but made the title
    modification.
    
    Plot Phi as a contour plot.
    
    Arguments
    ---------
    Phi : 2D array
          potential on lattice
    gamma : gamma value
    title : String type. Creates the title of out contour plot. 
    filename : string or None, optional (default: None)
          If `None` then show the figure and return the axes object.
          If a string is given (like "contour.png") it will only plot 
          to the filename and close the figure but return the filename.
    cmap : colormap
          pick one from matplotlib.cm          
    """
    fig = plt.figure(figsize=(5,4))
    ax = fig.add_subplot(111)

    x = np.arange(Phi.shape[0])
    y = np.arange(Phi.shape[1])
    X, Y = np.meshgrid(x, y)
    Z = Phi[X, Y]
    cset = ax.contourf(X, Y, Z, 20, cmap=cmap)
    
    ax.set_title(title)
    ax.set_xlabel('r')
    ax.set_ylabel(r'10$\times$z/$\Gamma$')
    ax.set_aspect(1)
    
    cb = fig.colorbar(cset, shrink=0.5, aspect=5)
    
def l_2(v, nr = Size, nz = Size):
    v_norm_l = np.sqrt(1/((nr+1)*(nz+1)) * np.sum(np.power(v,2)))
    return v_norm_l

def RHS(v):
    u = np.zeros((Size,75))
    for i in range(1,Size-1):
        for j in range(1,75-1):
            a = (v[i+1,j] - 2*v[i,j] + v[i-1,j])
            b = (v[i+1,j] - v[i-1,j]) / (2*i)
            c = v[i,j]/i**2
            d = (v[i,j+1] - 2*v[i,j] +v[i,j-1])
            #print(f'a: {a}, b: {b} , c: {c} , d: {d} ')
            u[i,j] = (1/Re) * (a+b-c+d)/dr**2
    #u[:,0] = r
    return(u)    

def Heun(v):
    predict = v + RHS(v) * dt
    
    corrected = v + 0.5 * dt * (RHS(v) + RHS(predict))
    return corrected
    
def check(v0, v1, tol = tolerance):
    if ((l_2(v0) - l_2(v1))/l_2(v0)) < tol:
        
        #print((l_2(v1) - l_2(v0))/l_2(v0))
        return True
        
    else:
        return False

def plotheatmap(v_k, k):
    # Clear the current plot figure
    plt.clf()
    plt.title(f"Time Evolution at time step {k}")
    plt.xlabel("r")
    plt.ylabel("z")

    # This is to plot v_k (v at time-step k)
    plt.contourf(v_k, cmap=plt.cm.turbo, vmin=0, vmax=1)
    plt.colorbar()

    return plt


v = np.zeros((Size,int(Size*Gamma)))
v[:,0] = r 

#matrix for animation
animation_matrix = np.zeros((1,50,int(50*Gamma))) # This is where I store the timesteps for the animation
animation_matrix[0,:,0] = r #this has it only with the rotating bottom and zero everywhere else


#Some testing I did that I kept
j = Heun(v)

v = np.array([v,j])

h = np.zeros((1,50,int(50*Gamma)))
k = Heun(v[-1])
#k[:,0] = r 

h[0] = k

v = np.vstack([v,h])
check(v[-1],[-2])
print(check(v[-1],v[-2]))


# Idea behind method #
#Think of a stack of two cards here. The old frame is the top card and the new frame is the bottom.
#Once I compare the two I throw away the old card and put a new card below the remaining card.
#I count every card and I make a stack of every 50th card 

#The two cards initiallly 
matrix_checker = np.zeros((2,50,int(50*Gamma)), dtype = np.float64)
matrix_checker[-1] = h[0]
matrix_checker[0] = j

#initializer for loop
l = np.zeros((1,50,int(50*Gamma)), dtype = np.float64)

# CHECKING LOOP #
while check(v[-1],v[-2]) == False:
    inter_num += 1

    g = Heun(matrix_checker[0])
    #g[:,0] = r
    l[0] = g
    matrix_checker[-1] = matrix_checker[0]
    matrix_checker[0] = l[0]
    if inter_num % 50 == 0:
        animation_matrix = np.vstack([animation_matrix,l]) 
    
    if inter_num == max_interation -3: #manual kill in case things blow up
        print('Rel error at the end is' , np.abs((l_2(matrix_checker[-1]) - l_2(matrix_checker[-2]))/l_2(matrix_checker[-1])))
        break

    
animation_matrix = np.vstack([animation_matrix,l]) #Final size here 151,50,75
plot_contour(v[-1], 1.5, 'TEST GRAPH')

#error between the exact and last numerical solution
print('Error between last time steps and exact: ', np.abs((l_2(w) - l_2(matrix_checker[-1])))/l_2(w))     

# ANIMATION #
def animate(k):
    plotheatmap(np.transpose(animation_matrix[k]), k)

anim = animation.FuncAnimation(plt.figure(), animate, interval=1, frames=143, repeat=True, save_count=1500) 

# Saving the gif file
save_start = datetime.now()
f = r"c://Users/[COMPUTER'S USER HERE]/Desktop/MAT_462_Time_Evo_MKX_Ver_II_Re_1.gif" 
writergif = animation.PillowWriter(fps=60) 
anim.save(f, writer=writergif)
print('Gif save run time: ', datetime.now()-save_start)


# Total run time
print('total run time: ', datetime.now() - start)
