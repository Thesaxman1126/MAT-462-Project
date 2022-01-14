## THIS CODE IS NOT FINISHED PLEASE USE THE MATLAB REWRITE WHILE THIS GETS FIXED ##
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 01:16:05 2021

@author: Nick
@class: MAT 462
Numerical Non-Linear flow 
"""

import numpy as np #Used for arrays and matrix like structures
import matplotlib.pyplot as plt #Used for plotting 
from mpl_toolkits.mplot3d import Axes3D #Used for some 3D-plotting 
from datetime import datetime #Used in runtime calculation
from scipy import linalg as la 


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

# Tri Diagonal Matrix Algorithm(a.k.a Thomas algorithm) solver
def TDMAsolver(a, b, c, d):
 
    nf = len(d) # number of equations
    ac, bc, cc, dc = (x.astype(float) for x in (a, b, c, d)) # copy arrays & cast to floats
    for it in range(1, nf):
        mc = ac[it-1]/bc[it-1]
        bc[it] = bc[it] - mc*cc[it-1]
        dc[it] = dc[it] - mc*dc[it-1]
 
    xc = bc
    xc[-1] = dc[-1]/bc[-1]
 
    for il in range(nf-2, -1, -1):
        xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]
 
    return xc

# MATRIX CONSTUCTORS FOR A_nn AND B_mm #
# A_nn 
def Construct_A_nn(Nmax,dr):
    A_nn = np.zeros((Nmax, Nmax), dtype=np.float64)
    #Main
    for i in range(Nmax):
    
        
        A_nn[i,i] = -2/(dr**2) #use i+1 to avoid division by 0 since numpy arrays index from 0 to n-1
    
    #Sub    
    for i in range(Nmax-1):
        
     
        A_nn[i+1,i] = (1+1/(2*(i+1)))/(dr**2)
         
    #Super
    for i in range(Nmax-1):
        
 
        A_nn[i,i+1] = (1-1/(2*(i+1)))/(dr**2)
        
    return A_nn

# B_mm
def Construct_B_mm(Mmax, GAMMA):
    #Z-DIFF tridiag matrix B_mm 
    B_mm = np.zeros((Mmax,Mmax), dtype = np.float64)
    
    #B_mm
    #Main
    for j in range(Mmax):
        
        dz = GAMMA/(Mmax)
        B_mm[j,j] = -2/(dz**2)
        
    #Sub
    for j in range(Mmax-1):
        
        dz = GAMMA/(Mmax)
        B_mm[j+1,j] = 1/(dz**2)
        
    #Super
    for j in range(Mmax-1):
        
        dz = GAMMA/(Mmax)
        B_mm[j,j+1] = 1/(dz**2)

    return B_mm

def Construct_F_nm(eta , dr):
    pass


# ARRAY FUNCTIONS #
# L2 Norm
def l_2(v, nr = 100, nz = int(100*1.5)):
    
    v_norm_l = np.sqrt(1/((nr+1)*(nz+1)) * np.sum(v)**2)
    return v_norm_l

# Right Hand Side made from part III (Shows up in RHS_v and RHS_eta)
def RHS_III(array, Size = 100):
    u = np.zeros((Size,int(1.5*Size)))
    for i in range(1,Size-1):
        for j in range(1,int(1.5*Size)-1):
            v_rr = (array[i+1,j] - 2*array[i,j] + array[i-1,j])
            v_r_on_r = (array[i+1,j] - array[i-1,j]) / (2*i)
            v_on_r2 = array[i,j]/i**2
            v_zz = (array[i,j+1] - 2*array[i,j] + array[i,j-1])
            
            u[i,j] = (1/Re) * (v_rr + v_r_on_r - v_on_r2 + v_zz)/dr**2
    #u[:,0] = r
    return u 

#RHS of the interior of eta_t
def RHS_eta(v, psi, eta, dr = 1/100, dz = 1/100, Size = 100):
    u = np.zeros((Size,int(1.5*Size)))
    for i in range(1,Size-1):
        for j in range(1,int(1.5*Size)-1):
            
            psi_r = (psi[i+1,j] - psi[i-1,j])/(2*dr)
            psi_z = (psi[i,j+1] - psi[i,j-1])/(2*dz)
            
            v_r = (v[i+1,j] - v[i-1,j])/(2*dr)
            v_z = (v[i,j+1] - v[i,j-1])/(2*dz)
            
            eta_r = (eta[i+1,j] - eta[i-1,j])/(2*dr)
            eta_z = (eta[i,j+1] - eta[i,j-1])/(2*dz)
            
            u[i,j] = 1/(i*dr)* (psi_z*eta_r - psi_r*eta_z) - eta[i,j]*psi_z/((i*dr)**2) + (2*v[i,j] * v_z)/(i*dr) 
    return u + RHS_III(eta)

#RHS of the interior of v_t
def RHS_v(v, psi, eta, dr = 1/100, dz = 1/100, Size = 100):
    u = np.zeros((Size,int(1.5*Size)))
    for i in range(1,Size-1):
        for j in range(1,int(1.5*Size)-1):
            
            psi_r = (psi[i+1,j] - psi[i-1,j])/(2*dr)
            psi_z = (psi[i,j+1] - psi[i,j-1])/(2*dz)
            
            v_r = (v[i+1,j] - v[i-1,j])/(2*dr)
            v_z = (v[i,j+1] - v[i,j-1])/(2*dz)
    
            u[i,j] = 1/(i*dr) * psi_z * (v_r + v[i,j]/(i*dr)) - 1/(i*dr) * psi_r * v_z  
    return u + RHS_III(v)
    

# Two step Heun Predictor/Corrector
def Heun_v(v, psi, eta, dt = 1e-4, dr = 1/100, dz = 1/100):
    predict_v = v + RHS_v(v, psi, eta) * dt
    
    corrected_v = v + 0.5 * dt * (RHS_v(v, psi, eta) + RHS_v(predict_v, psi, eta))
    return corrected_v

def Heun_eta(v, psi, eta, dt = 1e-4):
    predict_eta = eta + RHS_eta(v, psi, eta) * dt
    
    corrected_eta = eta + 0.5 * dt * (RHS_eta(v, psi, eta) + RHS_eta(v, psi, predict_eta))
    return corrected_eta


# Makes Contourplots
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

def solve_psi(B_mm, eig, I_mm, Z_nn, inverseZ_nn, F_nm, ):
    for index in range(nr-2):
        H_mn = np.transpose(F_nm).dot(np.transpose(inverseZ_nn))
        U_nm = np.zeros((nr-2, nz-2), dtype=np.float64)
        B_mmEigI = B_mm - eig[index]*I_mm
        (Permutation, Lower_diag, Upper_diag) = la.lu(B_mmEigI)
        
        y_vector = TDMAsolver((Lower_diag).diagonal(-1),(Lower_diag).diagonal(0),(Lower_diag).diagonal(1), H_mn[:,index])
        
        U_nm[index,:] = TDMAsolver((Upper_diag).diagonal(-1),(Upper_diag).diagonal(0),(Upper_diag).diagonal(1), y_vector)
        
    psi_nm = Z_nn.dot(U_nm)
    #print(f"Time to run this step was: {datetime.now() - start_time}")
    for r in range(nr-2):
        for z in range(nz-2):
            psi[r+1][z+1] = psi_nm[r][z]
    
    return psi

# CONSTANTS/INITIALIZERS #

start_time = datetime.now()

#Constants
gamma = 1.5  #Aspect ratio

dr = 1/100 #spacing between grid points for r 
nr = 100 #number of grid points on r-axis

dz = 1/100 #spacing between grid points for z
nz = int(gamma / dz) #number of grid points on z-axis

Re = 1e3 #Reynolds number

t = 0 #Initial time
dt = 1e-4 #Change in time

r = np.linspace(0,1,100)

max_interation = 1500 #max interation time

#Matrix initialization
# v matricies #
v = np.zeros((nr,nz))
v[:,0] = r 

pain = l_2(v)
# Psi matricies #
psi = np.zeros((nr,nz))

#psi_nm
psi_nm = np.zeros((nr-2,nz-2))

#eta
eta_ij = np.zeros((nr,nz))
#Test = RHS_eta(v,psi,eta_ij)

v = Heun_v(v, psi, eta_ij)
#v = Heun_v(v, psi, eta_ij)


eta_ij = Heun_eta(v, psi, eta_ij)



#F_nm
#Intialize F_nm
F_nm = np.zeros((nr-2,nz-2))
F_nm = RHS_eta(v,psi,eta_ij)[1:-1,1:-1]


# Construct Matrices #
#A_nn
A_nn = Construct_A_nn(nr-2,dr)

#Z_nn
eig, Z_nn = np.linalg.eig(A_nn)

#Z_nn^(-1)
inverseZ_nn = np.linalg.inv(Z_nn) 

#B_mm
B_mm = Construct_B_mm(nz-2, gamma)

#I_mm
I_mm = np.identity(nz-2)

#Solving the system (B - e *I) u = h#


psi = solve_psi(B_mm, eig, I_mm, Z_nn, inverseZ_nn, F_nm)

print((l_2(v) - pain)/l_2(v))
print()
test = 0
while test < 100:
    v = Heun_v(v, psi, eta_ij)
    #v = Heun_v(v, psi, eta_ij)
    
    ## update bounds on eta
    #right

    for j in range(nz):
        eta_ij[nr-1,j] = 1/(2*dr**2) * (psi[nr-3,j] - 8*psi[nr-2,j])
    
    #top
    for i in range(nr):
        eta_ij[i,nr-1] = 1/(2*(i+1)*dr * dz**2) * (psi[i,nz-3] - 8*psi[i,nz-2])
    
    #bottom
    for i in range(nr):
        eta_ij[i,0] = 1/(2*(i+1)*dr * dz**2) * (psi[i,2] - 8*psi[i,1])
    eta_ij = Heun_eta(v, psi, eta_ij)
    
    F_nm = eta_ij[1:-1,1:-1]
    psi = solve_psi(B_mm, eig, I_mm, Z_nn, inverseZ_nn, F_nm)
    
    

    test += 1
    
    
    
