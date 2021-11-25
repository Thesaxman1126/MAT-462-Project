# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 16:55:04 2021

@author: Nick Navas
@class: MAT 462 

Used to generate figures
    
"""
# IMPORTS #
import scipy.special #Used for Bessel I 
import numpy as np #Used for arrays and matrix like structures
import matplotlib.pyplot as plt #Used for plotting 
from mpl_toolkits.mplot3d import Axes3D #Used for some 3D-plotting 
from datetime import datetime #Used in runtime calculation

# FUNCTIONS #
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
 


#Analytic solution: used for testing/comparison/troubleshooting purposes
def v_analytic(r,z,gamma):
    """This here is used for comparison and testing purposes
    Arguments
    ---------
    r : coordinate in r
    z : coordinate in z
    gamma : aspect ratio
    """
    
    a= r*(1-z/gamma) 
    sums = 0
    n = 0
    for n in range(1,26):
        sums += ((scipy.special.iv(1,(n*np.pi*r)/gamma))/(n*scipy.special.iv(1,(n*np.pi)/gamma)))*np.sin(n*np.pi*z/gamma)    
    
    return a - (2/np.pi)*sums

#Contour plot: Creates countour plot
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
    
    if filename:
        fig.savefig(filename)
        plt.close(fig)
        return filename
    else:
        return ax

# BRAINS OF MY OPERATION #
###############################################################################

# CONSTANTS/INITIALIZERS #
Gamma = [0.5,1,2.5] #[0.5,1,2.5] #You can add to this for different gammas
Nmax = 100 #Size of the nxn matrix
start = datetime.now()

for GAMMA in Gamma: 
#while GAMMA == 1:
    loop_start = datetime.now()
    Mmax = int(Nmax * GAMMA) #Sixe of mxm matrix

    #Array for exact answer
    v_exact = np.zeros((Nmax+2,Mmax+2), dtype = np.float64 ) #Create array 
    
    #Popuate the v_exact array 
    for r in range(Nmax+2): 
        for z in range(Mmax+2):
            
            v_exact[r][z] = v_analytic(1/(Nmax+1)*r, 1/(Mmax+1)*z*GAMMA, GAMMA)
            
    #R-Diff tridiag matrix A_nn
    A_nn = np.zeros((Nmax, Nmax), dtype=np.float64)
    
    # WRITE DIAGONALS #
    #A_nn
    #Main
    for i in range(Nmax): #Primary
    
        dr = 1/(Nmax)
        A_nn[i,i] = -(2+(1/((i+1)**2)))/(dr**2) #use i+1 to avoid division by 0 since numpy arrays index from 0 to n-1
    
    #Sub    
    for i in range(Nmax-1):
        
        dr = 1/(Nmax)
        A_nn[i+1,i] = (1-1/(2*(i+1)))/(dr**2)
         
    #Super
    for i in range(Nmax-1):
        
        dr = 1/(Nmax)
        A_nn[i,i+1] = (1+1/(2*(i+1)))/(dr**2)
        
    #Z_nn
    eig, Z_nn = np.linalg.eig(A_nn)
    
    inverseZ_nn = np.linalg.inv(Z_nn) #Z^(-1)
    
    #F_nm
    F_nm = np.zeros((Nmax, Mmax), dtype=np.float64)
    
    for n in range(Nmax):
        
        dz = GAMMA / Mmax
        dr = 1/(Nmax)
        F_nm[n,0] = - (n+1) * (dr / (dz**(2)))
    
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
    
    
    #H_mn array
    H_mn = np.transpose(F_nm).dot(np.transpose(inverseZ_nn))
    
    #Identiy of size m x m
    I_mm = np.identity(Mmax)
    
    #U_nm array
    U_nm = np.zeros((Nmax, Mmax), dtype=np.float64)
    solve_time = datetime.now()
    for i in range(Nmax):
        
        U_nm[i,:] = TDMAsolver((B_mm + eig[i]*I_mm).diagonal(-1),(B_mm + eig[i]*I_mm).diagonal(0),(B_mm + eig[i]*I_mm).diagonal(1), H_mn[:,i])
    print(rf'MK 8 Solving Tridiagonal matrix time for gamma = {GAMMA}: ', datetime.now()-solve_time,'\n')
    #V_nm = Z_nn * U_nm 
    v_numerical = Z_nn.dot(U_nm)
    
    
    
    w = np.zeros((Nmax+2,Mmax+2), dtype = np.float64)
    for i_index in range(Nmax+2):
        w[i_index,0] = 1/(Nmax+1)*i_index
#TEST
    for r in range(Nmax):
        for z in range(Mmax):
            
            w[r+1,z+1] = v_numerical[r][z]
    #plot the functions using defined function above
    plot_contour(w, GAMMA, rf'Numerical Solution for $\Gamma$ = {GAMMA}')
    plot_contour(v_exact, GAMMA, rf'Analytic Solution for $\Gamma$ = {GAMMA}')
    
    #Create an error array
    Error_MAT  = np.zeros((Nmax+2,Mmax+2), dtype = np.float64 )
    
    #Populate error array with (v_exact[r][z] - v_numerical[r][z])^2
    for r in range(Nmax+2):
        for z in range(Mmax+2):
    
            Error_MAT[r][z] = np.abs((v_exact[r][z]-w[r][z])/v_exact[r][z])
    
    #l_2_norm. np.sum() sums over whole array
    l_2_norm_exact = np.sqrt(1/((Nmax+3) * (Mmax+3))*np.sum(np.power(v_exact,2)))
    l_2_norm_approx = np.sqrt(1/((Nmax+3) * (Mmax+3))*np.sum(np.power(w,2)))
    l_2_norm_error = (l_2_norm_exact - l_2_norm_approx)/l_2_norm_exact

    #Plot Error array in case you want a visual to see where the error is worse 

    
    plot_contour((Error_MAT), GAMMA, rf'Error Plot for $\Gamma$ = {GAMMA}')
    
    #print(rf"The value of v(1,0) for gamma = {GAMMA} is: ",v_numerical[-1][0],'\n') #Prints value at v(1,0).
    print(rf"The rel error for gamma = {GAMMA} is: ", l_2_norm_error,'\n') #Prints rel error
    print(rf"The l2_approx for gamma = {GAMMA} is: ", l_2_norm_approx,'\n') #Prints L_2 norm.
    print(rf'The run time for $\Gamma$ = {GAMMA}:', datetime.now() - loop_start, '\n')
    
print("Total run time:", datetime.now() - start)

###############################################################################
