%% Project Part IV
clear all; % Clear all variables
close all; % Close all figures 
clc; % Clear console
tic % Starts timer

%% Constants and givens

Re = 1; % Reynolds Number
Gamma = 1.5; % Aspect ratio
N = 100; nr = N; nz = Gamma*N; % Basis for grid
dr = 1/100; dz = dr;           % Spacing in r and z
dt = 1e-4/Re;                  % Time step 
r = linspace(0,1,N);           % This is the discrete rotating bottom

%% Array intialization for v, eta, psi

v = zeros(nr,nz); % Velocity
v(:,1) = r; % This is the boundary condition v(r,0,t) = r for r \in [0,1]
psi_nm = zeros(nr,nz); % Streamfn
eta_nm = zeros(nr,nz); % Vorticity 

%% A_nn Matrix

A_nn = zeros(nr-2);                                  % Initalizing A_nn matrix
a_super = zeros(1,nr-3); a_sub = a_super;            % Initalizing diagonals as vectors

for i=1:nr-3                                         % Populate superdiagonal
    a_super(i) = (1-1/(2*i))/(dr^2); 
end 

a_main = -2/(dr^2) + zeros(1,nr-2);                  % Populate main diagonal

for i=1:nr-3                                         % Populate subdiagonal
    a_sub(i) = (1+1/(2*i))/(dr^2); 
end 

A_nn = diag(a_main)+diag(a_sub,-1)+diag(a_super,1);  % Use the vectors to form diagonals in the array

%% B_nn Matrix

B_mm = zeros(nz-2);                                  % Initalize B_mm
b_main = -2/dz^2 + zeros(1,nz-2);                    % Define main diagonal values since it's const
b_super_sub = 1/dz^2 + zeros(1,nz-3);                % Define the super/sub diagonal values since they're const

B_mm = diag(b_main)+diag(b_super_sub,-1)+diag(b_super_sub,1);

%% Z_nn Martices and Eigenvalue Array

[Z_nn,E_nn] = eig(A_nn); % Z_nn eigenvector matrix, E_nn diaonal eigenvalue matrix
eigenvalues = diag(E_nn); % Stores eigenvales in a vector
inv_Z_nn = inv(Z_nn); % Inverse Z_nn

%% Misc Matricies used for A_nn Psi_nm + Psi_nm B_mm = F_nm

I_mm = eye(nz-2); % Identity matrix
F_nm = zeros(nr-2,nz-2); % Initialize F_nm
U_nm = zeros(nr-2,nz-2); % Initialize U_nm

%% LU decomp

% This will 3D arrays as a way to have L U and P for ALL eigenvalues
lower_diag = zeros(nz-2,nz-2,nr-2); % Lower diagonal matrix L
upper_diag = zeros(nz-2,nz-2,nr-2); % Upper diagonal matrix U

for i = 1:nr-2 %LU decompose ance store
    
    [lower_diag(:,:,i),upper_diag(:,:,i)] = lu(B_mm + eigenvalues(i) * I_mm);     
end

%% Huen loop
index_i = 0; index_j = 0; Rel_error = 1;
while Rel_error > 10e-5
    
    
end
toc
