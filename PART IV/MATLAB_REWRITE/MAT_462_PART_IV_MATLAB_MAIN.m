%% MAT 462 PART IV: MAIN METHOD

clear all; % Clear all variables
close all; % Close all figures 
clc; % Clear console
tic % Starts timer

%% Constants and givens

Re = 1000; % Reynolds Number
Gamma = 1.5; % Aspect ratio
N = 100; nr = N; nz = Gamma*N; % Basis for grid
dr = 1/100; dz = dr;           % Spacing in r and z
dt = 1e-1/Re;                  % Time step 
r = linspace(0,1,N);           % This is the discrete rotating bottom
z = linspace(1,Gamma,Gamma*N);
%% Array intialization for v, eta, psi

v_k = zeros(nr,nz); % Velocity
v_k(:,1) = r; % This is the boundary condition v(r,0,t) = r for r \in [0,1]
psi_k = zeros(nr,nz); % Streamfn
eta_k = zeros(nr,nz); % Vorticity 

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
lower_diag = zeros(nz-2,nz-2,nr-2); % Lower diagonal matrix L initializer
upper_diag = zeros(nz-2,nz-2,nr-2); % Upper diagonal matrix U initializer
permutation_mat = zeros(nz-2,nz-2,nr-2); % Perumtation matrix initializer
for i = 1:nr-2 %LU decompose 
    
    [lower_diag(:,:,i),upper_diag(:,:,i),permutation_mat(:,:,i)] = lu(B_mm + eigenvalues(i) * I_mm);     
end

%% Huen loop

% This will have hierarchy of solving in this order
% begin loop
% V -> ETA -> PSI -> check -> repeat (if check fails)
% end loop 

%$ Initializers
% indexes and error
i_index = 1; j_index = 1; Rel_error = 1; % note Rel_error is non zero to start since 0< 10e-5
time_step_counter = 0; % Used to count how many steps this will take

% Predictor and corrected matricies
% Predictors
v_predict = zeros(nr,nz);
eta_predict = zeros(nr,nz);

% Correctors
v_correct = zeros(nr,nz);
eta_correct = zeros(nr,nz);

% Main loop
while Rel_error > 1e-5
    %% v predict/correct
    v_predict = v_k + dt * RHS_v_v2(v_k, psi_k, dr, dz, Re); % Predict
    v_correct = v_k + dt/2 * (RHS_v_v2(v_k, psi_k, dr, dz, Re) + RHS_v_v2(v_predict, psi_k, dr, dz, Re)); % Correct
    
    % Don't store the v_correct as v_k yet to use since we need to use the
    % same v_k for eta_k
    
    %% eta predict/correct
    eta_predict = v_k + dt * RHS_eta_v2(v_k, psi_k, eta_k, dr, dz, Re); % Predict
    eta_correct = v_k + dt/2 * (RHS_eta_v2(v_k, psi_k, eta_k, dr, dz, Re) + RHS_eta_v2(v_k, psi_k, eta_predict, dr, dz, Re)); % Correct  
    
    %% Solve Poisson's Equation using the A Psi + Psi B = F method 
    % Populate F_nm
    for  i = 1:nr-2
        for j = 1:nz-2
            
            F_nm(i,j) = -i*dr*eta_k(i+1,j+1);
            
        end
    end
    % Construct H_mn matrix
    H_mn = transpose(F_nm)*transpose(inv_Z_nn);
    
    % LU solver loop
    % METHOD L*U*u = h -> let y = U*u -> L*y = h, solve for y -> U*u = y ->
    % solve for u
    for i = 1:nr-2
       
        y_vector = lower_diag(:,:,i)\(permutation_mat(:,:,i) * H_mn(:,i)); % Let y = U*u -> L*y - h, solve for
        U_nm(i,:) = upper_diag(:,:,i)\y_vector; % U*u = y -> solve for u
        
    end
    
    % Populate psi (NOTE THE BOUNDARY IS 0 ON ALL WALLS SO ONLY UPDATES INTERIOR)
    psi(2:nr-1,2:nz-1) = Z_nn * U_nm;
    
    %% Update eta at the walls
    
    % Side wall 
    for j = 1:nz
       
        eta_correct(nr,j) = (psi_k(nr-2,j) - 8*psi_k(nr-1,j))/(2*i*dr^3);
        
    end
    
    % Top and bottom walls
    for i = 1:nr
       
        eta_correct(i,nz) = (psi_k(i,nz-2)-psi_k(i,nz-1))/(2*i*dr*dz^2); % Top wall
        eta_correct(i,1) = (psi_k(i,2)-psi_k(i,1))/(2*i*dr*dz^2); % Bottom Wall
        
    end
    
    %% Check Error
    if rem(i_index,50) == 0 % check every 50 interations
        Rel_error = (l2_norm(v_correct,nr,nz) - l2_norm(v_k,nr,nz))/(l2_norm(v_correct,nr,nz));
        error(j_index) = Rel_error;
        
        %% Plot Contour Plots
        
        % Plot v(r,z)
        figure();
        contourf(transpose(v_k))
        set(gca,'DataAspectRatio',[1 1 1], 'XTick',0:0:0,'YTick',0:0:0);
        saveas(gcf,sprintf('%.0f.png',j_index))
        
        F(j_index) = getframe;
        j_index = j_index + 1; 
        
    end
    
    v_k = v_correct;
    eta_k = eta_correct;
    
    i_index = i_index + 1;
    time_step_counter = time_step_counter +1;
end
toc
