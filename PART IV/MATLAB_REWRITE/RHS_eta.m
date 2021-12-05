function detadt = RHS_eta(v, psi, eta, dr, dz,Re)
[nr, nz] = size(v);
detadt = zeros(size(v));
for i = 2:nr-1
    for j =2:nz-1
        % Derivatives of eta
        % First derivatives
        eta_z = (eta(i,j+1)- eta(i,j-1))/(2*dz); % z
        eta_r = (eta(i+1,j)- eta(i-1,j))/(2*dr); % r
        
        % Second Derivatives
        eta_zz = (eta(i+1,j) - 2*eta(i,j) + eta(i-1,j))/dr^2; % zz
        eta_rr = (eta(i,j+1) - 2*eta(i,j) + eta(i,j-1))/dz^2; % rr
        
        % Derivatives of psi 
        psi_z = (psi(i,j+1)- psi(i,j-1))/(2*dz); % z
        psi_r = (psi(i+1,j)- psi(i-1,j))/(2*dr); % r
        
        % Derivatives of v
        % First Derivatives 
        v_z = (v(i,j+1)- v(i,j-1))/(2*dr); % z
        v_r = (v(i+1,j)- v(i-1,j))/(2*dr); % r
       
        % Evaluating deta/dt
        detadt(i,j) = (psi_z*eta_r)/(i*dr) - (psi_r*eta_z)/(i*dr) - eta(i,j)/((i*dr)^2)*psi_z+...
                      2*(v(i,j)*v_z)/(i*dr) + 1/Re * (eta_rr + (eta_r)/(i*dr) - eta(i,j)/((i*dr)^2 + eta_zz) );
    end
end
end
