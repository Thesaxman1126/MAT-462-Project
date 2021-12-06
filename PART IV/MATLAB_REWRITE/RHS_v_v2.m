%% MAT 462 PART IV: RHS OF dv/dt 

function dvdt = RHS_v_v2(v,psi,dr,dz,Re)% This is the RHS of equation 34 discretetized 
[nr, nz] = size(v);
dvdt = zeros(size(v));
% This works on INTERIOR POINTS so 2->nr-1 and 2->nz-1
for i = 2:nr-1
    for j = 2:nz-1
       % Derivatives of psi
       psi_z = (psi(i,j+1)- psi(i,j-1))/(2*dz); % z
       psi_r = (psi(i+1,j)- psi(i-1,j))/(2*dr); % r
       
       % Derivatives of v
       % First Derivatives 
       v_z = (v(i,j+1)- v(i,j-1))/(2*dr); % z
       v_r = (v(i+1,j)- v(i-1,j))/(2*dr); % r
       
       % Second Derivatives
       v_zz = (v(i+1,j) - 2*v(i,j) + v(i-1,j))/dr^2; % zz
       v_rr = (v(i,j+1) - 2*v(i,j) + v(i,j-1))/dz^2; % rr
       
       % FULL RHS
       dvdt(i,j) = 1/(i*dr) * psi_z*(v_r + v(i,j)/(i*dr)) - 1/(i*dr) * psi_r*v_z +...
           1/(Re) * (v_rr + v_r/(i*dr) - v(i,j)/((i*dr)^2) + v_zz);
       
    end
end
end
