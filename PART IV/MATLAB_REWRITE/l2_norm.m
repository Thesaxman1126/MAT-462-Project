%% MAT 462 PART IV: L2 NORM 

function l2 = l2_norm(array,nr,nz)
    l2 = sqrt(1/((nr+1)*(nz+1)) * sum(array.^2, 'all')); % Definition of discrete l2 norm
end
