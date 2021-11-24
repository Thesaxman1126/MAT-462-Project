%% Used to make contour plots for figures 1 - 3 %%

clear all
N = 100; % Number of points evaluated
gamma1 = 0.5;
gamma2 = 1.0;
gamma3 = 2.5; % Set cylinder aspect ratio
% Equispaced grid
r = linspace(0,1,N);
z1 = linspace(0, gamma1, gamma1*N);
z2 = linspace(0, gamma2, gamma2*N);
z3 = linspace(0, gamma3, gamma3*N);


[R1,Z1] = meshgrid(r,z1);
[R2,Z2] = meshgrid(r,z2);
[R3,Z3] = meshgrid(r,z3);

V1 = v(R1,Z1,gamma1);
V2 = v(R2,Z2,gamma2);
V3 = v(R3,Z3,gamma3);

figure(1)
contourf(V1)
set(gca,'DataAspectRatio',[1 1 1], 'XTick',0:0:0,'YTick',0:0:0);

figure(2)
contourf(V2)
set(gca,'DataAspectRatio',[1 1 1], 'XTick',0:0:0,'YTick',0:0:0);

figure(3)
contourf(V3)
set(gca,'DataAspectRatio',[1 1 1], 'XTick',0:0:0,'YTick',0:0:0);

function v = v(r,z,gamma)
    v = f(r,z,gamma) + g(r,z,gamma,200*gamma);
end

function f = f(r,z,gamma)
    f = r.*(1-z./gamma);
end

function g = g(r,z,gamma,Nmax)
    C = -2/pi;
    g = 0;
    for n = 1:Nmax
        w = n*pi/gamma;
        zarg = w.*z;
        rarg = w.*r;
        g = g + C*(1/n)*sin(zarg).*besseli(1,rarg)/besseli(1,w);
    end
end
