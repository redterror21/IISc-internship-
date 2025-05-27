clc
clear

% Grid size
x = 32;
y = 32;
Lx = 1;
Ly = 1;
dx = Lx / x;
dy = Ly / y;

% Temperature field
T = zeros(y+2, x+2);  % include ghost cells


% SOR setup
Ap = zeros(y+2, x+2);
Ae = ones(y+2, x+2) / dx^2;
Aw = ones(y+2, x+2) / dx^2;
An = ones(y+2, x+2) / dy^2;
As = ones(y+2, x+2) / dy^2;

Aw(2:end-1,1)     = 0; % Left wall
Ae(2:end-1,end-1) = 0; % Right wall
An(end-1,2:end-1) = 0; % Top wall
As(end-1,2:end-1)     = 0; % Bottom wall

Ap = -(Ae + Aw + An + As);


%%
function [p,err]=sor_solver(p, S, Ap, Ae, Aw, An, As, x, y)
    pk = zeros(size(p));
    it = 0;
    err = 1e10;
    tol = 1e-8;
    maxit=10000;
    B = 1.5;

    while err > tol && it < maxit
        pk = p;
        for i =2:x+1
            for j =2:y+1
                             
               % Boundary Conditions (Dirichlet and Neumann)
               p(2,:) = 600;           % Top wall: cold
               p(:,end-1) = 900;       % Right wall: hot
                              
               % Neumann BCs (insulated)  dT/dx = 0
               p(:,1)   = p(:,2);      % Left
               p(end-1,:)=p(end-2,:);  % Bottom


               ap = Ap(j,i); ae = Ae(j,i); aw = Aw(j,i); an = An(j,i); aso = As(j,i);

               pe = p(j,i+1); pw = p(j,i-1); pn = p(j+1,i); ps = p(j-1,i);

               res = S(j,i) - (ae*pe + aw*pw + an*pn + aso*ps);
               p(j,i) = B * res / ap + (1-B) * pk(j,i);
            end
        end
        err = norm(p(:) - pk(:), 2);
        it = it+1;
    end
end
%%

% Solve Laplace's equation using SOR
S = zeros(y+2, x+2);  % RHS is zero

[T, err] = sor_solver(T, S, Ap, Ae, Aw, An, As, x, y);



% Plotting the solution
[X, Y] = meshgrid(linspace(0, Lx, x), linspace(0, Ly, y));
contourf(X, Y, flipud(T(2:end-1, 2:end-1)), 20);
colorbar;
title('Steady-State 2D Heat Conduction');
xlabel('x');
ylabel('y');
axis equal tight;
