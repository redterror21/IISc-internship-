clc
clear

% grid and domain properties while using staggered grid 
% defining no. of points inside boundary

%% 

x=64 ; % no. of actual grid points in consideration (no ghost cells)
y=64 ;
Lx=1 ;
Ly=0.2 ;
dx=Lx/x ;
dy=Ly/y ;


visc=0.01 ;
%% 

% Problem: Top wall moving, bottom wall fixed : couette flow
Ut=10;
Ub=0 ;

% choose timestep size based on advection-diffusion stability criteria % USING Stability conditions
% rule 1:
dt1 = 0.5/(visc*(1/(dx^2) + 1/(dy^2))) ;  % For Diffusive part von Neumann stability analysis (done on explicit time integration) to check numerical instability 
% to ensure diffusion doesn't propagate too quicly over timesteps

% rule 2:
%For advection part, wheere we let the viscous part of time to dominate  
CFL = 0.5;  % or less
u_max = max(abs([Ut,Ub]));
% v_max =max(abs(v));
% dt2=1000; % randomly put to make d1 as lowest term during diffusion flows
dt2 = CFL * min(dx/u_max );


dt = min(dt1, dt2) % so satisafies both advection and diffusion rule

%% 

% pre-allocated fields         
%time advance code
p = zeros(y+2, x+2);    %y+2,x+2 to add ghost cells
u = zeros(y+2, x+2);    % same mesh just displaced by half cell
v = zeros(y+2, x+2);

% from Solution Algorithm 
ut = zeros(y+2, x+2); % u tilda velocities values
vt = zeros(y+2, x+2); % v tilda velocities values
divut = zeros(y+2, x+2);


% poission solver SOR
% Build the coefficient arrays 

function [p,err]=sor_solver(p, S,Lx,Ly,x, y)
    dx=Lx/x ;
    dy=Ly/y ;
    Ae = ones(y+2, x+2) / dx^2;
    Aw = ones(y+2, x+2) / dx^2;
    An = ones(y+2, x+2) / dy^2;
    As = ones(y+2, x+2) / dy^2;
    Ap=-(Ae+Aw+An+As);
    
    it = 0;
    err = 1e10;
    tol = 1e-8;
    maxit=10000;
    B = 1.83;

    while err > tol && it < maxit
        pk = p;
        for i =2:x+1
            for j =2:y+1
               ap = Ap(j,i); ae = Ae(j,i); aw = Aw(j,i); an = An(j,i); as = As(j,i);

               pe = p(j,i+1); pw = p(j,i-1); pn = p(j+1,i); ps = p(j-1,i);

               res = S(j,i) - (ae*pe + aw*pw + an*pn + as*ps);
               p(j,i) = B * res / ap + (1-B) * pk(j,i);
            end
        end
        err = norm(p(:) - pk(:), 2);
        it = it+1;
    end
end
 
% velocities at cell centers
uc = 0.5 * (u(2:end-1, 2:end-1) + u(2:end-1, 3:end));
vc = 0.5 * (v(2:end-1, 2:end-1) + v(3:end, 2:end-1));
[X, Y] = meshgrid(dx/2:dx:Lx - dx/2, dy/2:dy:Ly - dy/2);

fig = figure('Name', 'Flow Simulation', 'NumberTitle', 'off');

subplot(1, 2, 1);
hQuiver = quiver(X, Y, uc, vc, 'AutoScaleFactor', 3);

title('Velocity Field');
xlabel('x'); ylabel('y');
axis equal tight;

% Initial divergence calculation
divu = zeros(size(p));
for i = 2:x+1
    for j = 2:y+1
        divu(j,i) = (u(j,i) - u(j,i-1)) / dx + (v(j,i) - v(j-1,i)) / dy;
    end
end

subplot(1, 2, 2);
hIm = imagesc(linspace(0,Lx,x), linspace(0,Ly,y), flipud(divu(2:end-1, 2:end-1)));
colorbar;
title('Divergence of Velocity');
xlabel('x'); ylabel('y');
axis equal tight;
grid on;

t=0;
tsteps=1000;

for n = 1:tsteps

    [u, v] = apply_boundary_conditions(u, v, Ut, Ub);


    % x momentum equation
    for i = 3:x+1 
        for j = 2:y+1 
            ue = 0.5 * (u(j, i+1) + u(j, i));
            uw = 0.5 * (u(j, i) + u(j, i-1));
            un = 0.5 * (u(j+1, i) + u(j, i));
            us = 0.5 * (u(j, i) + u(j-1, i));
            vn = 0.5 * (v(j+1, i-1) + v(j+1, i));
            vs = 0.5 * (v(j, i-1) + v(j, i));
            convection = -(ue^2 - uw^2)/dx - (un*vn - us*vs)/dy;
            diffusion = visc * ((u(j, i-1) - 2*u(j, i) + u(j, i+1))/dx^2 + (u(j-1, i) - 2*u(j, i) + u(j+1, i))/dy^2);
            % ut(j, i) = u(j, i) + dt * (convection + diffusion);
            ut(j, i) = u(j, i) + dt * (diffusion);
        end
    end

    % y momentum equation
    for i = 2:x+1
        for j = 3:y+1
            ve = 0.5 * (v(j, i+1) + v(j, i));
            vw = 0.5 * (v(j, i) + v(j, i-1));
            ue = 0.5 * (u(j, i+1) + u(j-1, i+1));
            uw = 0.5 * (u(j, i) + u(j-1, i));
            vn = 0.5 * (v(j+1, i) + v(j, i));
            vs = 0.5 * (v(j, i) + v(j-1, i));
            convection = -(ue*ve - uw*vw)/dx - (vn^2 - vs^2)/dy;
            diffusion = visc * ((v(j, i+1) - 2*v(j, i) + v(j, i-1))/dx^2 + (v(j+1, i) - 2*v(j, i) + v(j-1, i))/dy^2);
            vt(j, i) = v(j, i) + dt * (convection + diffusion);
        end
    end

    % computing divergence of intermediate velocity
    rho=1; %density
    divut(2:end-1, 2:end-1) = (ut(2:end-1, 3:end) - ut(2:end-1, 2:end-1)) / dx + (vt(3:end, 2:end-1) - vt(2:end-1, 2:end-1)) / dy;
    rhs = rho*divut / dt;

    [p, ~] = sor_solver(p, rhs,Lx,Ly, x, y);
    % zero-gradient boundary condition for pressure at boundaries
    p(:,1)     = p(:,end-1);        % left
    p(:,end)   = p(:,2);            % right
    

    
    %velocity correction
    u(2:end-1, 3:end-1) = ut(2:end-1, 3:end-1) - dt * (p(2:end-1, 3:end-1) - p(2:end-1, 2:end-2)) / dx;
    v(3:end-1, 2:end-1) = vt(3:end-1, 2:end-1) - dt * (p(3:end-1, 2:end-1) - p(2:end-2, 2:end-1)) / dy;

    uc = 0.5 * (u(2:end-1, 2:end-1) + u(2:end-1, 3:end));
    vc = 0.5 * (v(2:end-1, 2:end-1) + v(3:end, 2:end-1));
    
    
    if mod(n,20)==0

        % update velocity quiver plot
        set(hQuiver, 'UData', uc, 'VData', vc);
        subplot(1, 2, 1);
        title(['Velocity Field at time step ', num2str(n)]);

        % update divergence plot
        set(hIm, 'CData', flipud(divu(2:end-1, 2:end-1)));
        subplot(1, 2, 2);
        title(['Divergence at time step ', num2str(n)]);
    
        drawnow;
    end

    t = t + dt;
end

fprintf('Max divergence: %.2e\n', max(abs(divu(:))));


function [u, v] = apply_boundary_conditions(u, v, Ut, Ub)
    % Left wall
    u(:,1) = u(:,2);
    v(:,1) = v(:,2);

    % Right wall
    u(:,end) = u(:,end-1);
    v(:,end) = v(:,end-1);

    % Top wall
    u(end-1,:) = Ut;
    v(end,:) = 0;

    % Bottom wall
    u(1,:) = Ub;
    v(1,:) = 0;

end


