%% MAIN PROGRAM
clear; close all
% Viscous-Inviscid Interaction Solver
% 1. Inviscid solver using Potential flow method
% 2. Boundary layer approximation method with:
%       2.1. Laminar: Karman-Pohlhausen
%       2.2. Transition: Cebeci and Smith (1974)
%       2.3. Turbulent: Temporarily neglected

% Input values
airfoil_dir = 'naca2308.txt';
c = 1; % [m] chord length

mu = 1.7894e-5; % dynamic viscosity
rho = 1.225;    % [kg/m^3] density

%% Input geometry
[xb, yb, m, mp1, U_inf, alpha] = readData(airfoil_dir);
num_panel = m;
q = 0.5*rho*U_inf^2;
Re = rho*U_inf*c/mu;

%% Viscous-Inviscid Iteration
threshold = 1e-05;
err = 1000;
iter = 0;
% Initialization
old_dels = zeros(num_panel, 1);
g = zeros(num_panel+1, 1);

while err >= threshold
    %% Calculate inviscid term using Potential Flow method
    % TEMPORARY SUBROUTINE
    [x, y, gamma, vtan, cp] = panelMethod(xb,yb,m,mp1,alpha,g);
    
    % The results are U_in, cp, cl
    U_in = abs(U_inf*vtan');

    %% Calculate viscous term using Karman-Pohlhausen BL method
    [delta, deltas, t_wall, cf, trans, stag, uns] = boundLayer_v3...
        (U_in, xb', yb', x', y', c, rho, mu);

    % Calculate boundary condition for inviscid solver
    g(1:end-1) = 0.03.*deltas;
    g(end) = 0.03.*deltas(m);
    
    % Looping criteria
    err = sum(abs((deltas - old_dels)./old_dels)*100);
        
    old_dels(:) = deltas;
    
    % Plot error
    iter = iter + 1;
    num_iter(iter) = iter;
    err_iter(iter) = err;
    %figure(1)
    %plot(num_iter(1:iter),err_iter(1:iter),'r');

    fprintf('iteration %d ------ error = %.2d\n', iter, err);
end

%% Post-Processing
up = stag+1; low = stag-1;
% Plot error iteration
figure; grid on;
plot(num_iter, err_iter, 'r', 'linewidth', 2);
xlabel('Number of iteration'); ylabel('Total error [%]');
title('Error convergence');

%% Calculate lift and drag coefficient
cp_u = 0; cp_l = 0;
cf_u = 0; cf_l = 0;
% Upper airfoil
for i = up:m-1
    cp_u = cp_u + (cp(i+1) + cp(i))*(x(i+1)-x(i))/2;
    if i < uns(1)
        cf_u = cf_u + (cf(i+1) + cf(i))*(x(i+1)-x(i))/2;
    end
end
% Lower airfoil
for i = low:-1:2
    cp_l = cp_l + (cp(i-1) + cp(i))*(x(i-1)-x(i))/2;
    if i > uns(2)
        cf_l = cf_l + (cf(i-1) + cf(i))*(x(i-1)-x(i))/2;
    end
end
cl = (cp_l - cp_u)*cos(alpha*pi/180);
cd = (cf_u + cf_l);
% Displaying the results
disp('Results: ');
fprintf('Cl = %.2f\n', cl);
fprintf('Cd = %.4f\n', cd);

%% Plot pressure coefficient and velocity distribution
ind_u = stag:length(x);
ind_l = 1:stag;

figure; hold on; grid on
plot(x(ind_u), -cp(ind_u), 'b-');
plot(x(ind_l), -cp(ind_l), 'r-');
axis([0 1 min(-cp) max(-cp)]);
xlabel('x/c'); ylabel('-c_p');
title('Coefficient of pressure distribution');
legend('Upper', 'Lower');
hold off

figure; hold on; grid on
plot(x(ind_u), abs(vtan(ind_u)), 'b-');
plot(x(ind_l), abs(vtan(ind_l)), 'r-');
axis([0 1 min(vtan) max(vtan)]);
xlabel('x/c'); ylabel('v/U');
title('Velocity distribution');
legend('Upper', 'Lower');
hold off

% Plot shear at wall
figure; hold on; grid on
plot(x(ind_u), cf(ind_u), 'b-');
plot(x(ind_l), cf(ind_l), 'r-');
xlabel('x/c'); ylabel('cf');
title('Coefficient of friction distribution');
legend('Upper', 'Lower');
hold off

%% Plot airfoil with boundary layer
%boundary layer thickness
yBL = zeros(1,m);
yBL(up:end) = y(up:end) + delta(up:end)';
yBL(1:low) = y(1:low) - delta(1:low)';
yBL(stag) = y(stag) + delta(stag)*sin(alpha+pi);
xBL(stag) = x(stag) - delta(stag)*cos(alpha);

figure; grid on; hold on; axis equal
plot(x, y, 'k');
plot(x(uns(2)+1:uns(1)-1), yBL(uns(2)+1:uns(1)-1), 'r');
plot(xBL(stag), yBL(stag), 'b*');
legend('Airfoil', 'Boundary Layer');
xlabel('x/c'); ylabel('y');
title('Boundary Layer at Airfoil');