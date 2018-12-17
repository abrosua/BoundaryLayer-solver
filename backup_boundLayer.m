%% Subroutine for boundary layer calculation
% using Approximation method: Karman-Pohlhausen
% NOTE: *BL calculation only for laminar region
%       *transition checked using Cebeci and Smith (1974) method
%           which is an improvement of Michel's method
%       *turbulent region is neglected

function [delta_star, theta_star, t_wall, cf, cd] = boundLayer(U_in,...
    xp, yp, c, rho, U_inf, mu)
%% Initialize input variables
temp = 1;
% Obtain stagnation point
[U_stag, ind_stag] = min(abs(U_in));

% Rearrange for upper (u) and lower (l) variables
xc_u = xp(ind_stag:end); xc_l = xp(ind_stag:-1:1);
xc = [xc_u; xc_l];
yc_u = yp(ind_stag:end); yc_l = yp(ind_stag:-1:1);
yc = [yc_u; yc_l];

r_u = [xc_u, yc_u];
r_l = [xc_l, yc_l];
% Calculate dr for both upper and lower airfoil
ds_u = zeros(length(xc_u)-1, 1);
ds_l = zeros(length(xc_l)-1, 1);
for i = 1:length(xc_u)-1
    ds_u(i) = norm(r_u(i+1) - r_u(i));
end
for i = 1:length(xc_l)-1
    ds_l(i) = norm(r_l(i+1) - r_l(i));
end

V_u = U_in(ind_stag:end); V_l = U_in(ind_stag:-1:1);
V = [V_u, V_l];

%% Handle functions
lambda_fun = @(x) x*(37/315 - x/945 - (x^2)/9072)^2;
H_fun = @(x) (3/10 - x/120)/(37/315 - x/945 - (x^2)/9072);
l_fun = @(x) (2 + x/6)/(37/315 - x/945 - (x^2)/9072); 

%% Solution for UPPER airfoil
U_u = V_u(1:end-temp);
% Calculate dU/d
dU_u = zeros(length(U_u)-1,1);
for i=1:length(U_u)-1
    dU_u(i) = (U_u(i+1) - U_u(i))/(norm(r_u(i+1,:) - r_u(i,:))/c);
end
% Calculate d2U/d2
d2U_u = zeros(length(dU_u)-1,1);
for i=1:length(dU_u)-1
    d2U_u(i) = (dU_u(i+1) - dU_u(i))/(norm(r_u(i+1,:) - r_u(i,:))/c);
end

% Variables initialization: A, lambda, H, l, etc
A_u(1) = 7.052; % stable
lambda_u(1) = lambda_fun(A_u(1));
H_u(1) = H_fun(A_u(1));
l_u(1) = l_fun(A_u(1));
F_u(1) = 2*(l_u(1) - (2 + H_u(1))*lambda_u(1))/U_u(1);
Z_u(l) = lambda_u(1)/dU_u(1);
% For transition check
RHS_u(1) = 1;
LHS_u(1) = 0;

% Iteration to calculate handle variables for all position
for j = 2:length(dU_u)-1
    if abs(RHS_u(j-1) - LHS_u(j-1)) >= 1e-08
        % Update variables
        for i = 1:length(dU_u)-1
            Z_u(i+1) = Z_u(i) + F_u(i)*ds_u(i)/c;
            lambda_u(i+1) = dU_u(i+1)*Z_u(i+1);
            if lambda_u(i+1) <=0
                G = fzero(@(L) lambda_fun(L) - lambda_u(i+1), -8);
            elseif lambda_u(i+1) > 0
                G = fzero(@(L) lambda_fun(L) - lambda_u(i+1), 8);
            end
            A_u(i+1) = G;
            
            % Calculate other variables
            H_u(i+1) = H_fun(G);
            l_u(i+1) = l_fun(G);
            F_u(i+1) = 2*(l_u(i+1) - (2 + H_u(i+1))*lambda_u(i+1))/U_u(i+1);
            
            clear G
        end
        % Checking the transition point
        RHS_u(j) = abs(U_u(j))*Re*sqrt(Z_u(j)/Re);
        k = rho*abs(U_u(j))*U_inf*sum(ds_u(1:j))/mu;
        LHS_u(j) = 1.174*(1+22400/k)*k^0.46; % Cebeci and Smith (1974)
    else
        RHS_u(j) = 1;
        LHS_u(j) = 0;
    end
end

% Construct the BL result variables
delta_u = sqrt(A_u./(Re*dU_u(1:length(A_u))));
deltas_u = delta_u.*(3/10 - delta_u./120);
thetas_u = sqrt(Z_u./Re);
tw_u = l.*U_u(1:length(l)).*sqrt(Re./Z_u);

% Calculate cfx_u (coefficient of friction)
cfx_u = 2*tw_u/Re;

%% Solution for LOWER airfoil
U_l = V_l(1:end-temp);
% Calculate dU/d
dU_l = zeros(length(U_l)-1,1);
for i=1:length(U_l)-1
    dU_l(i) = (U_l(i+1) - U_l(i))/(norm(r_l(i+1,:) - r_l(i,:))/c);
end
% Calculate d2U/d2
d2U_l = zeros(length(dU_l)-1,1);
for i=1:length(dU_l)-1
    d2U_l(i) = (dU_l(i+1) - dU_l(i))/(norm(r_l(i+1,:) - r_l(i,:))/c);
end

% Variables initialization: A, lambda, H, l, etc
A_l(1) = 7.052; % stable
lambda_l(1) = lambda_fun(A_l(1));
H_l(1) = H_fun(A_l(1));
l_l(1) = l_fun(A_l(1));
F_l(1) = 2*(l_l(1) - (2 + H_l(1))*lambda_l(1))/U_l(1);
Z_l(l) = lambda_l(1)/dU_l(1);
% For transition check
RHS_l(1) = 1;
LHS_l(1) = 0;

% Iteration to calculate handle variables for all position
for j = 2:length(dU_l)-1
    if abs(RHS_l(j-1) - LHS_l(j-1)) >= 1e-08
        % Update variables
        for i = 1:length(dU_l)-1
            Z_l(i+1) = Z_l(i) + F_l(i)*ds_l(i)/c;
            lambda_l(i+1) = dU_l(i+1)*Z_l(i+1);
            if lambda_l(i+1) <=0
                G = fzero(@(L) lambda_fun(L) - lambda_l(i+1), -8);
            elseif lambda_l(i+1) > 0
                G = fzero(@(L) lambda_fun(L) - lambda_l(i+1), 8);
            end
            A_l(i+1) = G;
            
            % Calculate other variables
            H_l(i+1) = H_fun(G);
            l_l(i+1) = l_fun(G);
            F_l(i+1) = 2*(l_l(i+1) - (2 + H_l(i+1))*lambda_l(i+1))/U_l(i+1);
            
            clear G
        end
        % Checking the transition point
        RHS_l(j) = abs(U_l(j))*Re*sqrt(Z_l(j)/Re);
        k = rho*abs(U_l(j))*U_inf*sum(ds_l(1:j))/mu;
        LHS_l(j) = 1.174*(1+22400/k)*k^.46;
    else
        RHS_l(j) = 1;
        LHS_l(j) = 0;
    end
end

% Construct the BL result variables
delta_l = sqrt(A_l./(Re*dU_l(1:length(A_l))));
deltas_l = delta_l.*(3/10 - delta_l./120);
thetas_l = sqrt(Z_l./Re);
tw_l = l.*U_l(1:length(l)).*sqrt(Re./Z_l);

% Calculate cfx_l (coefficient of friction)
cfx_l = 2*tw_l/Re;

%% Calculate total Cf
if isempty(X) % X, MM and M are unknown! figure it out later!
    cf = sum(abs(cfx_u(1:MM)).*ds_u(1:MM)) + sum(abs(cfx_l).*ds_l(1:length(cfx_bal)));
    Drag = .5*rho*U_inf^2*c*span*cf;
elseif ~isempty(X) ~= 0
    cf = sum(abs(cfx_u(1:MM)).*ds_u(1:MM)) + sum(abs(cfx_l(1:M)).*ds_l(1:M));
    Drag = .5*rho*U_inf^2*c*span*cf;
end

%% Plotting the data

end