% Subroutine for boundary layer calculation
% using Approximation method: Karman-Pohlhausen
% NOTE: *BL calculation only for laminar region
%       *transition checked using Cebeci and Smith (1974) method
%           which is an improvement of Michel's method
%       *turbulent region is neglected

function [delta, deltas, thetas, tw, cf, trans_id, sp] = boundLayer_v2...
    (U_in, xb, yb, xp, yp, c, rho, U_inf, mu)
%% Initialize input variables
if size(U_in, 2) ~= 1
    U_in = U_in';
end
kmu = mu/rho;
Re = U_inf*c/kmu;

% Obtain stagnation point
[~, ind_stag] = min(abs(U_in));
if U_in(ind_stag+1) < U_in(ind_stag-1)
    up = ind_stag + 1;
    low = ind_stag;
else
    up = ind_stag;
    low = ind_stag + 1;
end

sp = [up low];
rp = [xp, yp];

% Calculate ds
rb = [xb, yb];
ds = zeros(length(xp), 1);
for i = 1:length(xp)
    ds(i) = norm(rb(i+1,:) - rb(i,:));
end

%% Handle functions
lambda_fun = @(x) x*(37/315 - x/945 - (x^2)/9072)^2;
H_fun = @(x) (3/10 - x/120)/(37/315 - x/945 - (x^2)/9072);
l_fun = @(x) (2 + x/6)/(37/315 - x/945 - (x^2)/9072);
F_fun = @(x,U) 2*(l_fun(x) - (2 + H_fun(x))*lambda_fun(x));

% ODE functions
syms Af(x)
Z = lambda_fun(Af);
dZ = diff(Z, x);

%% Array Initialization
A = zeros(length(U_in),1);
delta = zeros(length(U_in),1);
deltas = zeros(length(U_in),1);
thetas = zeros(length(U_in),1);
tw = zeros(length(U_in),1);
cf = zeros(length(U_in),1);

%% Solution for UPPER airfoil
id = up:length(U_in);
% Calculate dU_in/dx1
dU_u = zeros(length(U_in(id)),1);
for i=id
    j = i - up + 1;
    if i ~= length(U_in)
        dU_u(j) = (U_in(i+1) - U_in(i))/(norm(rp(i+1,:) - rp(i,:))/c);
    else
        dU_u(j) = (U_in(i) - U_in(i-1))/(norm(rp(i,:) - rp(i-1,:))/c);
    end
end

% Find maximum tangential velocity
[~, umax] = max(abs(U_in(id)));
BC_u = xp(umax);

% Variables initialization: ODE and Boundary Condition
cond = Af(BC_u) == 0;
ku = dU_u./U_in(id);

% For transition check
xt = 0;

% Iteration to calculate handle variables for all position
for i = 1:length(dU_u)
    j = i+up-1;
    
    % Solve ODE
    F_u = F_fun(Af)*ku(i);
    ode = dZ == F_u;
    ASol = dsolve(ode, cond);

    % Calculate A
    G = ASol(xp(j));
    A(j) = G;
     
    % Calculate xt
    xt = xt + ds(j);

    clear temp G F_u ode ASol
end

% Construct the BL result variables
delta(id) = sqrt(A(id)*kmu./dU_u);
deltas(id) = delta(id).*(3/10 - A(id)./120);
thetas(id) = delta(id).*(37/315 - A(id)./945 - (A(id).^2)./9072);
tw(id) = (mu*U_in(id)./delta(id)).*(2 + A(id)./6);
cf(id) = 2*kmu*l_fun(A(id))./(U_in(id).*thetas(id));

% Find the transition point
RHS = abs(U_in(id)).*thetas(id)./kmu;
k = abs(U_in(id)).*xt./kmu;
LHS = 1.174*(1+22400./k).*k.^0.46; % Cebeci and Smith (1974)
trans_u = abs(RHS - LHS);

trans_id(1) = find(trans_u <= 1e-05);

% Clear unused parameters
clear RHS LHS xt id
clear x Au Z ode cond ASol

%% Solution for LOWER airfoil
id = low:-1:1;
% Calculate dU_in/dx1
dU_l = zeros(length(U_in(id)),1);
j = 1;
for i=id
    if i ~= 1
        dU_l(j) = (U_in(i+1) - U_in(i))/(norm(rp(i+1,:) - rp(i,:))/c);
    else
        dU_l(j) = (U_in(i) - U_in(i-1))/(norm(rp(i,:) - rp(i-1,:))/c);
    end
    j = j + 1;
end

% Find maximum tangential velocity
[~, umax] = max(abs(U_in(id)));
BC_l = xp(umax);

% Variables initialization: ODE and Boundary Condition
syms Au(x)
Z = lambda_fun(Au)./dU_l;
ode = diff(Z, x) == F_fun(Au);
cond = Au(BC_l) == 0;
ASol = dsolve(ode, cond);

% For transition check
xt = 0;

% Iteration to calculate handle variables for all position
for i = 1:length(dU_l)
    j = low+1-i;
    
    temp = ASol(xp(j)); % temp
    G = temp(i); % temp

    % Calculate A
    A(j) = G;
    
    % Calculate xt
    xt = xt + ds(j);

    clear temp G
end

% Construct the BL result variables
delta(id) = sqrt(A(id)*kmu./dU_l);
deltas(id) = delta(id).*(3/10 - A(id)./120);
thetas(id) = delta(id).*(37/315 - A(id)./945 - (A(id).^2)./9072);
tw(id) = (mu*U_in(id)./delta(id)).*(2 + A(id)./6);
cf(id) = 2*kmu*l_fun(A(id))./(U_in(id).*thetas(id));

% Find the transition point
RHS = abs(U_in(id)).*thetas(id)./kmu;
k = abs(U_in(id)).*xt./kmu;
LHS = 1.174*(1+22400./k).*k.^0.46; % Cebeci and Smith (1974)
trans_l = abs(RHS - LHS);

trans_id(2) = find(trans_l <= 1e-05);

% Clear unused parameters
clear RHS LHS xt id
clear x Au Z ode cond ASol

end