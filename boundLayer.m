% Subroutine for boundary layer calculation
% using Approximation method: Karman-Pohlhausen
% NOTE: *BL calculation only for laminar region
%       *transition checked using Cebeci and Smith (1974) method
%           which is an improvement of Michel's method
%       *turbulent region is neglected

function [delta, deltas, thetas, tw, cf, trans_id, sp] = boundLayer...
    (U_in, xb, yb, xp, yp, c, rho, U_inf, mu)
%% Initialize input variables
if size(U_in, 2) ~= 1
    U_in = U_in';
end
kmu = mu/rho;
Re = U_inf*c/kmu;

% Obtain stagnation point
[~, ind_stag] = min(abs(U_in));
disp(ind_stag);

if abs(U_in(ind_stag+1)) < abs(U_in(ind_stag-1))
    up = ind_stag + 1;
    low = ind_stag;
else
    up = ind_stag;
    low = ind_stag - 1;
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
lambda_fun = @(x) x.*(37/315 - x/945 - (x.^2)/9072).^2;
H_fun = @(x) (3/10 - x/120)./(37/315 - x/945 - (x.^2)/9072);
l_fun = @(x) (2 + x/6)./(37/315 - x/945 - (x.^2)/9072);

%% Array Initialization
A = zeros(length(U_in),1);
delta = zeros(length(U_in),1);
deltas = zeros(length(U_in),1);
thetas = zeros(length(U_in),1);
tw = zeros(length(U_in),1);
cf = zeros(length(U_in),1);
trans_id = zeros(2,1);

%% Solution for UPPER airfoil
init = up;
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
d2U_u1 = (dU_u(2) - dU_u(1))/(norm(rp(init+1,:) - rp(init,:))/c);

% Variables initialization: A, lambda, H, l, etc
A(init) = 7.052; % stable
lambda(init) = lambda_fun(A(init));
H(init) = H_fun(A(init));
l(init) = l_fun(A(init));
F(init) = -.0652857*d2U_u1/dU_u(1)^2;
%F(init) = 2*(l(init) - (2 + H(init))*lambda(init))/U_in(init);
Z(init) = lambda(init)/dU_u(1);

% For transition check
xt = zeros(length(dU_u),1);

% Iteration to calculate handle variables for all position
for i = 1:length(dU_u)-1
    j = i+up-1;
    Z(j+1) = Z(j) + F(j)*ds(j)/c;
    lambda(j+1) = dU_u(i+1)*Z(j+1);
    if lambda(j+1) <= 0
        G = fzero(@(L) lambda_fun(L) - lambda(j+1), -8);
    elseif lambda(j+1) > 0
        G = fzero(@(L) lambda_fun(L) - lambda(j+1), 8);
    end
    A(j+1) = G;

    % Calculate other variables
    H(j+1) = H_fun(G);
    l(j+1) = l_fun(G);
    F(j+1) = 2*(l(j+1) - (2 + H(j+1))*lambda(j+1))/U_in(j+1);

    clear G

    % Calculate xt
    xt(i+1) = sum(ds(init+1:j+1));
end
    
% Construct the BL result variables
delta(id) = sqrt(abs(A(id)*kmu./dU_u));
deltas(id) = delta(id).*(3/10 - A(id)./120);
thetas(id) = delta(id).*(37/315 - A(id)./945 - (A(id).^2)./9072);
tw(id) = (mu*U_in(id)./delta(id)).*(2 + A(id)./6);
cf(id) = 2*kmu*l_fun(A(id))./(U_in(id).*thetas(id));

% Checking the transition point
RHS = abs(U_in(id)).*thetas(id)./kmu;
k = abs(U_in(id)).*xt./kmu;
LHS = 1.174*(1+22400./k).*k.^0.46; % Cebeci and Smith (1974)
trans_u = abs(RHS - LHS);

temp = find(trans_u <= 1e-05);
if isempty(temp)
    trans_id(1) = nan;
else
    trans_id(1) = temp;
end

clear RHS LHS xt id temp

%% Solution for LOWER airfoil
init = low;
id = low:-1:1;

% Calculate dU_in/dx1
dU_l = zeros(length(U_in(id)),1);
for i=id
    j = low - i + 1;
    if i ~= 1
        dU_l(j) = (U_in(i-1) - U_in(i))/(norm(rp(i-1,:) - rp(i,:))/c);
    else
        dU_l(j) = (U_in(i) - U_in(i+1))/(norm(rp(i,:) - rp(i+1,:))/c);
    end
end
d2U_l1 = (dU_l(2) - dU_l(1))/(norm(rp(init-1,:) - rp(init,:))/c);

% Variables initialization: A, lambda, H, l, etc
A(init) = 7.052; % stable
lambda(init) = lambda_fun(A(init));
H(init) = H_fun(A(init));
l(init) = l_fun(A(init));
F(init) = -.0652857*d2U_l1/dU_l(1)^2;
%F(init) = 2*(l(init) - (2 + H(init))*lambda(init))/U_in(init);
Z(init) = lambda(init)/dU_l(1);

% For transition check
xt = zeros(length(dU_l),1);

% Iteration to calculate handle variables for all position
for i = 1:length(dU_l)-1
    j = low-i+1;
    Z(j-1) = Z(j) + F(j)*ds(j)/c;
    lambda(j-1) = dU_l(i+1)*Z(j-1);
    if lambda(j-1) <=0
        G = fzero(@(L) lambda_fun(L) - lambda(j-1), -8);
    elseif lambda(j-1) > 0
        G = fzero(@(L) lambda_fun(L) - lambda(j-1), 8);
    end
    A(j-1) = G;

    % Calculate other variables
    H(j-1) = H_fun(G);
    l(j-1) = l_fun(G);
    F(j-1) = 2*(l(j-1) - (2 + H(j-1))*lambda(j-1))/U_in(j-1);

    clear G

    % Calculate xt
    xt(i+1) = sum(ds(j-1:init-1));

end

% Construct the BL result variables
delta(id) = sqrt(abs(A(id)*kmu./dU_l));
deltas(id) = delta(id).*(3/10 - A(id)./120);
thetas(id) = delta(id).*(37/315 - A(id)./945 - (A(id).^2)./9072);
tw(id) = (mu*U_in(id)./delta(id)).*(2 + A(id)./6);
cf(id) = 2*kmu*l_fun(A(id))./(U_in(id).*thetas(id));

% Checking the transition point
RHS = abs(U_in(id)).*thetas(id)./kmu;
k = abs(U_in(id)).*xt./kmu;
LHS = 1.174*(1+22400./k).*k.^0.46; % Cebeci and Smith (1974)
trans_l = abs(RHS - LHS);

temp = find(trans_l <= 1e-05);
if isempty(temp)
    trans_id(2) = nan;
else
    trans_id(2) = temp;
end

clear RHS LHS xt id temp

end