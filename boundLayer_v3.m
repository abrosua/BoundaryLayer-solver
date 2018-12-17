% Subroutine for boundary layer calculation
% using Approximation method: Karman-Pohlhausen
% NOTE: *BL calculation only for laminar region
%       *transition checked using Cebeci and Smith (1974) method
%           which is an improvement of Michel's method
%       *turbulent region is neglected

function [delta, deltas, thetas, tw, cf, trans_id, sp] = boundLayer_v3...
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
K = zeros(length(U_in),1);
L = zeros(length(U_in),1);
delta = zeros(length(U_in),1);
deltas = zeros(length(U_in),1);
thetas = zeros(length(U_in),1);
tw = zeros(length(U_in),1);
cf = zeros(length(U_in),1);
trans_id = zeros(2,1);

options = optimset('Display','off');

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

% For transition check
xt = zeros(length(dU_u),1);

tetau = 0;
% Iteration to calculate handle variables for all position
for i = 1:length(dU_u)
    j = i+up-1;
    
    tetau = tetau + sqrt((0.47*kmu/U_in(j)^6)*(U_in(j)^5 + U_in(j-1)^5)*...
        (abs((xp(j)-xp(j-1))))/2);
    K(j) = tetau;
    
    lambda = (tetau^2)*dU_u(i)/kmu;
    L(j) = lambda;
    G = fsolve(@(x)(x*((37/315)-(x/945)-(x^2/9072))^2 - lambda),0,options);
    
    A(j) = G;

    clear G
end

% Calculate xt
for i = 2:length(dU_u)
    j = i+up-1;
    for m = init:j
        xt(i) = xt(i) + abs(xp(m) - xp(m-1));
    end
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

% For transition check
xt = zeros(length(dU_l),1);

tetal = 0;
% Iteration to calculate handle variables for all position
for i = 1:length(dU_l)
    j = low-i+1;

    tetal = tetal + sqrt((0.47*kmu/U_in(j)^6)*(U_in(j)^5 + U_in(j+1)^5)*...
        (abs((xp(j)-xp(j+1))))/2);
    K(j) = tetal;
    
    lambda = (tetal^2)*dU_l(i)/kmu;
    L(j) = lambda;
    G = fsolve(@(x)(x*((37/315)-(x/945)-(x^2/9072))^2 - lambda),0,options);    
    
    A(j) = G;
    
    clear G
end

% Calculate xt
for i = 2:length(dU_l)
    j = low-i+1;
    for m = init:-1:j
        xt(i) = xt(i) + abs(xp(m) - xp(m+1));
    end
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