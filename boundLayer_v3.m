% Subroutine for boundary layer calculation
% using Approximation method: Karman-Pohlhausen
% NOTE: *BL calculation only for laminar region
%       *transition checked using Cebeci and Smith (1974) method
%           which is an improvement of Michel's method
%       *turbulent region is neglected

function [delta, deltas, thetas, tw, cf, trans_id, ind_stag] = boundLayer_v3...
    (U_in, xb, yb, xp, yp, c, rho, U_inf, mu)
%% Initialize input variables
if size(U_in, 2) ~= 1
    U_in = U_in';
end
kmu = mu/rho;
Re = U_inf*c/kmu;

% Obtain stagnation point
[~, ind_stag] = min(abs(U_in));
%disp(ind_stag);

up = ind_stag+1;
low = ind_stag-1;

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
x0 = 0.01;
factor = 1000;

%% Solution for UPPER airfoil
init = up;
id = up:length(U_in);

% Calculate dU_in/dx1
dU_u = zeros(length(U_in(id)),1);
for i=id
    j = i - up + 1;
    if i ~= length(U_in)
        dU_u(j) = (U_in(i+1) - U_in(i))/((xp(i+1) - xp(i))/c);
    else
        dU_u(j) = (U_in(i) - U_in(i-1))/((xp(i) - xp(i-1))/c);
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
    
    lambda = (tetau^2)*dU_u(i)/(factor*kmu);
    L(j) = lambda;
    G = fzero(@(x)(x*((37/315)-(x/945)-(x^2/9072))^2 - lambda),x0,options);
    if G < -12
        G = -12;
    elseif G > 12
        G = 12;
    end
    
    A(j) = G;

    clear G lambda
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
cf(id) = kmu*(2+A(id)/6)./(U_in(id).*delta(id));

% Checking the transition point
RHS = abs(U_in(id)).*thetas(id)./kmu;
k = abs(U_in(id)).*xt./kmu;
LHS = 1.174*(1+22400./k).*k.^0.46; % Cebeci and Smith (1974)
trans_u = abs(RHS - LHS);

temp = find(trans_u <= 1e-05);
if isempty(temp)
    trans_id(1) = length(xp)*1000;
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
        dU_l(j) = (U_in(i-1) - U_in(i))/((xp(i-1) - xp(i))/c);
    else
        dU_l(j) = (U_in(i) - U_in(i+1))/((xp(i) - xp(i+1))/c);
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
    
    lambda = (tetal^2)*dU_l(i)/(factor*kmu);
    L(j) = lambda;
    G = fzero(@(x)(x*((37/315)-(x/945)-(x^2/9072))^2 - lambda),x0,options);    
    if G < -12
        G = -12;
    elseif G > 12
        G = 12;
    end
    
    A(j) = G;
    
    clear G lambda
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
cf(id) = kmu*(2+A(id)/6)./(U_in(id).*delta(id));

% Checking the transition point
RHS = abs(U_in(id)).*thetas(id)./kmu;
k = abs(U_in(id)).*xt./kmu;
LHS = 1.174*(1+22400./k).*k.^0.46; % Cebeci and Smith (1974)
trans_l = abs(RHS - LHS);

temp = find(trans_l <= 1e-05);
if isempty(temp)
    trans_id(2) = length(xp)*1000;
else
    trans_id(2) = temp;
end

clear RHS LHS xt id temp

%% Properties calculation at stagnation point
i = ind_stag;
dU_stag = abs((U_in(i+1) - U_in(i))/((xp(i+1) - xp(i))/c));
thetas(i) = sqrt(0.075*kmu/dU_stag);
delta(i) = thetas(i)*315/37;
deltas(i) = delta(i)*3/10;
tw(i) = mu*U_in(i)*2/delta(i);
end