function [Z,z] = uncertDPGMMMPC(N,H,W)
% Robust Model Predictive Control (RMPC) or in reference [15]
% Fisrt step: bounding the uncertainty or disturbance data (47)
% Hyper-rectangular uncertainty boundary of Wt(D) in equ (47)
% N is the no. of multi-samples of disturbance
% H is the control horizon
% The uncertainty data satisifies Z * w <=z

% Dirichlet Process Gaussian Mixture Model- Data Driven Uncertainty set
addpath(genpath('\Users\lushuwen\Downloads\research\MPC\DPGMMMPC\DPMM'));
% Variational Dirichlet process Gaussian mixture model by Blei and Jordan
T = 10;                       % with T=10
opts = mkopts_bj(T);
result = vdpgm(W, opts);

% Extract posterior parameters (Notations are the same as in the paper)
tau = result.hp_posterior.gamma(1,:);
nu = result.hp_posterior.gamma(2,:);

% Find m, where m  is the no. of components of Gaussian mixture
gamma(1) = tau(1)/(tau(1)+nu(1));
i = 1;
gamma_star = 0.02;            % threshold value
while gamma(i) >= gamma_star && i < T
    i = i + 1;
    gamma(i) = tau(i)/(tau(i)+nu(i));
    for j = 1:i-1
        gamma(i) = gamma(i) * nu(j)/(tau(j)+nu(j));
    end
end
m = length(gamma) - 1;      % m = 1

% Extract other parameters
lambda = result.hp_posterior.xi(1:m);
omega = result.hp_posterior.eta(1:m);
mu = result.hp_posterior.m(:,1:m);
% Probability of belonging to each Gaussian component
p = result.hp_posterior.q_of_z(:,1:m);  
% inv_psi is the inverse matrix of psi
inv_psi = result.hp_posterior.inv_B(:,:,1:m);

% Label kth data point belongs to which component
comp = zeros(N,1);       
for k = 1:N
    [~,comp_idx] = max(p(k,:));
    if p(k,comp_idx) > 0.67
        comp(k) = comp_idx;
    else
        comp(k) = 0;
    end
end

% s in eqn 7
s = sqrt((lambda+1)./(lambda.*(omega+1-H)));
% Gamma satisfies the probability guarantee
Gamma = [6.80 6.80 6.80];       % m = 1, only the first is used

% Find all constraints satisfied by w, say Z*w<=z
Z = zeros(2^H+2*H,H,m);
z = zeros(2^H+2*H,m);
% Find A, b such that A*w<=b in Gaussian mixture model
A = zeros(2^H,H,m);
b = zeros(2^H,m);
% Find upper bound ub and lower bound lb: l infinity norm
ub = zeros(H,m);
lb = zeros(H,m);
matrix = create_matrix(H);
for i = 1:m
    ub(:,i) = max(W(:,comp==i),[],2);
    lb(:,i) = min(W(:,comp==i),[],2);
    in = sqrtm(inv_psi(:,:,i))/s(i);
    A(:,:,i) = matrix * in;
    b(:,i) = Gamma(i) + matrix * in * mu;
    % Combine A*w<=b and lb<=w<=ub into Z and z
    Z(:,:,i) = [A;kron([1;-1],eye(H))];
    z(:,i) = [b;ub;-lb];
end