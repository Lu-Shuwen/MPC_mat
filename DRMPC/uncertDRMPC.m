function [alpha,SV,theta_tilde,W_SV_hat] = uncertDRMPC(Nt,Nc,H,ww,Q,L,nu)
% Data Driven Model Predictive Control (DRMPC) or in reference [25]
% Fisrt step: bounding the uncertainty or disturbance data (21)
% N is the no. of multi-samples of disturbance
% H is the control horizon
% alpha's are Lagrange multipliers needed for spporting vectors
% to establish the enclosing circle in the feature space
% SV is the index set of supporting vectors
% BSV is the index set of boundary support vectors
% Nt is different for each t-th row block of constraints in (49)
% Since epsilon_t and beta_t are the same for each t,
% they share the same gammaMin and gammaMax for simplicity, Nt = N

% clc;
% close all; 
% clear;
% 
% dbstop if error
% 
% % Nt is No. of samples for each row block, Nt = N
% Nt = 300;             % training dataset
% Nc = 59;              % calibration dataset
% 
% % Control horizon
% H = 5;
% nu = 0.05;
% load('disturbancedata.mat','ww');
% [~,~,~,Q,L] = Est(Nt,Nc,H,ww);

% Uncertainty boundary of W(mu,D) in equ (26) (28)
N = Nt + Nc;
W = zeros(H,N);
for i = 1:N
    W(:,i) = ww((i-1)*H+1:(i-1)*H+H)';
                        % N columns meaning no. of scenarios
                        % H rows meaning each scenario has H horizons
end

W_hat = [tanh(W);W];

% Solve for Lagrange multipliers based on a dual quadratic program equ (21)
% Initialize for Kernel function
K = zeros(Nt,Nt);
for i = 1:Nt
    for j = 1:Nt
        K(i,j) = L - norm(Q*(W_hat(:,i)-W_hat(:,j)),1);
    end
end

% Lagrange multipliers
alpha = sdpvar(Nt,1);

% Define constraints 
Cons = [0 <= alpha(:) <= 1/(Nt*nu), sum(alpha) == 1];

% Define an objective
Obj = alpha'* K * alpha - alpha'*diag(K);

% Solve the problem
sol = optimize(Cons,Obj,'');

% Analyze error flags
if sol.problem == 0
% Extract value
    alpha = value(alpha);
else
    disp('One more step close to success');
    sol.info
    yalmiperror(sol.problem)
end 

% index set of supporting vectors and boundary support vectors
SV = find(alpha > 1e-7);
BSV = find(alpha > 1e-7 & 1/(Nt*nu)-alpha > 1e-7 );

% theta in equ (27)
theta = 0;
for i = 1:length(SV)
        theta = theta + alpha(SV(i))*norm(Q*(W_hat(:,BSV(1))-W_hat(:,SV(i))),1);
end

% calibration of theta in (31)
f_w = zeros(1,Nc);
for k = 1:Nc
    for i = 1:length(SV)
        f_w(k) = f_w(k) + alpha(SV(i))*norm(Q*(W_hat(:,Nt+k)-W_hat(:,SV(i))),1);
    end
end
theta_tilde = max(f_w);
W_SV_hat = [tanh(W(:,SV));W(:,SV)];
