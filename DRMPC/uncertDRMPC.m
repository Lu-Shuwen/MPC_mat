function [alpha,SV,theta_tilde,W_SV_hat] = uncertDRMPC(Nt,Nc,Q,L,nu,W_hat)
% Data Driven Model Predictive Control (DRMPC) or in reference [25]
% Fisrt step: bounding the uncertainty or disturbance data eqn (21)
% Find uncertainty set W(mu,D) in equ (26) (28)
% Nt is the no. of training dataets of disturbance vector
% Nc is the no. of calibration datasets of disturbance vector
% alpha's are Lagrange multipliers needed for spporting vectors
% to establish the enclosing circle in the feature space
% SV is the index set of supporting vectors
% BSV is the index set of boundary support vectors


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
% supporting vectors
W_SV_hat = W_hat(:,SV);