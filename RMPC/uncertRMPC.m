function [gammaMin,gammaMax] = uncertRMPC(N,H,W)
% Robust Model Predictive Control (RMPC) or in reference [15]
% Fisrt step: bounding the uncertainty or disturbance data (47)
% Hyper-rectangular uncertainty boundary of Wt(D) in equ (47)
% N is the no. of multi-samples of disturbance
% H is the control horizon
% gammaMin, gammaMax are hyper-rectangular constraints of sampling data

% Solve for Hyper-rectangle captured scenarios based on equ (47)
gammaMin = sdpvar(H,1);
gammaMax = sdpvar(H,1);

% Define constraints 
Cons = [];
for i = 1:N
    Cons = [Cons, gammaMin <= W(:,i), gammaMax >= W(:,i)];
end

% Define an objective
Obj = norm(gammaMax-gammaMin,1);
 
% Solve the problem
sol = optimize(Cons,Obj,'');

% Analyze error flags
if sol.problem == 0
% Extract value
    gammaMin = value(gammaMin);
    gammaMax = value(gammaMax);
else
    disp('One more step close to success');
    sol.info
    yalmiperror(sol.problem)
end 