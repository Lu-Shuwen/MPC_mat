function [gammaMin,gammaMax] = uncertRMPC(N,H)
% Robust Model Predictive Control (RMPC) or in reference [15]
% Fisrt step: bounding the uncertainty or disturbance data (47)
% N is the no. of multi-samples of disturbance
% H is the control horizon
% gammaMin, gammaMax are hyper-rectangular constraints of sampling data
% Nt is different for each t-th row block of constraints in (49)
% Since epsilon_t and beta_t are the same for each t,
% they share the same gammaMin and gammaMax for simplicity

load('disturbancedata.mat','ww');


    % Hyper-rectangular uncertainty boundary of Wt(D) in equ (47)
    W = zeros(H,N);         
    for i = 1:N
        W(:,i) = ww((i-1)*H+1:(i-1)*H+H)';
                            % N columns meaning no. of scenarios
                            % H rows meaning each scenario has H horizons
    end
    
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
    
    % Set some options for YALMIP and solver
%     opt = ...
%         sdpsettings('verbose',1,'solver','quadprog','quadprog.maxiter',100);
%     
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