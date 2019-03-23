function [h,M,uu,xx] = RMPC(H,x0,...
                  PhiPhiTExp,wPhiTExp,...
                  AA,BBu,BBw,QQ,RR,FF,ff,GG,gg,...
                  gammaMin,gammaMax)
% Robust Model Predictive Control (RMPC) or in reference [15]
% Second step: Solve RMPC problem with robust constraints solved in equ (47)
% H is the control horizon

h  = sdpvar(H,1);
M  = sdpvar(H,H);
% uu = sdpvar(H,1);
% w  = sdpvar(H,1);

% Simulate for all w belongs to W
% Sample size to simulate all w's
S = 100;
w = gammaMin+(gammaMax - gammaMin).*rand(H,S);

% Input constraints equ (13)
Cons = [];

% State constraints equ (49)
for i = 1:S
    Cons = [Cons, GG*(h + M*tanh(w(:,i))) <= gg,
    FF * (AA*x0 + BBu*M*tanh(w(:,i)) + BBu*h + BBw*w(:,i)) <= ff];
end
        
        
% Objective function J(M,h) equ (16)
Obj = (h'*(RR+BBu'*QQ*BBu)*h + ...
          trace((RR+BBu'*QQ*BBu)*M*PhiPhiTExp*M') + ...
          2*h'*BBu'*QQ*AA*x0 + 2*trace(M'*BBu'*QQ*BBw*wPhiTExp))

% Solve the problem
sol = optimize(Cons,Obj,'');

% Analyze error flags
if sol.problem == 0
% Extract value
    h = value(h);
    M = value(M);
else
    disp('One more step close to success');
    sol.info
    yalmiperror(sol.problem)
end 
% control variable equ (9)
uu = h;
% Update states: equ (14) and (12)
xx = AA*x0 + BBu*uu;
