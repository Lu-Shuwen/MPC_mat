function [h,M,uu,xx] = RMPC(H,x0,...
                  PhiPhiTExp,wPhiTExp,...
                  AA,BBu,BBw,QQ,RR,FF,ff,GG,gg,...
                  gammaMin,gammaMax)
% Robust Model Predictive Control (RMPC) or in reference [15]
% Second step: Solve RMPC problem with robust constraints solved in equ (47)
% N is the no. of multi-samples of disturbance
% H is the control horizon
% x0 is current states, the initial condition of next control horizon

h  = sdpvar(H,1);
M  = sdpvar(H,H);
uu = sdpvar(H,1);
w  = sdpvar(H,1);

% control variable equ (9)
for t = 1:H
    uu(t) = h(t);
    for j = 1:t-1
        uu(t) = uu(t)+M(t,j)*tanh(w(j));
    end
end

% Input constraints equ (13)
Cons = [GG*uu <= gg];

% State constraints equ (49)
for t = 1:H
    Cons = [Cons, gammaMin <= w, w <= gammaMax,...
    FF(4*t-3:4*t,:) * (AA*x0 + BBu*M*tanh(w) + BBu*h + BBw*w) <= ff(4*t-3:4*t)];
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
    uu = value(uu);
    w = value(w);
else
    disp('One more step close to success');
    sol.info
    yalmiperror(sol.problem)
end 

% Update states: equ (14) and (12)
xx = AA*x0 + BBu*uu + BBw*w;