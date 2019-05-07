function [h,M,u,x] = CEMPC(H,x0,w,...
                  PhiPhiTExp,wPhiTExp,...
                  AA,BBu,BBw,QQ,RR,FF,ff,GG,gg)
% Certainty-Equivalent Model Predictive Control (CEMPC)
% Second step: Solve RMPC problem with robust constraints solved in equ (47)
% H is the control horizon

h  = sdpvar(H,1);
M  = sdpvar(H,H);
xx = sdpvar(4*H,1);
uu = h;

for i = 1:H
    M(i,i:end) = 0;
end

% Input constraints equ (13) and state constraints equ (49) in duality
Cons = [xx == AA * x0 + BBu * uu,...
        FF * xx <= ff,...
        GG * uu <= gg];
        
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

% update control variable equ (9)
u = h(1);
% Update states: equ (14) and (12)
x = AA(1:4,:)*x0 + BBu(1:4,1)*u + BBw(1:4,1)*w;