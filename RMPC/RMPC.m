function [h,M,u,x] = RMPC(H,x0,w,...
                  PhiPhiTExp,wPhiTExp,...
                  AA,BBu,BBw,QQ,RR,FF,ff,GG,gg,...
                  gammaMin,gammaMax)
% Robust Model Predictive Control (RMPC) or in reference [33]
% Second step: Solve RMPC problem with robust constraints solved in equ (47)
% H is the control horizon

% Robust counterpart derivation (6800 slides)
h  = sdpvar(H,1);
M  = sdpvar(H,H);
Lambda = sdpvar(4*H,6*H);           % dual variable in matrix form
for i = 1:H
    M(i,i:end) = 0;
end

% lifted uncertainty w_hat in equ (39) satisfying max(c'* w_hat) <= b
c = [FF * [BBu*M  BBw]
     GG * [M zeros(H)]];
b = [ff - FF * (AA*x0 + BBu*h)
     gg - GG * h];
% General polyhedral uncertainty set
%  |phi(w)| <= 1 and gammaMin <= w <= gammaMax
W = blkdiag(kron([1;-1],eye(H)),kron([1;-1],eye(H)));
v = [ones(2*H,1);gammaMax;-gammaMin];

% Input constraints equ (13) and state constraints equ (49) in duality
Cons = [Lambda'*v <= b, Lambda >= 0, Lambda'*W == c];
        
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
