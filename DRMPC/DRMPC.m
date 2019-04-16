function [h,M,u,x] = DRMPC(H,x0,w,...
                  PhiPhiTExp,wPhiTExp,Q,...
                  AA,BBu,BBw,QQ,RR,FF,ff,GG,gg,...
                  alpha,SV,theta_tilde,W_SV_hat)
% Data-Dirven Model Predictive Control (DRMPC)
% Second step: Solve DRMPC problem with robust constraints solved in equ (47)
% H is the control horizon


len = length(SV);

% control variables in eqn (16)
h  = sdpvar(H,1);
M  = sdpvar(H,H);
for i = 1:H
    M(i,i:end) = 0;
end

% Objective function J(M,h) equ (16)
Obj = (h'*(RR+BBu'*QQ*BBu)*h + ...
          trace((RR+BBu'*QQ*BBu)*M*PhiPhiTExp*M') + ...
          2*h'*BBu'*QQ*AA*x0 + 2*trace(M'*BBu'*QQ*BBw*wPhiTExp))
      
% dual variables in matrix form for state constraints
Heta  = sdpvar(4*H,1);
Lambda = sdpvar(2*H,len,4*H);           
Mu = sdpvar(2*H,len,4*H);
% dual variables in matrix form for input constraints
Z = sdpvar(2*H,2*H,'full');

Cons = [];

% lifted uncertainty w_hat in equ (39) satisfying max(c1* w_hat) <= b1
c1 = FF * [BBu*M BBw];
b1 = ff - FF * (AA*x0 + BBu*h);
% State constraints equ (49)
for k = 1:length(b1)
    Cons = [Cons,trace((Mu(:,:,k)-Lambda(:,:,k))'*Q*W_SV_hat)+Heta(k)*theta_tilde<= b1(k),...
        sum(Q'*(Lambda(:,:,k)-Mu(:,:,k)),2)  == -c1(k,:)',...
        Lambda(:,:,k)+ Mu(:,:,k) == Heta(k)*ones(2*H,1)*alpha(SV)'];
end

% max(c2* phi(W)) <= b2
c2 = GG *  M;
b2 = gg - GG * h;
% where phi(W) belongs to a polytopic uncertainty set, |phi(w)| <= 1
W = kron([1; -1],eye(H));
V = ones(2*H,1);
% Input constraints equ (13) robust counterpart derivation (6800 slides)
Cons = [Cons, Z'*V <= b2, Z'*W == c2,...
         Z(:) >= 0, Lambda(:) >= 0, Mu(:) >= 0, Heta(:) >= 0];
    

% Solve the problem
sol = optimize(Cons,Obj,'');

% Analyze error flags
if sol.problem == 0
% Extract value
    h = value(h);
%     M = value(M);
else
    disp('One more step close to success');
    sol.info
    yalmiperror(sol.problem)
end 

% update control variable equ (9)
u = h(1);
% Update states: equ (14) and (12)
x = AA(1:4,:)*x0 + BBu(1:4,1)*u + BBw(1:4,1)*w;
