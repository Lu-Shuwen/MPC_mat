% function [h,M,u,x] = DRMPC(H,x0,w,...
%                   PhiPhiTExp,wPhiTExp,Q,...
%                   AA,BBu,BBw,QQ,RR,FF,ff,GG,gg,...
%                   alpha,SV,theta_tilde,w_hat)
% Data-Dirven Model Predictive Control (DRMPC)
% Second step: Solve DRMPC problem with robust constraints solved in equ (47)
% H is the control horizon

clc;
close all; 
clear;
dbstop if error
% Nt is No. of samples for each row block, Nt = N
Nt = 300;
Nc = 59;
% Control horizon
H = 5;
nu = 0.05;
load('disturbancedata.mat','ww');
[~,PhiPhiTExp,wPhiTExp,Q,L] = Est(Nt,Nc,H,ww);
[AA,BBu,BBw,QQ,RR,FF,ff,GG,gg] = model(H);
[alpha,SV,theta_tilde,W_SV_hat] = uncertDRMPC(Nt,Nc,H,ww,Q,L,nu);
x0 = [0.2 1 -0.1 0.1]';
w = ww(2000);

len = length(SV);
% Robust counterpart derivation (6800 slides)
h  = sdpvar(H,1);
M  = sdpvar(H,H);

% for i = 1:H
%     M(i,i:end) = 0;
% end

lambda = sdpvar(2*H,1);
Heta  = sdpvar(6*H,1);
Lambda = sdpvar(2*H,len,6*H);           % dual variable in matrix form
Mu = sdpvar(2*H,len,6*H);

% Uncertainty parameter u is lifted uncertainty w_hat in equ (39)
W = [kron([1; -1],eye(H)),zeros(2*H,H)];
V = ones(2*H,1);
 
% Input constraints equ (13) in ref [31] and state constraints equ (49)
% max(c* w_hat) <= b

c = [FF * [BBu*M BBw];
     GG * [M zeros(H)] ];
b = [ff - FF * (AA*x0 + BBu*h)
     gg - GG * h];

Cons = [];
for k = 1:length(b)
    Cons = [Cons,trace((Mu(:,:,k)-Lambda(:,:,k))'*Q*W_SV_hat)+Heta(k)*theta_tilde...
                + V'* lambda       <= b(k);
        sum(Q'*(Lambda(:,:,k)-Mu(:,:,k)),2) + W'*lambda  == c(k,:)',...
        Lambda(:,:,k)+ Mu(:,:,k) == Heta(k)*ones(2*H,1)*alpha(SV)'];
end
Cons = [Cons, lambda(:) >= 0, Lambda(:) >= 0, Mu(:) >= 0, Heta(:) >= 0];
    
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
