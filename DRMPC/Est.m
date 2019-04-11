function [PhiExp,PhiPhiTExp,wPhiTExp,Q,L] = Est(Nt,Nc,H,ww)
% Estimate the expectation of phi(w)phi(w)^T and wphi(w)^T from data

% clc;
% clear;
% % Nt is No. of samples for each row block, Nt = N
% Nt = 300;
% Nc = 59;
% % Control horizon
% H = 5;
% 
% load('disturbancedata.mat','ww');

% Initialize
N = Nt+Nc;
PhiPhiTSum = zeros(H,H);
wPhiTSum = zeros(H,H);

W = zeros(H,N);

for i = 1:N
    W(:,i) = ww((i-1)*H+1:(i-1)*H+H)';
end

% Check whether the expectation of phi(w) is zero 
PhiSum = sum(tanh(W),2);
% Estimate the two expectations in equ (16)
for i = 1:N
    
    PhiPhiTSum = PhiPhiTSum + tanh(W(:,i))*tanh(W(:,i))';
    wPhiTSum = wPhiTSum + W(:,i)*tanh(W(:,i))';
end
PhiExp = PhiSum/N;
PhiPhiTExp = PhiPhiTSum/N;
wPhiTExp = wPhiTSum/N;

W_hat = [tanh(W(:,1:Nt));W(:,1:Nt)];
% Calculate weighting matrix for DRMPC
temp = zeros(2*H,2*H);
for i = 1:Nt
    temp = temp+W_hat(:,i)*W_hat(:,i)';
end
sigma = 1/(Nt-1)*(temp - 1/Nt*sum(W_hat(:,1:Nt),2)*sum(W_hat(:,1:Nt),2)');
Q = sigma^(-0.5);

% L is sufficiently large, estimate lower bound by [25] equ (31)
l = zeros(Nt,1);
for k = 1:2*H
    % l(k) needs to be larger than max(Q(:,k)'*W)-min(Q(:,k)'*W)
    l(k) = max(Q(:,k)'*W_hat(:,1:Nt))-min(Q(:,k)'*W_hat(:,1:Nt));
end
% L = sum(l), add H to compensate for l(k) strictly larger
L = sum(l) + 2*H;