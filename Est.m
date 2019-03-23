function [PhiExp,PhiPhiTExp,wPhiTExp,Q] = Est(N,H)
% Estimate the expectation of phi(w)phi(w)^T and wphi(w)^T from data

load('disturbancedata.mat','ww')

% Initialize

PhiPhiTSum = zeros(H,H);
wPhiTSum = zeros(H,H);

len = length(ww)/H;
W = zeros(H,len);

for i = 1:len
    W(:,i) = ww((i-1)*H+1:(i-1)*H+H)';
end

% Check whether the expectation of phi(w) is zero 
PhiSum = sum(tanh(W),2);
% Estimate the two expectations in equ (16)
for i = 1:length(W(1,:))
    
    PhiPhiTSum = PhiPhiTSum + tanh(W(:,i))*tanh(W(:,i))';
    wPhiTSum = wPhiTSum + W(:,i)*tanh(W(:,i))';
end
PhiExp = PhiSum/floor(len);
PhiPhiTExp = PhiPhiTSum/floor(len);
wPhiTExp = wPhiTSum/len;

% Calculate weighting matrix for DRMPC
temp = zeros(H,H);
for i = 1:N
    temp = temp+W(:,i)*W(:,i)';
end
sigma = 1/(N-1)*(temp - 1/N*sum(W(:,1:N),2)*sum(W(:,1:N),2)');
Q = sigma^(-0.5);