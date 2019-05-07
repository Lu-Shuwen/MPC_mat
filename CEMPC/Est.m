function [PhiExp,PhiPhiTExp,wPhiTExp,W_hat] = Est(N,H,ww)
% Estimate the expectation of phi(w)phi(w)^T and wphi(w)^T from data


% Initialize

PhiPhiTSum = zeros(H,H);
wPhiTSum = zeros(H,H);

W = zeros(H,N);

for i = 1:N
    W(:,i) = ww(i:i+H-1)';
end

% % Saturating of data
phi = (1 - exp(-200 * W)) ./ (1 + exp(-200 * W));
W_hat = [phi;W];
% Check whether the expectation of phi(w) is zero 
PhiSum = sum(phi,2);
% PhiSum = sum(tanh(W),2);

% Estimate the two expectations in equ (16)
for i = 1:N
    PhiPhiTSum = PhiPhiTSum + phi(:,i)*phi(:,i)';
    wPhiTSum = wPhiTSum + W(:,i)*phi(:,i)';
    
%     PhiPhiTSum = PhiPhiTSum + tanh(W(:,i))*tanh(W(:,i))';
%     wPhiTSum = wPhiTSum + W(:,i)*tanh(W(:,i))';
end
PhiExp = PhiSum/N;
PhiPhiTExp = PhiPhiTSum/N;
wPhiTExp = wPhiTSum/N;