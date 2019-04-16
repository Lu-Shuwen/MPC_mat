function [PhiExp,PhiPhiTExp,wPhiTExp,Q,L,W_hat] = Est(Nt,Nc,H,ww)
% Estimate the expectation of phi(w)phi(w)^T and wphi(w)^T from data


% Initialize
N = Nt+Nc;
PhiPhiTSum = zeros(H,H);
wPhiTSum = zeros(H,H);

W = zeros(H,N);
for i = 1:N
    W(:,i) = ww(i:i+H-1)';
end

% % Saturating of data
phi = (1 - exp(-200 * W)) ./ (1 + exp(-200 * W));
% lifted lifted uncertainty W_hat in eqn (39)
W_hat = [phi;W];
% W_hat = [tanh(W);W];

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

% Covariance Matrix
temp = zeros(2*H,2*H);
for i = 1:Nt
    temp = temp+W_hat(:,i)*W_hat(:,i)';
end
sigma = 1/(Nt-1)*(temp - 1/Nt*sum(W_hat(:,1:Nt),2)*sum(W_hat(:,1:Nt),2)');
% Calculate weighting matrix for DRMPC
Q = sigma^(-0.5);

% L is sufficiently large, estimate lower bound by [25] equ (31)
l = max(Q*W_hat,[],2)-min(Q*W_hat,[],2);
% L = sum(l), add 2*H to compensate for l(k) strictly larger
L = sum(l) + 2*H;