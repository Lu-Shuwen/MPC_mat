function [AA,BBu,BBw,QQ,RR,FF,ff,GG,gg] = model(H)
% Shuwen Lu - 03/14/2019

% This function builds the model of the two-mass-spring system
% and also the contraints matrices


% clc;
% close all; 
% clear;
% 
% dbstop if error
% 
% % Nt is No. of samples for each row block, Nt = N
% Nt = 300;
% Nc = 59;
% % Control horizon
% H = 5;
% nu = 0.05;
% 
% load('disturbancedata.mat','ww');


% Parameters
K  = 1;
m1 = 0.5;
m2 = 2;

% Model of the system
A  = [1        0             0.1    0;
      0        1             0      0.1;
     -K/m1     0.1*K/m1      1      0;
      K/m2    -0.1*K/m2      0      1];
Bu = [0        0             0.1/m1 0]';
Bw = [1.0      0.5           0.3    0.4]';

[nrBu,ncBu] = size(Bu);         % No. of rows and columns of Bu
[nrBw,ncBw] = size(Bw);         % No. of rows and columns of Bw

% System matrices in equ (3)
AA = [];                        % Initialize
BBu = zeros(H*nrBu,H*ncBu);
BBw = zeros(H*nrBw,H*ncBw);

for i = 1:H
    AA = [AA; A^i];
    BBu = BBu + kron(diag(ones(1,H-(i-1)),-(i-1)),A^(i-1)*Bu);
    BBw = BBw + kron(diag(ones(1,H-(i-1)),-(i-1)),A^(i-1)*Bw);
end

% Coefficient matrices in objective function
Q  = 5*eye(4);
R  = 1;
Qf = eye(4);

QQ = blkdiag(kron(eye(H-1),Q),Qf);
RR = kron(eye(H),R);

% State constraints Fx<=f and robust input constraints
F = [0  0  1  0;
     0  0  0  1;
     0  0 -1  0;
     0  0  0 -1];
f = 0.38*ones(4,1);
FF = kron(eye(H),F);
ff = kron(ones(H,1),f);

G = [ 1;
     -1];
g = 1.6*ones(2,1);
GG = kron(eye(H),G);
gg = kron(ones(H,1),g);