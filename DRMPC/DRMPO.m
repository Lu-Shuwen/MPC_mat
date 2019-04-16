% Implement the paper https://doi.org/10.1016/j.jprocont.2018.12.013
% Robust Model Predictive Control (RMPC) or in reference [15]
clc;
close all; 
clear;

dbstop if error

% Nt is the no. of training dataets of disturbance vector
Nt = 300;
% Nc is the no. of calibration datasets of disturbance vector
Nc = 59;
% Control horizon
H = 5;
% Probability guarantee
nu = 0.05;

load('disturbancedata.mat','ww');

[~,PhiPhiTExp,wPhiTExp,Q,L,W_hat] = Est(Nt,Nc,H,ww);
[AA,BBu,BBw,QQ,RR,FF,ff,GG,gg] = model(H);
[alpha,SV,theta_tilde,W_SV_hat] = uncertDRMPC(Nt,Nc,Q,L,nu,W_hat);

x0 = [0.2 1 -0.1 0.1]';

% Sampling time
Ts = 0.1;
% Time horizon
t = 10;
% Time step
T = t/Ts;

uu = zeros(1,T);
xx = zeros(4,T);

for k = 1:T
    [~,~,u,x] = DRMPC(H,x0,ww(2000+k-1),...
                  PhiPhiTExp,wPhiTExp,Q,...
                  AA,BBu,BBw,QQ,RR,FF,ff,GG,gg,...
                  alpha,SV,theta_tilde,W_SV_hat);
    % Only the first step is implemented
    uu(k) = u;
    x0 = x;
    xx(:,k) = x0;
end

% Save data
save('DRMPC.mat','uu','xx');
% Plot states and control variables
time = Ts:Ts:t;
figure(1)
plot(time,xx(1,:));
title('x(1) vs time');
xlabel('time (s)');
ylabel('x(1)');

figure(2)
plot(time,xx(2,:));
title('x(2) vs time');
xlabel('time (s)');
ylabel('x(2)');

figure(3)
plot(time,xx(3,:));
title('x(3) vs time');
xlabel('time (s)');
ylabel('x(3)');

figure(4)
plot(time,xx(4,:));
title('x(4) vs time');
xlabel('time (s)');
ylabel('x(4)');

figure(5)
plot(time,uu);
title('u vs time');
xlabel('time (s)');
ylabel('u');