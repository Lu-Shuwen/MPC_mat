% Implement the paper https://doi.org/10.1016/j.jprocont.2018.12.013
% Robust Model Predictive Control (RMPC) or in reference [15]
clc;
close all; 
clear;

dbstop if error

% Nt is No. of samples for each row block, Nt = N
N = 311;
% Control horizon
H = 5;

load('disturbancedata.mat','ww');

[~,PhiPhiTExp,wPhiTExp,W] = Est(N,H,ww(1:H*N));
[AA,BBu,BBw,QQ,RR,FF,ff,GG,gg] = model(H);
[gammaMin,gammaMax] = uncertRMPC(N,H,W);

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
    [h,M,u,x] = RMPC(H,x0,ww(2000+k-1),...          % Initial conditions
                  PhiPhiTExp,wPhiTExp,...           % Estimate expectations
                  AA,BBu,BBw,QQ,RR,FF,ff,GG,gg,...  % System model
                  gammaMin,gammaMax);               % uncertainty
    % Only the first step is implemented
    uu(k) = u;
    x0 = x(1:4);
    xx(:,k) = x0';
end

% Save data
save('RMPC.mat','uu','xx');
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