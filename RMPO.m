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

[~,PhiPhiTExp,wPhiTExp,~] = Est(N,H);
[AA,BBu,BBw,QQ,RR,FF,ff,GG,gg] = model(H);
[gammaMin,gammaMax] = uncertRMPC(N,H);

% Sampling time
Ts = 0.1;
% Time horizon
t = 10;
% Time step
T = t/Ts;

x0 = [0.2 1 -0.1 0.1]';
u = zeros(1,T);
x = zeros(4,T);

for k = 1:T
    [h,M,uu,xx] = RMPC(H,x0,...                     % Initial conditions
                  PhiPhiTExp,wPhiTExp,...           % Estimate expectations
                  AA,BBu,BBw,QQ,RR,FF,ff,GG,gg,...  % System model
                  gammaMin,gammaMax);               % uncertainty
    % Only the first step is implemented
    u(k) = uu(1);
    x0 = xx(1:4,1);
    x(:,k) = x0';
end

% Plot states and control variables
time = Ts:Ts:t;
figure(1)
plot(time,x(1,:));
title('x(1) vs time');
xlabel('time (s)');
ylabel('x(1)');

figure(2)
plot(time,x(2,:));
title('x(2) vs time');
xlabel('time (s)');
ylabel('x(2)');

figure(3)
plot(time,x(3,:));
title('x(3) vs time');
xlabel('time (s)');
ylabel('x(3)');

figure(4)
plot(time,x(4,:));
title('x(4) vs time');
xlabel('time (s)');
ylabel('x(4)');

figure(5)
plot(time,u);
title('u vs time');
xlabel('time (s)');
ylabel('u');