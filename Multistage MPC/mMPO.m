% Scenario tree based multi-stage (MPC) optimization with robust horizon 1
% Shuwen Lu - June 24, 2019

clc;
close all; 
clear;

dbstop if error

% Prediction horizon
Np = 5;

% Disturbance w is the nomally distributed uncertain parameter
% w ~ N(mu,sigma^2), where mu = 0 and sigma = 1/200
mu = 0;
sigma = 1/200;
% 3 candidate disturbances truncated at 2*sigma = 0.01
d = [-2*sigma 0 2*sigma] + mu;

% Probability weight
omega = [0.1587 0.6826 0.1587]';

load('disturbancedata.mat','ww');

% Get parameters for the model
[AA,BBu,BBw,QQ,RR,FF,ff,GG,gg] = model(Np);

% Initial conditions
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
    [u0,x1] = mMPC(Np,x0,ww(2000+k-1),...           % Initial conditions
            AA,BBu,BBw,QQ,RR,FF,ff,GG,gg,...        % System model
            d,omega);      % Candidate uncertainty and probability weight
      
    % Only the first step is implemented
    uu(k) = u0;
    x0 = x1;
    xx(:,k) = x0;
end

% Save data
save('mMPC.mat','uu','xx');
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