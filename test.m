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

% I change the name of variables to be u1 x1 to distinguish from model test 
load('ctrlperformance_RMPC_038_H5.mat','u1','x1');
load('disturbancedata.mat','ww');

% x2 stores for the corresponding states under control u1
x2 = zeros(4,100);
x0 = [0.2 1 -0.1 0.1]';

for k = 1:100
    xx = A*x0 + Bu*u1(k,:)+Bw*ww(k);
    x2(:,k) = xx;
    x0 = xx;
end
save('model_test.mat','x2');

time = 0.1:0.1:10;
figure(1)
plot(time,x2(1,:));
title('x(1) vs time');
xlabel('time (s)');
ylabel('x(1)');

figure(2)
plot(time,x2(2,:));
title('x(2) vs time');
xlabel('time (s)');
ylabel('x(2)');

figure(3)
plot(time,x2(3,:));
title('x(3) vs time');
xlabel('time (s)');
ylabel('x(3)');

figure(4)
plot(time,x2(4,:));
title('x(4) vs time');
xlabel('time (s)');
ylabel('x(4)');