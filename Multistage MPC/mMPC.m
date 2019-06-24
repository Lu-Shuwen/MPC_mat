function [u0,x1] = mMPC(Np,x0,w,...
                  AA,BBu,BBw,QQ,RR,FF,ff,GG,gg,...
                  d,omega)
% Scenario tree based multi-stage (MPC) optimization at each sampling time
% Shuwen Lu - June 24, 2019

% # of candidate disturbances
N = length(d);
% Candidate disturbance vector over a planning horizon
dd = kron(d,ones(Np,1));

% Decision variables
u = sdpvar(Np,N);
% State varibles
x = sdpvar(4*Np,N);
% Scenario cost
J = sdpvar(3,1);
% Constraints
Cons = [];
for j = 1:N
    x(:,j) = AA*x0 + BBu*u(:,j) + BBw*dd(:,j);
    J(j) = x(:,j)'*QQ*x(:,j) + u(:,j)'*RR*u(:,j);
    Cons = [Cons, FF * x(:,j) <= ff, GG * u(:,j) <= gg];
end
Cons = [Cons, u(1,1) == u(1,2), u(1,1) == u(1,3)];
% Objective function
Obj = omega'*J;

% Solve the problem
sol = optimize(Cons,Obj,'');

% Analyze error flags
if sol.problem == 0
% Extract value
    u = value(u);
    x = value(x);
else
    disp('One more step close to success');
    sol.info
    yalmiperror(sol.problem)
end 

% Update control variable
u0 = u(1,1);
% Update states: real state with real disturbance w
x1 = AA(1:4,:)*x0 + BBu(1:4,1)*u0 + BBw(1:4,1)*w;