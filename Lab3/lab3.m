close all;
clear all;
clc;
%% Part 1: Symbolic Linearization
% declare symbolic variables
syms y y_dot theta theta_dot t u M m l g real;

% define nonlinear systems symbolically
xdot1 = y_dot;
xdot2 = (-m*l*sin(theta)*theta_dot^2+m*g*sin(theta)*cos(theta)+u)/(M+m*sin(theta)^2);
xdot3 = theta_dot;
xdot4 = (-m*l*sin(theta)*cos(theta)*theta_dot^2+(M+m)*g*sin(theta)+u*cos(theta))/l/(M+m*sin(theta)^2);
x = [y;y_dot;theta;theta_dot];
xdot = [xdot1;xdot2;xdot3;xdot4];

% define structure parameters
parameters.M = 1.0731;
parameters.m = 0.23;
parameters.l = 0.3302;
parameters.g = 9.8;

% LINEARIZE NL system at given equilibria
% compute symbolic jacobians
A = jacobian(xdot,x);
B = jacobian(xdot,u);

% define equilibrium states
xbar = [0;0;0;0];

% evaluate jacobians at equilibria
A_evaluated = subs(A,x,xbar)
B_evaluated = subs(B,x,xbar)

% set numerical values for m,M,g,l
g = parameters.g;
l = parameters.l;
M = parameters.M;
m = parameters.m;

% evaluate symbolic matrices and convert them to double
A_double = double(eval(A_evaluated))
B_double = double(eval(B_evaluated))

%% Part 2: Controllability and Pole Assignment
% verify controllability
Q = ctrb(A_double, B_double)
Q_rank = rank(Q)

% pole placement_1
p_1 = [-1;-2;-3;-4];
K1 = -place(A_double, B_double, p_1) % convention u = Kx

% setting (ODE options, time range)
options = odeset('RelTol',1e-7,'AbsTol',1e-7);
Tspan = linspace(0,10,1e3);
x0 = [-0.5;0;-pi/4;0]; % initial condition

Acl_1 = A_double + B_double*K1; % convention u = Kx
sys_1 = @(t,x)Acl_1*x;
[t,x_1]=ode45(sys_1,Tspan,x0,options);

% compute u = Kx
u_1 = K1*x_1.';
u_1 = u_1.';

% pole placement_2
p_2 = [-1;-2;-3;-20];
K2 = -place(A_double, B_double, p_2) % convention u = Kx

Acl_2 = A_double + B_double*K2; % convention u = Kx
sys_2 = @(t,x)Acl_2*x;
[t,x_2]=ode45(sys_2,Tspan,x0,options);

% compute u = Kx
u_2 = K2*x_2.';
u_2 = u_2.';

figure('Name','Effect of Pole Assignment')
subplot(511)
plot(t,x_1(:,1),t,x_2(:,1))
legend('Pole 1','Pole 2')
ylabel('y')

subplot(512)
plot(t,x_1(:,2),t,x_2(:,2))
legend('Pole 1','Pole 2')
ylabel('y dot')

subplot(513)
plot(t,x_1(:,3),t,x_2(:,3))
legend('Pole 1','Pole 2')
ylabel('theta')

subplot(514)
plot(t,x_1(:,4),t,x_2(:,4))
legend('Pole 1','Pole 2')
ylabel('theta dot')

subplot(515)
plot(t,u_1,t,u_2)
legend('Pole 1','Pole 2')
ylabel('u')
xlabel('time')

%% Part3: LQR: changing q1
% setting (ODE time range)
Tspan = linspace(0,30,1e3);

R = 0.5; % fixed R
Q_1 = [0.1 0 0 0; 0 0 0 0; 0 0 5 0; 0 0 0 0]; % set1
Q_2 = [0.005 0 0 0; 0 0 0 0; 0 0 5 0; 0 0 0 0]; % set2

K1 = -lqr(A_double, B_double, Q_1, R);
K2 = -lqr(A_double, B_double, Q_2, R);

Acl_1 = A_double + B_double*K1; % convention u = Kx
sys_1 = @(t,x)Acl_1*x;
[t,x_1]=ode45(sys_1,Tspan,x0,options);

% compute u = Kx
u_1 = K1*x_1.';
u_1 = u_1.';

Acl_2 = A_double + B_double*K2; % convention u = Kx
sys_2 = @(t,x)Acl_2*x;
[t,x_2]=ode45(sys_2,Tspan,x0,options);

% compute u = Kx
u_2 = K2*x_2.';
u_2 = u_2.';

figure('Name','Effect of q1')
subplot(511)
plot(t,x_1(:,1),t,x_2(:,1))
legend('q1 = 0.1','q1 = 0.005')
ylabel('y')

subplot(512)
plot(t,x_1(:,2),t,x_2(:,2))
legend('q1 = 0.1','q1 = 0.005')
ylabel('y dot')

subplot(513)
plot(t,x_1(:,3),t,x_2(:,3))
legend('q1 = 0.1','q1 = 0.005')
ylabel('theta')

subplot(514)
plot(t,x_1(:,4),t,x_2(:,4))
legend('q1 = 0.1','q1 = 0.005')
ylabel('theta dot')

subplot(515)
plot(t,u_1,t,u_2)
legend('q1 = 0.1','q1 = 0.005')
ylabel('u')
xlabel('time')

%% LQR: changing q2
R = 0.5; % fixed R
Q_1 = [0.05 0 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 0]; % set1
Q_2 = [0.05 0 0 0; 0 0 0 0; 0 0 2000 0; 0 0 0 0]; % set2

K1 = -lqr(A_double, B_double, Q_1, R);
K2 = -lqr(A_double, B_double, Q_2, R);

Acl_1 = A_double + B_double*K1; % convention u = Kx
sys_1 = @(t,x)Acl_1*x;
[t,x_1]=ode45(sys_1,Tspan,x0,options);

% compute u = Kx
u_1 = K1*x_1.';
u_1 = u_1.';

Acl_2 = A_double + B_double*K2; % convention u = Kx
sys_2 = @(t,x)Acl_2*x;
[t,x_2]=ode45(sys_2,Tspan,x0,options);

% compute u = Kx
u_2 = K2*x_2.';
u_2 = u_2.';

figure('Name','Effect of q2')
subplot(511)
plot(t,x_1(:,1),t,x_2(:,1))
legend('q2 = 1','q2 = 2000')
ylabel('y')

subplot(512)
plot(t,x_1(:,2),t,x_2(:,2))
legend('q2 = 1','q2 = 2000')
ylabel('y dot')

subplot(513)
plot(t,x_1(:,3),t,x_2(:,3))
legend('q2 = 1','q2 = 2000')
ylabel('theta')

subplot(514)
plot(t,x_1(:,4),t,x_2(:,4))
legend('q2 = 1','q2 = 2000')
ylabel('theta dot')

subplot(515)
plot(t,u_1,t,u_2)
legend('q2 = 1','q2 = 2000')
ylabel('u')
xlabel('time')

%% LQR: changing R
R1 = 0.005; % R1
R2 = 10; % R2
Q = [0.05 0 0 0; 0 0 0 0; 0 0 5 0; 0 0 0 0]; % fixed Q

K1 = -lqr(A_double, B_double, Q, R1);
K2 = -lqr(A_double, B_double, Q, R2);

Acl_1 = A_double + B_double*K1; % convention u = Kx
sys_1 = @(t,x)Acl_1*x;
[t,x_1]=ode45(sys_1,Tspan,x0,options);

% compute u = Kx
u_1 = K1*x_1.';
u_1 = u_1.';

Acl_2 = A_double + B_double*K2; % convention u = Kx
sys_2 = @(t,x)Acl_2*x;
[t,x_2]=ode45(sys_2,Tspan,x0,options);

% compute u = Kx
u_2 = K2*x_2.';
u_2 = u_2.';

figure('Name','Effect of R')
subplot(511)
plot(t,x_1(:,1),t,x_2(:,1))
legend('R = 0.005','R = 10')
ylabel('y')

subplot(512)
plot(t,x_1(:,2),t,x_2(:,2))
legend('R = 0.005','R = 10')
ylabel('y dot')

subplot(513)
plot(t,x_1(:,3),t,x_2(:,3))
legend('R = 0.005','R = 10')
ylabel('theta')

subplot(514)
plot(t,x_1(:,4),t,x_2(:,4))
legend('R = 0.005','R = 10')
ylabel('theta dot')

subplot(515)
plot(t,u_1,t,u_2)
legend('R = 0.005','R = 10')
ylabel('u')
xlabel('time')

%% Part 4: Nonlinear Comparison
xdot_evaluated = subs(xdot,u,K1*x);
xdot_evaluated = eval(xdot_evaluated);
symvar(xdot_evaluated)
%%
% define controlled nonlinear system
syms t;
controlled_sys = matlabFunction(xdot_evaluated,'Vars',{t,x});

% set initial conditions
x0 = [-1;0;pi/4;0];
x0_2 = [-11;0;pi/4;0];
%x0 = [-12;0;pi/4;0];

% setting (ODE time range)
Tspan = linspace(0,10,1e3);

[t,X]=ode45(controlled_sys,Tspan,x0,options); % for nonlinear system, y0 = -1
u = K1*X.';
u = u.';

[t,X_2]=ode45(controlled_sys,Tspan,x0_2,options); % for nonlinear system, y0 = -11
u_2 = K1*X_2.';
u_2 = u_2.';

[t,x_1]=ode45(sys_1,Tspan,x0,options); % for linearized system, y0 = -1
u_1 = K1*x_1.';
u_1 = u_1.';

figure('Name','Nonlinear Comparison')
subplot(511)
plot(t,x_1(:,1),t,X(:,1))
legend('linearized','nonlinear')
ylabel('y')

subplot(512)
plot(t,x_1(:,2),t,X(:,2))
legend('linearized','nonlinear')
ylabel('y dot')

subplot(513)
plot(t,x_1(:,3),t,X(:,3))
legend('linearized','nonlinear')
ylabel('theta')

subplot(514)
plot(t,x_1(:,4),t,X(:,4))
legend('linearized','nonlinear')
ylabel('theta dot')

subplot(515)
plot(t,u_1,t,u)
legend('linearized','nonlinear')
ylabel('u')
xlabel('time')

% CHANGING Y0
figure('Name','Changing y0 Effect')
subplot(511)
plot(t,X(:,1),t,X_2(:,1))
legend('y0 = -1','y0 = -11')
ylabel('y')

subplot(512)
plot(t,X(:,2),t,X_2(:,2))
legend('y0 = -1','y0 = -11')
ylabel('y dot')

subplot(513)
plot(t,X(:,3),t,X_2(:,3))
legend('y0 = -1','y0 = -11')
ylabel('theta')

subplot(514)
plot(t,X(:,4),t,X_2(:,4))
legend('y0 = -1','y0 = -11')
ylabel('theta dot')

subplot(515)
plot(t,u,t,u_2)
legend('y0 = -1','y0 = -11')
ylabel('u')
xlabel('time')

