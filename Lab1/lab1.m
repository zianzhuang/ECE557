%% Lab 1 - Numerical Simulation of Dynamical Systems and Symbolic Linearization
%% Numerical Integration
close all;
clear all;
clc;

% define structure parameters
parameters.M = 0.2;  % Kg
parameters.l = 0.15; % m
parameters.g = 9.81; % m/sec^2

% set an initial condition
x0 = [0; sqrt(parameters.g/parameters.l)];
x0_2 = [0; 1.99*sqrt(parameters.g/parameters.l)];

% set ODE options
options = odeset('RelTol',1e-7,'AbsTol',1e-7);

% set time range
Tspan = linspace(0,10,1e3); % from 0 to 10

% numeric integration for x0_1
[t,x]=ode45(@pendulum,Tspan,x0,options,parameters);

% extract output
x1 = x(:,1); % theta
x2 = x(:,2); % thetadot

% plotting state vector vs. t for x0_1
figure('Name', 'Output 1: state vector vs. t (for x0_1)');
subplot(211)
plot(t,x1)
xlabel('Time','Fontsize',12)
ylabel('\theta','Fontsize',15)
subplot(212)
plot(t,x2)
xlabel('Time','Fontsize',12)
ylabel('$\dot{\theta}$','Interpreter','latex','Fontsize',15)

% plotting orbit for x0_1
figure('Name', 'Output 1: orbit (for x0_1)');
plot(x1,x2)
xlabel('\theta','Fontsize',15)
ylabel('$\dot{\theta}$','Interpreter','latex','Fontsize',15)

% numeric integration for x0_2
[t,x_2]=ode45(@pendulum,Tspan,x0_2,options,parameters);

% extract output
x1_2 = x_2(:,1); % theta
x2_2 = x_2(:,2); % thetadot

% plotting state vector vs. t for x0_2
figure('Name', 'Output 1: state vector vs. t (for x0_2)');
subplot(211)
plot(t,x1_2)
xlabel('Time','Fontsize',12)
ylabel('\theta','Fontsize',15)
subplot(212)
plot(t,x2_2)
xlabel('Time','Fontsize',12)
ylabel('$\dot{\theta}$','Interpreter','latex','Fontsize',15)

% plotting orbit for x0_2
figure('Name', 'Output 1: orbit (for x0_2)');
plot(x1_2,x2_2)
xlabel('\theta','Fontsize',15)
ylabel('$\dot{\theta}$','Interpreter','latex','Fontsize',15)
%% Symbolic Linearization
close all;
clc;

% declare symbolic variables
syms x1 x2 t u m l g real;

% define symbolic NL system
xdot1 = x2;
xdot2 = -g/l*sin(x1)-1/(m*l)*cos(x1)*u;
xdot = [xdot1;xdot2];
y = x1;
x = [x1;x2];

% compute symbolic jacobians
A = jacobian(xdot,x); % f/x; xdot = f(x,u)
B = jacobian(xdot,u); % f/u; xdot = f(x,u)
C = jacobian(y,x);    % h/x; y = h(x,u)
D = jacobian(y,u);    % h/u; y = h(x,u)

% define equilibrium states(pendulunm pointing downward)
xbar = [0;0];
ubar = 0;

% evaluate jacobians at equilibrium
A_evaluated = subs(A,x,xbar);
A_evaluated = subs(A_evaluated,u,ubar);
B_evaluated = subs(B,x,xbar);
B_evaluated = subs(B_evaluated,u,ubar);
C_evaluated = subs(C,x,xbar);
C_evaluated = subs(C_evaluated,u,ubar);
D_evaluated = subs(D,x,xbar);
D_evaluated = subs(D_evaluated,u,ubar);

% Output 2: display state space matrices at ([0;0],0)
A_evaluated % same as linearized LTI system
B_evaluated
C_evaluated
D_evaluated

% declare symbolic variables
syms theta real;

% define equilibrium states
xbar_new = [theta;0];
ubar_new = -m*g*tan(theta);

% evaluate jacobians at equilibrium
A_new = subs(A,x,xbar_new);
A_new = subs(A_new,u,ubar_new);
B_new = subs(B,x,xbar_new);
B_new = subs(B_new,u,ubar_new);

% Output 2: display state space matrices (A,B) at ([theta;0],-m*g*tan(theta))
A_new % same as Chapter 5 in textbook
B_new
%% Compare Symbolic Expression & Numerical Integration Solutions
close all;
clc;

% declare symbolic variables
syms z1 z2 real;
z = [z1;z2];

% define symbolic NL system at equilibrium state ([0;0],0)
zdot = A_evaluated*(z-xbar) + B_evaluated*(u-ubar);

% define symbolic augmented ODE
Xdot = [xdot;zdot];

% substitute parameters and input signal
Xdot_evaluated = subs(Xdot,m,parameters.M);
Xdot_evaluated = subs(Xdot_evaluated,g,parameters.g);
Xdot_evaluated = subs(Xdot_evaluated,l,parameters.l);
Xdot_evaluated = subs(Xdot_evaluated,u,0);

% check symbolic variables
symvar(Xdot_evaluated)

% create matlab function with augumented ODE
augmented_pend = matlabFunction(Xdot_evaluated,'Vars',{t,[x;z]});
%% numerical integration
% set an initial condition
x0 = [0; sqrt(parameters.g/parameters.l)];
z0 = x0;
X0 = [x0;z0];
x0_2 = [0; 1.99*sqrt(parameters.g/parameters.l)];
z0_2 = x0_2;
X0_2 = [x0_2;z0_2];

% numeric integration for x0_1
[t,X]=ode45(augmented_pend,Tspan,X0,options);

% extract output
X1 = X(:,1); % x1
X2 = X(:,2); % x2
X3 = X(:,3); % z1
X4 = X(:,4); % z2

% plotting state vector vs. t for x0_1
figure('Name', 'Output 3: state vector vs. t (for x0_1)');
subplot(211)
plot(t,X1)
hold on
plot(t,X3)
xlabel('Time','Fontsize',12)
ylabel('\theta','Fontsize',15)
legend('nonlinear','linearized')
subplot(212)
plot(t,X2)
hold on
plot(t,X4)
xlabel('Time','Fontsize',12)
ylabel('$\dot{\theta}$','Interpreter','latex','Fontsize',15)
legend('nonlinear','linearized')

% plotting orbit for x0_1
figure('Name', 'Output 3: orbit (for x0_1)');
plot(X1,X2)
hold on
plot(X3,X4)
xlabel('\theta','Fontsize',15)
ylabel('$\dot{\theta}$','Interpreter','latex','Fontsize',15)
legend('nonlinear','linearized')

% numeric integration for x0_2
[t,X_2]=ode45(augmented_pend,Tspan,X0_2,options);

% extract output
X1_2 = X_2(:,1); % x1_2
X2_2 = X_2(:,2); % x2_2
X3_2 = X_2(:,3); % z1_2
X4_2 = X_2(:,4); % z2_2

% plotting state vector vs. t for x0_2
figure('Name', 'Output 3: state vector vs. t (for x0_2)');
subplot(211)
plot(t,X1_2)
hold on
plot(t,X3_2)
xlabel('Time','Fontsize',12)
ylabel('\theta','Fontsize',15)
legend('nonlinear','linearized')
subplot(212)
plot(t,X2_2)
hold on
plot(t,X4_2)
xlabel('Time','Fontsize',12)
ylabel('$\dot{\theta}$','Interpreter','latex','Fontsize',15)
legend('nonlinear','linearized')

% plotting orbit for x0_2
figure('Name', 'Output 3: orbit (for x0_2)');
plot(X1_2,X2_2)
hold on
plot(X3_2,X4_2)
xlabel('\theta','Fontsize',15)
ylabel('$\dot{\theta}$','Interpreter','latex','Fontsize',15)
legend('nonlinear','linearized')
%% LTI Representations
close all;
clc;

% set numerical values for m,g,l
g = parameters.g;
l = parameters.l;
m = parameters.M;

% evaluate symbolic matrices
A_double = double(eval(A_evaluated));
B_double = double(eval(B_evaluated));
C_double = double(eval(C_evaluated));
D_double = double(eval(D_evaluated));

% define LTI object system
sys = ss(A_double,B_double,C_double,D_double);

% transfer function of linearized pendulum
TF = tf(sys) % same as the transfer function in textbook

% eigenvalues(V) and eigenvectors(D) of matrix A
[V,D]=eig(A_double)
%% Pendulum Stabilization
G = tf([-33.33],[1 0 65.40]);
C = tf(-30*[1 -10],[1 1000]); %-DC gain

zpk(1+C*G);

% find zeros
roots([1 1.011 55.56]);

[F,G,H,L]=ssdata(C);

Matrices.F = F;
Matrices.G = G;
Matrices.H = H;
Matrices.L = L;

z0 = 0; % controller intitial state
Z0 = [x0;z0];
[t,Z]=ode45(@controlled_pendulum,Tspan,Z0,options,parameters,Matrices);
% plotting state vector vs. t for x0_1
figure('Name', 'Output 5: state vector vs. t (for x0_1)');
subplot(311)
plot(t,Z(:,1))
xlabel('Time','Fontsize',12)
ylabel('\theta','Fontsize',15)
subplot(312)
plot(t,Z(:,2))
xlabel('Time','Fontsize',12)
ylabel('$\dot{\theta}$','Interpreter','latex','Fontsize',15)
subplot(313)
plot(t,Z(:,3))
xlabel('Time','Fontsize',12)
ylabel('\theta - controller','Fontsize',15)