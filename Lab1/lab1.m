%% Lab 1 - Numerical Simulation of Dynamical Systems and Symbolic Linearization
%% OUTPUT1: Numerical Integration for NL ODES
close all;
clear all;
clc;

% define pendulum parameters
parameters.M = 0.2;
parameters.l = 0.15;
parameters.g = 9.81;

% define initial conditions for two cases
x0 = [0; sqrt(parameters.g/parameters.l)];
x0_2 = [0; 1.99*sqrt(parameters.g/parameters.l)]; % higher angular velocity

% setting (ODE options, time range)
options = odeset('RelTol',1e-7,'AbsTol',1e-7);
Tspan = linspace(0,10,1e3);

% SOLVE ODEs for nonlinear systems
[t,x]=ode45(@pendulum,Tspan,x0,options,parameters); % for x0_1

% extract output for plotting
x1 = x(:,1);
x2 = x(:,2);

% PLOT state vs. t for x0_1
figure('Name', 'Output 1: state vs. t (for x0_1)');
subplot(211);
plot(t,x1);
xlabel('Time','Fontsize',12);
ylabel('\theta','Fontsize',15);
title('Output 1: state vs. t (for x0_1)');
subplot(212);
plot(t,x2);
xlabel('Time','Fontsize',12);
ylabel('$\dot{\theta}$','Interpreter','latex','Fontsize',15);

% PLOT orbit for x0_1
figure('Name', 'Output 1: orbit (for x0_1)');
plot(x1,x2);
xlabel('\theta','Fontsize',15);
ylabel('$\dot{\theta}$','Interpreter','latex','Fontsize',15);
title('Output 1: orbit (for x0_1)');

% SOLVE ODEs for nonlinear systems
[t,x_2]=ode45(@pendulum,Tspan,x0_2,options,parameters); % for x0_2

% extract output for plotting
x1_2 = x_2(:,1);
x2_2 = x_2(:,2);

% PLOT state vs. t for x0_2
figure('Name', 'Output 1: state vs. t (for x0_2)');
subplot(211);
plot(t,x1_2);
xlabel('Time','Fontsize',12);
ylabel('\theta','Fontsize',15);
title('Output 1: state vs. t (for x0_2)');
subplot(212);
plot(t,x2_2);
xlabel('Time','Fontsize',12);
ylabel('$\dot{\theta}$','Interpreter','latex','Fontsize',15);

% PLOT orbit for x0_2
figure('Name', 'Output 1: orbit (for x0_2)');
plot(x1_2,x2_2);
xlabel('\theta','Fontsize',15);
ylabel('$\dot{\theta}$','Interpreter','latex','Fontsize',15);
title('Output 1: orbit (for x0_2)');

%% OUTPUT2: Symbolic Linearization for System at Given Equilibria
close all;
clc;

% CREATE symbolic experessions for NL system
% declare symbolic variables
syms x1 x2 t u m l g real;

% define nonlinear systems symbolically        % same as pendulum.m, but here we use symbolic approach
xdot = [x2; -(g/l)*sin(x1)-1/(m*l)*cos(x1)*u]; % f(x,u) = xdot
y = x1;                                        % h(x,u) = y
x = [x1; x2];

% LINEARIZE NL system at given equilibria
% compute symbolic jacobians
A = jacobian(xdot,x); % df/dx; f(x,u) = xdot
B = jacobian(xdot,u); % df/du; f(x,u) = xdot
C = jacobian(y,x);    % dh/dx; h(x,u) = y
D = jacobian(y,u);    % dh/du; h(x,u) = y

% define equilibrium states (pendulunm pointing downward) ([0; 0], 0)
xbar = [0; 0];
ubar = 0;

% evaluate jacobians at equilibria ([0; 0], 0)
A_evaluated = subs(A,x,xbar);
A_evaluated = subs(A_evaluated,u,ubar);
B_evaluated = subs(B,x,xbar);
B_evaluated = subs(B_evaluated,u,ubar);
C_evaluated = subs(C,x,xbar);
C_evaluated = subs(C_evaluated,u,ubar);
D_evaluated = subs(D,x,xbar);
D_evaluated = subs(D_evaluated,u,ubar);

% Output 2a: state space matrices for linearzied LTI system at ([0;0],0)
A_evaluated % verified with Chapter 3 in textbook
B_evaluated
C_evaluated
D_evaluated

% define a new equilibrium states ([theta; 0],-m*g*tan(theta))
syms theta real;
xbar_new = [theta; 0];
ubar_new = -m*g*tan(theta);

% evaluate jacobians at new equilibria ([theta; 0],-m*g*tan(theta))
A_new = subs(A,x,xbar_new);
A_new = subs(A_new,u,ubar_new);
B_new = subs(B,x,xbar_new);
B_new = subs(B_new,u,ubar_new);

% Output 2b: state space matrices for linearzied LTI system at ([theta; 0],-m*g*tan(theta))
A_new % verified with Chapter 5 in textbook
B_new

%% OUTPUT3-1: Symbolic Expression for Augmented ODES (NL; LTI)
close all;
clc;

% DEFINE augmented ODEs
% declare symbolic variables
syms z1 z2 real;
z = [z1;z2]; % state for linearized LTI system

% define linearized systems symbolically at equilibrium state ([0;0],0)
zdot = A_evaluated*(z-xbar) + B_evaluated*(u-ubar);

% define symbolic augmented ODE
Xdot = [xdot;zdot]; % IN FORM OF [NL; LTI]

% substitute parameters and input signal
Xdot_evaluated = subs(Xdot,m,parameters.M);
Xdot_evaluated = subs(Xdot_evaluated,g,parameters.g);
Xdot_evaluated = subs(Xdot_evaluated,l,parameters.l);
Xdot_evaluated = subs(Xdot_evaluated,u,0); % SET control input

% check symbolic variables
symvar(Xdot_evaluated) % only x and z variables are not specified

% create a matlab function with augmented ODE
augmented_pend = matlabFunction(Xdot_evaluated,'Vars',{t,[x;z]});

%% OUTPUT3-2: Numerical Integration for Augmented ODES (NL; LTI)
% define initial conditions for two cases
z0 = x0;
X0 = [x0;z0]; % close to equilibria ([0; 0], 0)
z0_2 = x0_2;
X0_2 = [x0_2;z0_2]; % away from equilibria ([0; 0], 0)

% SOLVE augmented ODEs
[t,X]=ode45(augmented_pend,Tspan,X0,options); % for x0_1

% extract output for plotting
X1 = X(:,1); % x1
X2 = X(:,2); % x2
Z1 = X(:,3); % z1
Z2 = X(:,4); % z2

% PLOT state vs. t for x0_1
figure('Name', 'Output 3: state vs. t (for x0_1)');
subplot(211);
plot(t,X1);
hold on;
plot(t,Z1);
xlabel('Time','Fontsize',12);
ylabel('\theta','Fontsize',15);
legend('nonlinear','linearized');
title('Output 3: state vs. t (for x0_1)');
subplot(212);
plot(t,X2);
hold on;
plot(t,Z2);
xlabel('Time','Fontsize',12);
ylabel('$\dot{\theta}$','Interpreter','latex','Fontsize',15);
legend('nonlinear','linearized');

% PLOT orbit for x0_1
figure('Name', 'Output 3: orbit (for x0_1)');
plot(X1,X2);
hold on;
plot(Z1,Z2);
xlabel('\theta','Fontsize',15);
ylabel('$\dot{\theta}$','Interpreter','latex','Fontsize',15);
legend('nonlinear','linearized');
title('Output 3: orbit (for x0_1)');

% SOLVE augmented ODEs
[t,X_2]=ode45(augmented_pend,Tspan,X0_2,options); % for x0_2

% extract output for plotting
X1_2 = X_2(:,1); % x1_2
X2_2 = X_2(:,2); % x2_2
Z1_2 = X_2(:,3); % z1_2
Z2_2 = X_2(:,4); % z2_2

% PLOT state vs. t for x0_2
figure('Name', 'Output 3: state vector vs. t (for x0_2)');
subplot(211);
plot(t,X1_2);
hold on;
plot(t,Z1_2);
xlabel('Time','Fontsize',12);
ylabel('\theta','Fontsize',15);
legend('nonlinear','linearized');
title('Output 3: state vs. t (for x0_2)');
subplot(212);
plot(t,X2_2);
hold on;
plot(t,Z2_2);
xlabel('Time','Fontsize',12);
ylabel('$\dot{\theta}$','Interpreter','latex','Fontsize',15);
legend('nonlinear','linearized');

% PLOT orbit for x0_2
figure('Name', 'Output 3: orbit (for x0_2)');
plot(X1_2,X2_2);
hold on;
plot(Z1_2,Z2_2);
xlabel('\theta','Fontsize',15);
ylabel('$\dot{\theta}$','Interpreter','latex','Fontsize',15);
legend('nonlinear','linearized');
title('Output 3: orbit (for x0_2)');

%% OUTPUT4: LTI Representations
close all;
clc;

% set numerical values for m,g,l
g = parameters.g;
l = parameters.l;
m = parameters.M;

% evaluate symbolic matrices and convert them to double
A_double = double(eval(A_evaluated));
B_double = double(eval(B_evaluated));
C_double = double(eval(C_evaluated));
D_double = double(eval(D_evaluated));

% define LTI object system
sys = ss(A_double,B_double,C_double,D_double);

% transfer function of linearized pendulum
TF = tf(sys) % verified with Chapter 3 in textbook

% eigenvalues(D) and eigenvectors(V) of matrix A
[V,D] = eig(A_double)

% poles of transfer function
P = pole(sys) % same as D

%% OUTPUT5:Pendulum Stabilization
close all;
clc;

G = tf([-33.33],[1 0 65.40]); % transfer function of linearized system
C = tf([14.513],[1 10.88]) % transfer function of designed controller

% find poles of closed-loop system transfer function
zpk(1+C*G)
ZERO1 = roots([1 6.249]) % lie in C-
ZERO2 = roots([1 4.631 36.46]) % lie in C-

% convert controller to state space form and extract state matrices
[F,G,H,L]=ssdata(C);

% construct structural parameters
Matrices.F = F;
Matrices.G = G;
Matrices.H = H;
Matrices.L = L;

% TRY different initial conditions
x0 = [0; sqrt(parameters.g/parameters.l)]; % initial condition worked
% x0 = [3; sqrt(parameters.g/parameters.l)];
% x0 = [3.1; sqrt(parameters.g/parameters.l)]; % initial condition failed

% define controller initial state
z0 = 0; 
Z0 = [x0;z0];

% SOLVE ODEs for closed-loop system
[t,Z]=ode45(@controlled_pendulum,Tspan,Z0,options,parameters,Matrices);

% PLOT state vs. t for x0
figure('Name', 'Output 5: state vs. t (for x0)');
subplot(311)
plot(t,Z(:,1))
xlabel('Time','Fontsize',12)
ylabel('\theta','Fontsize',15)
title('Output 5: state vs. t (for x0)');
subplot(312)
plot(t,Z(:,2))
xlabel('Time','Fontsize',12)
ylabel('$\dot{\theta}$','Interpreter','latex','Fontsize',15)
subplot(313)
plot(t,Z(:,3))
xlabel('Time','Fontsize',12)
ylabel('\theta - controller','Fontsize',15)