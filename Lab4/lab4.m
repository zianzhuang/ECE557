%% Lab 4 - Output Feedback Stabilization
%% Symbolic Linearization
close all;
clear all;
clc;

% declare symbolic variables
syms z z_dot theta theta_dot t u M m l g real;

% define nonlinear systems symbolically
x = [z;z_dot;theta;theta_dot];
y = [z;theta];
xdot1 = z_dot;
xdot2 = (-m*l*sin(theta)*theta_dot^2+m*g*sin(theta)*cos(theta)+u)/(M+m*sin(theta)^2);
xdot3 = theta_dot;
xdot4 = (-m*l*sin(theta)*cos(theta)*theta_dot^2+(M+m)*g*sin(theta)+u*cos(theta))/l/(M+m*sin(theta)^2);
xdot = [xdot1;xdot2;xdot3;xdot4];

% define structure parameters
parameters.M = 1.0731;
parameters.m = 0.23;
parameters.l = 0.3302;
parameters.g = 9.8;

% compute symbolic jacobians
A = jacobian(xdot,x);
B = jacobian(xdot,u);
C = jacobian(y,x);

% define equilibrium states
xbar = [0;0;0;0];

% evaluate jacobians at equilibria
A_evaluated = subs(A,x,xbar);
B_evaluated = subs(B,x,xbar);
C_evaluated = subs(C,x,xbar);

% set numerical values for m,M,g,l
g = parameters.g;
l = parameters.l;
M = parameters.M;
m = parameters.m;

% evaluate symbolic matrices and convert them to double
A_double = double(eval(A_evaluated));
B_double = double(eval(B_evaluated));
C_double = double(eval(C_evaluated));

%% Part 1-1: Observability and Simulation Setup
% verify controllability
Q = obsv(A_double, C_double)
Q_rank = rank(Q)

% design observer gains
pole_obs_1 = [-10 -11 -12 -13];
K_obs_1 = -place(A_double.', C_double.', pole_obs_1); % convention u = Kx
L1 = K_obs_1.'

pole_obs_2 = [-40 -41 -42 -43];
K_obs_2 = -place(A_double.', C_double.', pole_obs_2); % convention u = Kx
L2 = K_obs_2.'

% design full state feedback controller
pole_original = [-5.1 -5.2 -5.3 -5.4];
K_original = -place(A_double, B_double, pole_original); % convention u = Kx

% construct structural parameters
Matrices.A = A_double;
Matrices.B = B_double;
Matrices.C = C_double;
Gain.K = K_original;

% define initial conditions
x0 = [-0.5; 0; -pi/4; 0];
x_hat0 = [0; 0; 0; 0];
Z0 = [x0; x_hat0];
      
% setting (ODE options, time range)
options = odeset('RelTol',1e-7,'AbsTol',1e-7);
Tspan = linspace(0,10,1e3);

%% Part 1-2: Noiseless State Estimation for Linearized System 
% Intergrate for L1
Gain.L = L1;
[t,Z_1]=ode45(@linearobs,Tspan,Z0,options,Matrices,Gain);
LSL1 = Z_1; % Linear State Feedback Controller for L1

% Intergrate for L2
Gain.L = L2;
[t,Z_2]=ode45(@linearobs,Tspan,Z0,options,Matrices,Gain);
LSL2 = Z_2; % Linear State Feedback Controller for L2

figure('Name','State Estimation Error for Noiseless Linearized System')
subplot(411)
plot(t,Z_1(:,5)-Z_1(:,1),t,Z_2(:,5)-Z_2(:,1))
legend('L1','L2')
ylabel('Error in x~ 1')

subplot(412)
plot(t,Z_1(:,6)-Z_1(:,2),t,Z_2(:,6)-Z_2(:,2))
legend('L1','L2')
ylabel('Error in x~ 2')

subplot(413)
plot(t,Z_1(:,7)-Z_1(:,3),t,Z_2(:,7)-Z_2(:,3))
legend('L1','L2')
ylabel('Error in x~ 3')

subplot(414)
plot(t,Z_1(:,8)-Z_1(:,4),t,Z_2(:,8)-Z_2(:,4))
legend('L1','L2')
ylabel('Error in x~ 4')
xlabel('Time')

%% Part 1-3: Noiseless State Estimation for Nonlinear System
% Intergrate for L1
Gain.L = L1;
[t,Z_1]=ode45(@nonlinearobs,Tspan,Z0,options,Matrices,Gain,parameters);
NLSL1 = Z_1; % Nonlinear State Feedback Controller for L1

% Intergrate for L2
Gain.L = L2;
[t,Z_2]=ode45(@nonlinearobs,Tspan,Z0,options,Matrices,Gain,parameters);
NLSL2 = Z_2; % Nonlinear State Feedback Controller for L2

figure('Name','State Estimation Error for Noiseless Nonlinear System')
subplot(411)
plot(t,Z_1(:,5)-Z_1(:,1),t,Z_2(:,5)-Z_2(:,1))
legend('L1','L2')
ylabel('Error in x~ 1')

subplot(412)
plot(t,Z_1(:,6)-Z_1(:,2),t,Z_2(:,6)-Z_2(:,2))
legend('L1','L2')
ylabel('Error in x~ 2')

subplot(413)
plot(t,Z_1(:,7)-Z_1(:,3),t,Z_2(:,7)-Z_2(:,3))
legend('L1','L2')
ylabel('Error in x~ 3')

subplot(414)
plot(t,Z_1(:,8)-Z_1(:,4),t,Z_2(:,8)-Z_2(:,4))
legend('L1','L2')
ylabel('Error in x~ 4')
xlabel('Time')

%% Part 2-1: Generate Random Noise
cov = [0.005 0; 0 0.001];
L_chol = chol(cov,'lower'); % cholesky decomposition
noise_x = L_chol*randn(2,1000);

%% Part 2-2: State Estimation for Linearized System with Noise
% Intergrate for L1
Gain.L = L1;
[t,Z_1]=ode45(@linearobsnoise,Tspan,Z0,options,Matrices,Gain,noise_x,Tspan);

% Intergrate for L2
Gain.L = L2;
[t,Z_2]=ode45(@linearobsnoise,Tspan,Z0,options,Matrices,Gain,noise_x,Tspan);

figure('Name','State Estimation Error for Linearized System with Noise')
subplot(411)
plot(t,Z_1(:,5)-Z_1(:,1),t,Z_2(:,5)-Z_2(:,1))
legend('L1','L2')
ylabel('Error in x~ 1')

subplot(412)
plot(t,Z_1(:,6)-Z_1(:,2),t,Z_2(:,6)-Z_2(:,2))
legend('L1','L2')
ylabel('Error in x~ 2')

subplot(413)
plot(t,Z_1(:,7)-Z_1(:,3),t,Z_2(:,7)-Z_2(:,3))
legend('L1','L2')
ylabel('Error in x~ 3')

subplot(414)
plot(t,Z_1(:,8)-Z_1(:,4),t,Z_2(:,8)-Z_2(:,4))
legend('L1','L2')
ylabel('Error in x~ 4')
xlabel('Time')

% compute mean squared estimation error
MSE_L1 = immse(Z_1(501:1000,1:4), Z_1(501:1000,5:8)) % for L1
MSE_L2 = immse(Z_2(501:1000,1:4), Z_2(501:1000,5:8)) % for L2

%% Part 3-1: Noiseless Output Feedback Control for Linearized System 
% Intergrate for L1
Gain.L = L1;
[t,Z_1]=ode45(@linearcontrol,Tspan,Z0,options,Matrices,Gain);

% Intergrate for L2
Gain.L = L2;
[t,Z_2]=ode45(@linearcontrol,Tspan,Z0,options,Matrices,Gain);

figure('Name','State Feedback vs. Output Feedback Control for Noiseless Linearized System')
subplot(411)
plot(t,Z_1(:,1),t,Z_2(:,1),t,LSL1(:,1),t,LSL2(:,1))
legend('OutputFB L1','OutputFB L2','StateFB L1','StateFB L2')
ylabel('z')

subplot(412)
plot(t,Z_1(:,2),t,Z_2(:,2),t,LSL1(:,2),t,LSL2(:,2))
legend('OutputFB L1','OutputFB L2','StateFB L1','StateFB L2')
ylabel('z dot')

subplot(413)
plot(t,Z_1(:,3),t,Z_2(:,3),t,LSL1(:,3),t,LSL2(:,3))
legend('OutputFB L1','OutputFB L2','StateFB L1','StateFB L2')
ylabel('theta')

subplot(414)
plot(t,Z_1(:,4),t,Z_2(:,4),t,LSL1(:,4),t,LSL2(:,4))
legend('OutputFB L1','OutputFB L2','StateFB L1','StateFB L2')
ylabel('theta dot')
xlabel('Time')

%% Part 3-2: Noiseless Output Feedback Control for Nonlinear System 
% Intergrate for L1
Gain.L = L1;
[t,Z_1]=ode45(@nonlinearcontrol,Tspan,Z0,options,Matrices,Gain,parameters);

% Intergrate for L2
Gain.L = L2;
[t,Z_2]=ode45(@nonlinearcontrol,Tspan,Z0,options,Matrices,Gain,parameters);

figure('Name','State Feedback vs. Output Feedback Control for Noiseless Nonlinear System')
subplot(411)
plot(t,Z_1(:,1),t,Z_2(:,1),t,NLSL1(:,1),t,NLSL2(:,1))
legend('OutputFB L1','OutputFB L2','StateFB L1','StateFB L2')
ylabel('z')

subplot(412)
plot(t,Z_1(:,2),t,Z_2(:,2),t,NLSL1(:,2),t,NLSL2(:,2))
legend('OutputFB L1','OutputFB L2','StateFB L1','StateFB L2')
ylabel('z dot')

subplot(413)
plot(t,Z_1(:,3),t,Z_2(:,3),t,NLSL1(:,3),t,NLSL2(:,3))
legend('OutputFB L1','OutputFB L2','StateFB L1','StateFB L2')
ylabel('theta')

subplot(414)
plot(t,Z_1(:,4),t,Z_2(:,4),t,NLSL1(:,4),t,NLSL2(:,4))
legend('OutputFB L1','OutputFB L2','StateFB L1','StateFB L2')
ylabel('theta dot')
xlabel('Time')

%% Part 4: Output Feedback Control for Nonlinear System with Noise
% Intergrate for L1
Gain.L = L1;
[t,NZ_1]=ode45(@NLnoisycontrol,Tspan,Z0,options,Matrices,Gain,parameters,noise_x,Tspan);

% Intergrate for L2
Gain.L = L2;
[t,NZ_2]=ode45(@NLnoisycontrol,Tspan,Z0,options,Matrices,Gain,parameters,noise_x,Tspan);

figure('Name','Noise Impact on Output Feedback Control for Nonlinear System')
subplot(411)
plot(t,NZ_1(:,1),t,NZ_2(:,1),t,Z_1(:,1),t,Z_2(:,1))
legend('Noisy L1','Noisy L2','Noiseless L1','Noiseless L2')
ylabel('z')

subplot(412)
plot(t,NZ_1(:,2),t,NZ_2(:,2),t,Z_1(:,2),t,Z_2(:,2))
legend('Noisy L1','Noisy L2','Noiseless L1','Noiseless L2')
ylabel('z dot')

subplot(413)
plot(t,NZ_1(:,3),t,NZ_2(:,3),t,Z_1(:,3),t,Z_2(:,3))
legend('Noisy L1','Noisy L2','Noiseless L1','Noiseless L2')
ylabel('theta')

subplot(414)
plot(t,NZ_1(:,4),t,NZ_2(:,4),t,Z_1(:,4),t,Z_2(:,4))
legend('Noisy L1','Noisy L2','Noiseless L1','Noiseless L2')
ylabel('theta dot')
xlabel('Time')