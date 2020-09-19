close all;
clear all; 
% create variable: parameters
parameters.M = 0.2;
parameters.l = 0.15;
parameters.g = 9.81;

% choose an initial condition
% x0 = [0; sqrt(parameters.g/parameters.l)];
x0 = [0; 1.99*sqrt(parameters.g/parameters.l)];

% set relative and absolute tolerances
options = odeset('RelTol',1e-7,'AbsTol',1e-7);

% set time span
Tspan = linspace(0,10,1e3);

% integration
[t,x]=ode45(@pendulum,Tspan,x0,options,parameters);

% extract output
x1 = x(:,1);
x2 = x(:,2);

% plot
figure;
subplot(211)
plot(t,x1)
xlabel('Time','Fontsize',15)
ylabel('\theta','Fontsize',15)
subplot(212)
plot(t,x2)
xlabel('Time','Fontsize',15)
ylabel('$\dot{\theta}$','Interpreter','latex','Fontsize',15)

figure;
plot(x1,x2)
xlabel('\theta','Fontsize',15)
ylabel('$\dot{\theta}$','Interpreter','latex','Fontsize',15)