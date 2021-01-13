function Zdot = NLnoisycontrol(t,Z,Matrices,Gain,parameters,noise_x,Tspan)
% extract structral parameters
A = Matrices.A;
B = Matrices.B;
C = Matrices.C;
K = Gain.K;
L = Gain.L;
g = parameters.g;
l = parameters.l;
M = parameters.M;
m = parameters.m;

% anonoymous function for noise interpolation
W = @(t) interp1(Tspan, noise_x.', t);

% construct ODEs
x = Z(1:4,:);
x_hat = Z(5:8,:);

u = K*x_hat;
y = C*x + W(t).';
x_hat_dot = (A+L*C+B*K)*x_hat - L*y;

z = Z(1,:);
z_dot = Z(2,:);
theta = Z(3,:);
theta_dot = Z(4,:);
xdot1 = z_dot;
xdot2 = (-m*l*sin(theta)*theta_dot^2+m*g*sin(theta)*cos(theta)+u)/(M+m*sin(theta)^2);
xdot3 = theta_dot;
xdot4 = (-m*l*sin(theta)*cos(theta)*theta_dot^2+(M+m)*g*sin(theta)+u*cos(theta))/l/(M+m*sin(theta)^2);
x_dot = [xdot1;xdot2;xdot3;xdot4];

Zdot = [x_dot(1);x_dot(2);x_dot(3);x_dot(4);x_hat_dot(1);x_hat_dot(2);x_hat_dot(3);x_hat_dot(4)];
end

