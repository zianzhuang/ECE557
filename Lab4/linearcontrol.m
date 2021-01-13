function Zdot = linearcontrol(t,Z,Matrices,Gain)
% extract structral parameters
A = Matrices.A;
B = Matrices.B;
C = Matrices.C;
K = Gain.K;
L = Gain.L;

% construct ODEs
x = Z(1:4,:);
x_hat = Z(5:8,:);

x_dot = A*x + B*K*x_hat;
y = C*x;
x_hat_dot = (A+L*C+B*K)*x_hat - L*y;

Zdot = [x_dot(1);x_dot(2);x_dot(3);x_dot(4);x_hat_dot(1);x_hat_dot(2);x_hat_dot(3);x_hat_dot(4)];
end

