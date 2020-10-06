function Zdot = pendulum(t,Z,parameters, Matrices)
    M = parameters.M;
    g = parameters.g;
    l = parameters.l;
    F = Matrices.F;
    G = Matrices.G;
    H = Matrices.H;
    L = Matrices.L;
    x1 = Z(1); 
    x2 = Z(2);
    z = Z(3);
    y = x1;
    u = H*z - L*y;
    xdot = [x2;(-g/l)*sin(x1)-1/(M*l)*cos(x1)*u]; % pendulum state
    zdot = F*z - G*y; % controller state
    Zdot = [xdot;zdot];
end